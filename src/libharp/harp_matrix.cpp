// @COPYRIGHT@

#include <harp_internal.hpp>


using namespace std;
using namespace harp;


harp::sparse_block::sparse_block ( matrix_sparse const & orig ) {

  global_rows = orig.Height();

  local_firstrow = orig.FirstLocalRow();
  local_rows = orig.LocalHeight();
  local_vals = orig.NumLocalEntries();

  local_row.resize ( local_vals );
  local_col.resize ( local_vals );
  data.resize ( local_vals );

  local_row.reserve ( local_vals );
  local_col.reserve ( local_vals );
  data.reserve ( local_vals );

  local_row_offset.resize ( local_rows );
  local_row_offset.reserve ( local_rows );

  local_row_nnz.resize ( local_rows );
  local_row_nnz.reserve ( local_rows );

  for ( int i = 0; i < local_rows; ++i ) {
    local_row_offset[i] = orig.LocalEntryOffset ( i );
    local_row_nnz[i] = orig.NumConnections ( i );
  }

  for ( int i = 0; i < local_vals; ++i ) {
    local_row[i] = orig.Row ( i );
    local_col[i] = orig.Col ( i );
    data[i] = orig.Value ( i );
  }

}


harp::sparse_block::sparse_block ( char * packed, size_t nbytes ) {

  char * cur = packed;

  global_rows = *( (int*)cur );
  cur += sizeof(int);

  local_firstrow = *( (int*)cur );
  cur += sizeof(int);

  local_rows = *( (int*)cur );
  cur += sizeof(int);

  local_vals = *( (int*)cur );
  cur += sizeof(int);

  local_row.resize ( local_vals );
  local_col.resize ( local_vals );
  data.resize ( local_vals );

  local_row.reserve ( local_vals );
  local_col.reserve ( local_vals );
  data.reserve ( local_vals );

  local_row_offset.resize ( local_rows );
  local_row_offset.reserve ( local_rows );

  local_row_nnz.resize ( local_rows );
  local_row_nnz.reserve ( local_rows );

  memcpy ( (void*)(&(local_row[0])), (void*)cur, (size_t)local_vals * sizeof(int) );
  cur += (size_t)local_vals * sizeof(int);

  memcpy ( (void*)(&(local_col[0])), (void*)cur, (size_t)local_vals * sizeof(int) );
  cur += (size_t)local_vals * sizeof(int);

  memcpy ( (void*)(&(local_row_offset[0])), (void*)cur, (size_t)local_rows * sizeof(int) );
  cur += (size_t)local_rows * sizeof(int);

  memcpy ( (void*)(&(local_row_nnz[0])), (void*)cur, (size_t)local_rows * sizeof(int) );
  cur += (size_t)local_rows * sizeof(int);

  memcpy ( (void*)(&(data[0])), (void*)cur, (size_t)local_vals * sizeof(double) );

}


char * harp::sparse_block::pack ( size_t & nbytes ) {
  char * buf;

  nbytes = 4 * sizeof(int) + (size_t)local_rows * ( 2 * sizeof(int) ) + (size_t)local_vals * ( sizeof(double) + 2 * sizeof(int) );

  buf = (char*) malloc ( nbytes );
  if ( ! buf ) {
    HARP_THROW( "cannot allocate sparse_block buffer" );
  }

  char * cur = buf;

  *( (int*)cur ) = global_rows;
  cur += sizeof(int);

  *( (int*)cur ) = local_firstrow;
  cur += sizeof(int);

  *( (int*)cur ) = local_rows;
  cur += sizeof(int);

  *( (int*)cur ) = local_vals;
  cur += sizeof(int);

  memcpy ( (void*)cur, (void*)(&(local_row[0])), (size_t)local_vals * sizeof(int) );
  cur += (size_t)local_vals * sizeof(int);

  memcpy ( (void*)cur, (void*)(&(local_col[0])), (size_t)local_vals * sizeof(int) );
  cur += (size_t)local_vals * sizeof(int);

  memcpy ( (void*)cur, (void*)(&(local_row_offset[0])), (size_t)local_rows * sizeof(int) );
  cur += (size_t)local_rows * sizeof(int);

  memcpy ( (void*)cur, (void*)(&(local_row_nnz[0])), (size_t)local_rows * sizeof(int) );
  cur += (size_t)local_rows * sizeof(int);

  memcpy ( (void*)cur, (void*)(&(data[0])), (size_t)local_vals * sizeof(double) );

  return buf;
}


void harp::dist_matrix_zero ( matrix_dist & mat ) {

  size_t hlocal = mat.LocalHeight();
  size_t wlocal = mat.LocalWidth();

  for ( size_t i = 0; i < wlocal; ++i ) {
    for ( size_t j = 0; j < hlocal; ++j ) {
      mat.SetLocal ( j, i, 0.0 );
    }
  }

  return;
}


void harp::local_matrix_zero ( matrix_local & mat ) {

  size_t hlocal = mat.Height();
  size_t wlocal = mat.Width();

  for ( size_t i = 0; i < wlocal; ++i ) {
    for ( size_t j = 0; j < hlocal; ++j ) {
      mat.Set ( j, i, 0.0 );
    }
  }

  return;
}


void harp::eigen_decompose ( matrix_dist & invcov, matrix_dist & D, matrix_dist & W ) {

  D.ResizeTo ( invcov.Height(), 1 );
  W.ResizeTo ( invcov.Height(), invcov.Height() );

  matrix_dist temp ( invcov );

  elem::DistMatrix < double, elem::VR, elem::STAR > eigvals;

  elem::HermitianEig ( elem::LOWER, temp, eigvals, W );

  elem::SortEig( eigvals, W );

  D = eigvals;

  return;
}


void harp::eigen_compose ( eigen_op op, matrix_dist & D, matrix_dist & W, matrix_dist & out ) {

  double threshold = 1.0e-15;

  out.ResizeTo ( W.Height(), W.Height() );

  // Get full local copy of eigenvalues

  matrix_local scaled ( D.Height(), 1 );
  local_matrix_zero ( scaled );

  elem::AxpyInterface < double > globloc;
  globloc.Attach( elem::GLOBAL_TO_LOCAL, D );
  globloc.Axpy ( 1.0, scaled, 0, 0 );
  globloc.Detach();

  // scale eigenvalues and apply threshold

  double max = 0.0;
  double min;
  double val;

  for ( size_t i = 0; i < scaled.Height(); ++i ) {
    val = scaled.Get( i, 0 );
    if ( val > max ) {
      max = val;
    }
  }

  min = max * threshold;

  switch ( op ) {
    case EIG_SQRT:
      for ( size_t i = 0; i < scaled.Height(); ++i ) {
        val = scaled.Get( i, 0 );
        val = ( val < min ) ? sqrt( min ) : sqrt( val );
        scaled.Set( i, 0, val );
      }
      break;
    case EIG_INVSQRT:
      for ( size_t i = 0; i < scaled.Height(); ++i ) {
        val = scaled.Get( i, 0 );
        val = ( val < min ) ? ( 1.0 / sqrt(min) ) : ( 1.0 / sqrt(val) );
        scaled.Set( i, 0, val );
      }
      break;
    default:
      for ( size_t i = 0; i < scaled.Height(); ++i ) {
        val = scaled.Get( i, 0 );
        val = ( val < min ) ? min : val;
        scaled.Set( i, 0, val );
      }
      break;
  }

  // Compute temp = W * op(D), by modifying our local elements. 

  matrix_dist temp ( W );

  size_t hlocal = temp.LocalHeight();
  size_t wlocal = temp.LocalWidth();

  size_t coloff = temp.RowShift();
  size_t colstride = temp.RowStride();
  size_t col;

  double mval;

  for ( size_t i = 0; i < wlocal; ++i ) {
    for ( size_t j = 0; j < hlocal; ++j ) {
      col = coloff + i * colstride;
      mval = temp.GetLocal ( j, i );
      mval *= scaled.Get ( col, 0 );
      temp.SetLocal ( j, i, mval );
    }
  }

  // compute out = ( W * op(D) ) * W^T

  elem::Gemm ( elem::NORMAL, elem::TRANSPOSE, 1.0, temp, W, 0.0, out );

  return;
}


void harp::column_norm ( matrix_dist & mat, matrix_dist & S ) {

  S.ResizeTo ( mat.Height(), 1 );
  dist_matrix_zero ( S );

  // Sum local columns

  matrix_local Sloc ( S.Height(), 1 );
  local_matrix_zero ( Sloc );

  size_t hlocal = mat.LocalHeight();
  size_t wlocal = mat.LocalWidth();

  size_t rowoff = mat.ColShift();
  size_t rowstride = mat.ColStride();
  size_t row;

  double mval;

  for ( size_t i = 0; i < wlocal; ++i ) {
    for ( size_t j = 0; j < hlocal; ++j ) {
      row = rowoff + j * rowstride;
      mval = mat.GetLocal ( j, i ) + Sloc.Get ( row, 0 );
      Sloc.Set ( row, 0, mval );
    }
  }

  // Reduce

  elem::AxpyInterface < double > locglob;
  locglob.Attach( elem::LOCAL_TO_GLOBAL, S );
  locglob.Axpy ( 1.0, Sloc, 0, 0 );
  locglob.Detach();

  // Invert

  hlocal = S.LocalHeight();
  wlocal = S.LocalWidth();

  for ( size_t i = 0; i < wlocal; ++i ) {
    for ( size_t j = 0; j < hlocal; ++j ) {
      mval = S.GetLocal ( j, i );
      mval = 1.0 / mval;
      S.SetLocal ( j, i, mval );
    }
  }

  return;
}


void harp::apply_norm ( matrix_dist & S, matrix_dist & mat ) {

  // Get local copy of S

  matrix_local Sloc ( S.Height(), 1 );
  local_matrix_zero ( Sloc );

  elem::AxpyInterface < double > globloc;
  globloc.Attach( elem::GLOBAL_TO_LOCAL, S );
  globloc.Axpy ( 1.0, Sloc, 0, 0 );
  globloc.Detach();

  // Apply to local matrix

  size_t hlocal = mat.LocalHeight();
  size_t wlocal = mat.LocalWidth();

  size_t rowoff = mat.ColShift();
  size_t rowstride = mat.ColStride();
  size_t row;

  double mval;

  for ( size_t i = 0; i < wlocal; ++i ) {
    for ( size_t j = 0; j < hlocal; ++j ) {
      row = rowoff + j * rowstride;
      mval = mat.GetLocal ( j, i );
      mval *= Sloc.Get ( row, 0 );
      mat.SetLocal ( j, i, mval );
    }
  }

  return;
}


void harp::apply_inverse_norm ( matrix_dist & S, matrix_dist & mat ) {

  // Get local copy of S

  matrix_local Sloc ( S.Height(), 1 );
  local_matrix_zero ( Sloc );

  elem::AxpyInterface < double > globloc;
  globloc.Attach( elem::GLOBAL_TO_LOCAL, S );
  globloc.Axpy ( 1.0, Sloc, 0, 0 );
  globloc.Detach();

  // Apply to local matrix

  size_t hlocal = mat.LocalHeight();
  size_t wlocal = mat.LocalWidth();

  size_t rowoff = mat.ColShift();
  size_t rowstride = mat.ColStride();
  size_t row;

  double mval;

  for ( size_t i = 0; i < wlocal; ++i ) {
    for ( size_t j = 0; j < hlocal; ++j ) {
      row = rowoff + j * rowstride;
      mval = mat.GetLocal ( j, i );
      mval /= Sloc.Get ( row, 0 );
      mat.SetLocal ( j, i, mval );
    }
  }

  return;
}


void harp::norm ( matrix_dist & D, matrix_dist & W, matrix_dist & S ) {

  matrix_dist temp ( W );
  dist_matrix_zero ( temp );  

  eigen_compose ( EIG_SQRT, D, W, temp );

  column_norm ( temp, S );
  
  return;
}


