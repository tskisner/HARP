// @COPYRIGHT@

#include <harp_mpi_math_internal.hpp>


using namespace std;
using namespace harp;



harp::mpi_matrix_sparse::mpi_matrix_sparse ( boost::mpi::communicator const & comm, size_t rows, size_t cols ) {
  comm_ = comm;
  rows_ = rows;
  cols_ = cols;

  // compute row offsets of our block

  mpi_dist_1D ( comm, rows, block_.rows, block_.firstrow );

}


void harp::mpi_matrix_zero ( mpi_matrix & mat ) {

  size_t hlocal = mat.LocalHeight();
  size_t wlocal = mat.LocalWidth();

  for ( size_t i = 0; i < wlocal; ++i ) {
    for ( size_t j = 0; j < hlocal; ++j ) {
      mat.SetLocal ( j, i, 0.0 );
    }
  }

  return;
}


void harp::local_matrix_zero ( elem_matrix_local & mat ) {

  size_t hlocal = mat.Height();
  size_t wlocal = mat.Width();

  for ( size_t i = 0; i < wlocal; ++i ) {
    for ( size_t j = 0; j < hlocal; ++j ) {
      mat.Set ( j, i, 0.0 );
    }
  }

  return;
}


void harp::mpi_eigen_decompose ( mpi_matrix const & invcov, mpi_matrix & D, mpi_matrix & W ) {

  D.ResizeTo ( invcov.Height(), 1 );
  W.ResizeTo ( invcov.Height(), invcov.Height() );

  mpi_matrix temp ( invcov );

  elem::DistMatrix < double, elem::VR, elem::STAR > eigvals ( invcov.Grid() );

  elem::HermitianEig ( elem::LOWER, temp, eigvals, W );

  // I don't think we need this, since we never need to order the
  // eigenpairs.
  //elem::SortEig( eigvals, W );

  D = eigvals;

  return;
}


void harp::mpi_eigen_compose ( eigen_op op, mpi_matrix const & D, mpi_matrix const & W, mpi_matrix & out ) {

  double threshold = 1.0e-12;

  out.ResizeTo ( W.Height(), W.Height() );

  // Get full local copy of eigenvalues

  elem_matrix_local scaled ( D.Height(), 1 );
  local_matrix_zero ( scaled );

  elem::AxpyInterface < double > globloc;
  globloc.Attach( elem::GLOBAL_TO_LOCAL, D );
  globloc.Axpy ( 1.0, scaled, 0, 0 );
  globloc.Detach();

  // scale eigenvalues and apply threshold

  double max = 0.0;
  double min;
  double val;
  double invmin;

  for ( size_t i = 0; i < scaled.Height(); ++i ) {
    val = scaled.Get( i, 0 );
    if ( val > max ) {
      max = val;
    }
  }

  switch ( op ) {
    case EIG_SQRT:
      min = sqrt(max) * threshold;
      for ( size_t i = 0; i < scaled.Height(); ++i ) {
        val = sqrt ( scaled.Get( i, 0 ) );
        val = ( val < min ) ? min : val;
        scaled.Set( i, 0, val );
      }
      break;
    case EIG_INVSQRT:
      min = sqrt(max) * threshold;
      invmin = 1.0 / min;
      for ( size_t i = 0; i < scaled.Height(); ++i ) {
        val = sqrt ( scaled.Get( i, 0 ) );
        val = ( val < min ) ? invmin : ( 1.0 / val );
        scaled.Set( i, 0, val );
      }
      break;
    default:
      min = max * threshold;
      for ( size_t i = 0; i < scaled.Height(); ++i ) {
        val = scaled.Get( i, 0 );
        val = ( val < min ) ? min : val;
        scaled.Set( i, 0, val );
      }
      break;
  }

  // Compute temp = W * op(D), by modifying our local elements. 

  mpi_matrix temp ( W );

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


void harp::mpi_column_norm ( mpi_matrix const & mat, mpi_matrix & S ) {

  S.ResizeTo ( mat.Height(), 1 );
  mpi_matrix_zero ( S );

  // Sum local columns

  elem_matrix_local Sloc ( S.Height(), 1 );
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


void harp::mpi_apply_norm ( mpi_matrix const & S, mpi_matrix & mat ) {

  // Get local copy of S

  elem_matrix_local Sloc ( S.Height(), 1 );
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


void harp::mpi_apply_inverse_norm ( mpi_matrix const & S, mpi_matrix & mat ) {

  // Get local copy of S

  elem_matrix_local Sloc ( S.Height(), 1 );
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


void harp::mpi_norm ( mpi_matrix const & D, mpi_matrix const & W, mpi_matrix & S ) {

  mpi_matrix temp ( W );
  mpi_matrix_zero ( temp );  

  mpi_eigen_compose ( EIG_SQRT, D, W, temp );

  mpi_column_norm ( temp, S );
  
  return;
}


void harp::mpi_gang_distribute ( mpi_matrix const & mat, mpi_matrix & gmat ) {

  if ( ( mat.Height() != gmat.Height() ) || ( mat.Width() != gmat.Width() ) ) {
    HARP_THROW( "matrix dimensions do not match" );
  }

  int gang_np;
  int gang_myp;

  MPI_Comm_size ( gmat.Grid().Comm(), &gang_np );
  MPI_Comm_rank ( gmat.Grid().Comm(), &gang_myp );

  // local column buffer

  elem_matrix_local colmat ( mat.Height(), 1 );

  // compute local gang-matrix elements

  elem::AxpyInterface < double > globloc;
  elem::AxpyInterface < double > locglob;

  // iterate over all columns

  mpi_matrix_zero ( gmat );

  for ( size_t col = 0; col < mat.Width(); ++col ) {

    local_matrix_zero ( colmat );

    // get full local copy of the column

    globloc.Attach( elem::GLOBAL_TO_LOCAL, mat );
    if ( gang_myp == 0 ) {
      globloc.Axpy ( 1.0, colmat, 0, col );
    }
    globloc.Detach();

    // distribute to gang

    locglob.Attach( elem::LOCAL_TO_GLOBAL, gmat );
    if ( gang_myp == 0 ) {
      locglob.Axpy ( 1.0, colmat, 0, col );
    }
    locglob.Detach();

  }

  return;
}


void harp::mpi_gang_accum ( mpi_matrix const & gmat, mpi_matrix & mat ) {

  if ( ( mat.Height() != gmat.Height() ) || ( mat.Width() != gmat.Width() ) ) {
    HARP_THROW( "matrix dimensions do not match" );
  }

  int gang_np;
  int gang_myp;

  MPI_Comm_size ( gmat.Grid().Comm(), &gang_np );
  MPI_Comm_rank ( gmat.Grid().Comm(), &gang_myp );

  // local column buffer

  elem_matrix_local colmat ( mat.Height(), 1 );

  // compute local gang-matrix elements

  elem::AxpyInterface < double > globloc;
  elem::AxpyInterface < double > locglob;

  // iterate over all columns

  mpi_matrix_zero ( mat );

  for ( size_t col = 0; col < mat.Width(); ++col ) {

    local_matrix_zero ( colmat );

    // get full local copy of the gang column

    globloc.Attach( elem::GLOBAL_TO_LOCAL, gmat );
    if ( gang_myp == 0 ) {
      globloc.Axpy ( 1.0, colmat, 0, col );
    }
    globloc.Detach();

    // each gang root process accumulates to global column

    locglob.Attach( elem::LOCAL_TO_GLOBAL, mat );
    if ( gang_myp == 0 ) {
      locglob.Axpy ( 1.0, colmat, 0, col );
    }
    locglob.Detach();

  }

  return;
}







