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

  block_.vals = 0;
}


void harp::mpi_matrix_sparse::clear ( ) {
  block_.vals = 0;
  block_.row.resize(0);
  block_.col.resize(0);
  block_.row_offset.resize(0);
  block_.row_nnz.resize(0);
  block_.data.resize(0);
  return;
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


void harp::mpi_eigen_decompose ( mpi_matrix const & invcov, mpi_matrix & D, mpi_matrix & W, bool regularize ) {

  D.Resize ( invcov.Height(), 1 );
  W.Resize ( invcov.Height(), invcov.Height() );

  mpi_matrix temp ( invcov );

  elem::DistMatrix < double, elem::VR, elem::STAR > eigvals ( invcov.Height(), 1, invcov.Grid() );

  elem::HermitianEig ( elem::LOWER, temp, eigvals, W );

  // I don't think we need this, since we never need to order the
  // eigenpairs.
  //elem::SortEig( eigvals, W );

  D = eigvals;

  if ( regularize ) {

    // find local min / max eigenvalues

    double min = 1.0e100;
    double max = -1.0e100;

    size_t hlocal = D.LocalHeight();
    size_t wlocal = D.LocalWidth();

    size_t rowoff = D.ColShift();
    size_t rowstride = D.ColStride();
    size_t row;

    double val;

    if ( wlocal > 0 ) {
      for ( size_t j = 0; j < hlocal; ++j ) {
        // the global element of sub_data
        row = rowoff + j * rowstride;

        val = D.GetLocal ( j, 0 );

        if ( val < min ) {
          min = val;
        }
        if ( val > max ) {
          max = val;
        }
      }
    }

    // find global min / max

    elem::mpi::AllReduce ( &min, 1, elem::mpi::MIN, D.Grid().Comm() );

    elem::mpi::AllReduce ( &max, 1, elem::mpi::MAX, D.Grid().Comm() );

    // compute condition number

    double rcond = min / max;

    // pick some delta that is bigger than machine precision, but still tiny
    double epsilon = 10.0 * std::numeric_limits < double > :: epsilon();

    // modify our local elements

    if ( rcond < epsilon ) {

      double reg = max * epsilon - min;

      if ( wlocal > 0 ) {
        for ( size_t j = 0; j < hlocal; ++j ) {
          // the global element of sub_data
          row = rowoff + j * rowstride;

          val = D.GetLocal ( j, 0 );
          val += reg;
          D.SetLocal ( j, 0, val );
        }
      }

    }

  }

  return;
}


void harp::mpi_eigen_compose ( eigen_op op, mpi_matrix const & D, mpi_matrix const & W, mpi_matrix & out ) {

  double threshold = 1.0e-12;

  out.Resize ( W.Height(), W.Height() );

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
    case EIG_INV:
      min = max * threshold;
      invmin = 1.0 / min;
      for ( size_t i = 0; i < scaled.Height(); ++i ) {
        val = scaled.Get( i, 0 );
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

  S.Resize ( mat.Height(), 1 );
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


void harp::mpi_sparse_mv_trans ( mpi_matrix_sparse const & AT, mpi_matrix const & in, elem_matrix_local & out ) {

  int myp = AT.comm().rank();
  int np = AT.comm().size();

  // check consistent sizes

  size_t nrows = AT.rows();
  size_t ncols = AT.cols();

  if ( in.Height() != nrows ) {
    std::ostringstream o;
    o << "number of rows in input vector (" << in.Height() << ") does not match number of rows in transposed matrix (" << nrows << ")";
    HARP_MPI_ABORT( myp, o.str().c_str() );
  }

  // get local chunk of input vector which goes with our range of sparse matrix rows

  size_t local_firstrow = AT.block().firstrow;
  size_t local_nrows = AT.block().rows;
  size_t local_nvals = AT.block().vals;

  elem_matrix_local local_in;
  if ( local_nrows > 0 ) {
    local_in.Resize ( local_nrows, 1 );
  }
  local_matrix_zero ( local_in );

  elem::AxpyInterface < double > globloc;
  globloc.Attach( elem::GLOBAL_TO_LOCAL, in );
  if ( local_nrows > 0 ) {
    globloc.Axpy ( 1.0, local_in, local_firstrow, 0 );
  }
  globloc.Detach();

  // compute local output contribution

  out.Resize ( ncols, 1 );
  local_matrix_zero ( out );

  double val;
  double inval;
  double outval;
  size_t row;
  size_t col;

  for ( size_t loc = 0; loc < local_nvals; ++loc ) {
    row = AT.block().row[ loc ];
    col = AT.block().col[ loc ];
    val = AT.block().data[ loc ];
    inval = local_in.Get ( row, 0 );
    outval = out.Get ( col, 0 );
    out.Set ( col, 0, outval + inval * val );
  }

  // accumulate to global solution

  mpi_matrix globout ( ncols, 1, in.Grid() );
  mpi_matrix_zero ( globout );

  elem::AxpyInterface < double > locglob;
  locglob.Attach( elem::LOCAL_TO_GLOBAL, globout );
  locglob.Axpy ( 1.0, out, 0, 0 );
  locglob.Detach();

  // get local copy for output
  
  local_matrix_zero ( out );
  globloc.Attach( elem::GLOBAL_TO_LOCAL, globout );
  globloc.Axpy ( 1.0, out, 0, 0 );
  globloc.Detach();

  return;
}


void harp::mpi_gang_distribute ( mpi_matrix const & mat, mpi_matrix & gmat ) {

  int gang_np;
  int gang_myp;

  MPI_Comm_size ( gmat.Grid().Comm(), &gang_np );
  MPI_Comm_rank ( gmat.Grid().Comm(), &gang_myp );

  if ( ( mat.Height() != gmat.Height() ) || ( mat.Width() != gmat.Width() ) ) {
    HARP_MPI_ABORT( gang_myp, "matrix dimensions do not match" );
  }

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

  int gang_np;
  int gang_myp;

  MPI_Comm_size ( gmat.Grid().Comm(), &gang_np );
  MPI_Comm_rank ( gmat.Grid().Comm(), &gang_myp );

  if ( ( mat.Height() != gmat.Height() ) || ( mat.Width() != gmat.Width() ) ) {
    HARP_MPI_ABORT( gang_myp, "matrix dimensions do not match" );
  }

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







