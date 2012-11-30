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



void harp::spec_project ( matrix_sparse const & psf, matrix_dist const & in, matrix_local & out ) {

  int np;
  int myp;

  MPI_Comm_size ( psf.Comm(), &np );
  MPI_Comm_rank ( psf.Comm(), &myp );

  // check consistent sizes

  size_t nbins = psf.Height();
  size_t npix = psf.Width();

  if ( out.Height() != npix ) {
    std::ostringstream o;
    o << "number of rows in output vector (" << out.Height() << ") does not match number of pixels in PSF (" << npix << ")";
    HARP_THROW( o.str().c_str() );
  }

  if ( in.Height() != nbins ) {
    std::ostringstream o;
    o << "number of rows in input vector (" << in.Height() << ") does not match number of bins in PSF (" << nbins << ")";
    HARP_THROW( o.str().c_str() );
  }

  // get local chunk of input vector which goes with our range of sparse matrix rows

  size_t local_firstrow = psf.FirstLocalRow();
  size_t local_rows = psf.LocalHeight();

  matrix_local local_in ( local_rows, 1 );
  local_matrix_zero ( local_in );

  elem::AxpyInterface < double > globloc;
  globloc.Attach( elem::GLOBAL_TO_LOCAL, in );
  globloc.Axpy ( 1.0, local_in, local_firstrow, 0 );
  globloc.Detach();

  // compute local output contribution

  matrix_local local_out ( npix, 1 );
  local_matrix_zero ( local_out );

  size_t col;
  double val;

  for ( size_t row = 0; row < local_rows; ++row ) {

    size_t off = psf.LocalEntryOffset ( row );
    size_t nnz = psf.NumConnections ( row );

    for ( size_t j = 0; j < nnz; ++j ) {
      col = psf.Col ( off + j );

      val = local_out.Get ( col, 0 );
      
      val += psf.Value ( off + j ) * local_in.Get ( row, 0 );

      local_out.Set ( col, 0, val );
    }

  }

  matrix_dist globout ( npix, 1 );
  dist_matrix_zero ( globout );

  elem::AxpyInterface < double > locglob;
  locglob.Attach( elem::LOCAL_TO_GLOBAL, globout );
  locglob.Axpy ( 1.0, local_out, 0, 0 );
  locglob.Detach();

  local_matrix_zero ( out );
  globloc.Attach( elem::GLOBAL_TO_LOCAL, globout );
  globloc.Axpy ( 1.0, out, 0, 0 );
  globloc.Detach();

  return;
}


void harp::noise_weighted_spec ( matrix_sparse const & psf, matrix_local const & invnoise, matrix_local const & img, matrix_dist & z ) {

  size_t nbins = psf.Height();
  size_t npix = psf.Width();

  if ( invnoise.Height() != npix ) {
    std::ostringstream o;
    o << "number of rows in inverse noise covariance (" << invnoise.Height() << ") does not match number of pixels in PSF (" << npix << ")";
    HARP_THROW( o.str().c_str() );
  }

  if ( img.Height() != npix ) {
    std::ostringstream o;
    o << "number of elements in image vector (" << img.Height() << ") does not match number of pixels in PSF (" << npix << ")";
    HARP_THROW( o.str().c_str() );
  }

  z.ResizeTo ( nbins, 1 );
  dist_matrix_zero ( z );

  // apply noise covariance to image

  matrix_local weight ( npix, 1 );
  for ( size_t i = 0; i < npix; ++i ) {
    weight.Set ( i, 0, invnoise.Get(i,0) * img.Get(i,0) );
  }

  // accumulate local pieces of z

  size_t first_loc_row = psf.FirstLocalRow();
  size_t loc_height = psf.LocalHeight();

  elem::AxpyInterface < double > locglob;
  locglob.Attach( elem::LOCAL_TO_GLOBAL, z );

  elem::Matrix < double > local_z ( loc_height, 1 );

  double val;
  size_t rowoff;
  size_t nnz;
  size_t col;

  for ( size_t loc_bin = 0; loc_bin < loc_height; ++loc_bin ) {
    rowoff = psf.LocalEntryOffset ( loc_bin );
    nnz = psf.NumConnections ( loc_bin );
    for ( size_t j = 0; j < nnz; ++j ) {
      val = psf.Value ( rowoff + j );
      col = psf.Col ( rowoff + j );
      local_z.Set ( loc_bin, 0, val * weight.Get(col,0) );
    }
  }

  locglob.Axpy ( 1.0, local_z, first_loc_row, 0 );

  locglob.Detach();

  return;
}


void harp::inverse_covariance ( matrix_sparse const & psf, matrix_local const & invnoise, matrix_dist & invcov ) {

  int np;
  int myp;

  MPI_Comm_size ( psf.Comm(), &np );
  MPI_Comm_rank ( psf.Comm(), &myp );

  // check consistent sizes

  size_t nbins = psf.Height();
  size_t npix = psf.Width();

  if ( invnoise.Height() != npix ) {
    std::ostringstream o;
    o << "number of rows in inverse noise covariance (" << invnoise.Height() << ") does not match number of pixels in PSF (" << npix << ")";
    HARP_THROW( o.str().c_str() );
  }

  invcov.ResizeTo ( nbins, nbins );
  dist_matrix_zero ( invcov );

  // First, accumulate our own diagonal subblock

  size_t local_firstrow = psf.FirstLocalRow();
  size_t local_rows = psf.LocalHeight();

  //cerr << "DBG:  proc " << myp << " has rows " << local_firstrow << " - " << local_firstrow + local_rows - 1 << endl;

  elem::Matrix < double > local_inv ( local_rows, local_rows );

  for ( size_t i = 0; i < local_rows; ++i ) {
    for ( size_t j = 0; j < local_rows; ++j ) {
      local_inv.Set ( j, i, 0.0 );
    }
  }

  //cerr << "DBG:  local " << local_rows << " x " << local_rows << " matrix allocated" << endl;

  elem::AxpyInterface < double > locglob;
  locglob.Attach( elem::LOCAL_TO_GLOBAL, invcov );

  double val;
  size_t lhs_off;
  size_t rhs_off;
  size_t lhs_nnz;
  size_t rhs_nnz;
  size_t lhs_col;
  size_t rhs_col;
  size_t j, k;

  //cerr << "DBG:  accumulate our own subblock" << endl;

  for ( size_t lhs_row = 0; lhs_row < local_rows; ++lhs_row ) {

    for ( size_t rhs_row = 0; rhs_row <= lhs_row; ++rhs_row ) {

      lhs_off = psf.LocalEntryOffset ( lhs_row );
      rhs_off = psf.LocalEntryOffset ( rhs_row );

      lhs_nnz = psf.NumConnections ( lhs_row );
      rhs_nnz = psf.NumConnections ( rhs_row );

      val = 0.0;

      k = 0;
      rhs_col = psf.Col ( rhs_off );

      for ( j = 0; j < lhs_nnz; ++j ) {

        lhs_col = psf.Col ( lhs_off + j );
        
        //cerr << " lhs/rhs = " << lhs_col << "/" << rhs_col << endl;
        while ( ( rhs_col < lhs_col ) && ( k < rhs_nnz - 1 ) ) {
          ++k;
          rhs_col = psf.Col ( rhs_off + k );
          //cerr << "  scroll rhs col to " << rhs_col << endl;
        }

        if ( rhs_col == lhs_col ) {
          //cerr << "  accumulate " << psf.Value ( lhs_off + j ) * psf.Value ( rhs_off + k ) << " to (" << lhs_row << ", " << rhs_row << ")" << endl; 
          val += psf.Value ( lhs_off + j ) * psf.Value ( rhs_off + k );
        }

      }

      local_inv.Set ( lhs_row, rhs_row, val );

      //cerr << "local_inv (" << lhs_row << ", " << rhs_row << ") set to " << val << endl;

    }
  }

  locglob.Axpy ( 1.0, local_inv, local_firstrow, local_firstrow );

  locglob.Detach();
  
  // In order to build up the final matrix, we must do (# procs / 2) communications of 
  // the distributed rows of the sparse PSF.  Send our data to the previous rank process
  // and receive from the next rank process.  Repeat this until all matrix blocks have 
  // been filled in.

  int nshift = (int)( np / 2 );

  int to_proc;
  if ( myp > 0 ) {
    to_proc = myp - 1;
  } else {
    to_proc = np - 1;
  }

  int from_proc;
  if ( myp == np - 1 ) {
    from_proc = 0;
  } else {
    from_proc = myp + 1;
  }

  char * sendbuf = NULL;
  size_t sendbytes;
  char * recvbuf = NULL;
  size_t recvbytes;
  MPI_Request send_size_request;
  MPI_Request send_request;
  MPI_Status status;

  //cerr << "DBG:  " << nshift << " shifts for " << np << " processes" << endl;

  for ( int shift = 0; shift < nshift; ++shift ) {

    if ( shift == 0 ) {
      // first shift, send our own data

      //cerr << "proc " << myp << " allocing sparse block for own data on shift 0" << endl;

      sparse_block * myblock = new sparse_block ( psf );

      sendbuf = myblock->pack ( sendbytes );

      delete ( myblock );

    } else {
      // pass along the buffer

      sendbytes = recvbytes;

      sendbuf = (char*)malloc ( sendbytes );
      if ( ! sendbuf ) {
        HARP_THROW( "cannot allocate send buffer" );
      }

      memcpy ( (void*)sendbuf, (void*)recvbuf, sendbytes );

      free ( recvbuf );
      
    }

    int temp = MPI_Barrier ( psf.Comm() );

    int send_size_key = (shift * 2 * np) + 2 * myp;
    int send_data_key = (shift * 2 * np) + 2 * myp + 1;
    int recv_size_key = (shift * 2 * np) + 2 * from_proc;
    int recv_data_key = (shift * 2 * np) + 2 * from_proc + 1;

    //cerr << "proc " << myp << " sending data size to " << to_proc << " with key " << send_size_key << endl;

    int ret = MPI_Isend ( (void*)(&sendbytes), 1, MPI_UNSIGNED_LONG, to_proc, send_size_key, psf.Comm(), &send_size_request );
    mpi_check ( psf.Comm(), ret );

    //cerr << "proc " << myp << " sending data to " << to_proc << " with key " << send_data_key << endl;

    ret = MPI_Isend ( (void*)sendbuf, sendbytes, MPI_CHAR, to_proc, send_data_key, psf.Comm(), &send_request );
    mpi_check ( psf.Comm(), ret );

    // receive block from sender

    //cerr << "proc " << myp << " receiving data size from " << from_proc << " with key " << recv_size_key << endl;

    ret = MPI_Recv ( (void*)(&recvbytes), 1, MPI_UNSIGNED_LONG, from_proc, recv_size_key, psf.Comm(), &status );
    mpi_check ( psf.Comm(), ret );

    //cerr << "proc " << myp << " allocing receive data" << endl;

    recvbuf = (char*)malloc ( recvbytes );
    if ( ! recvbuf ) {
      HARP_THROW( "cannot allocate receive buffer" );
    }

    //cerr << "proc " << myp << " receiving data from " << from_proc << " with key " << recv_data_key << endl;

    ret = MPI_Recv ( (void*)recvbuf, recvbytes, MPI_CHAR, from_proc, recv_data_key, psf.Comm(), &status );
    mpi_check ( psf.Comm(), ret );

    // reconstruct sparse_block

    //cerr << "proc " << myp << " reconstructing sparse block" << endl;

    sparse_block * other_block = new sparse_block ( recvbuf, recvbytes );

    // compute block

    //cerr << "proc " << myp << " computing block" << endl;

    locglob.Attach( elem::LOCAL_TO_GLOBAL, invcov );

    size_t axpy_row;
    size_t axpy_col;

    if ( from_proc > myp ) {
      // we are computing the transposed block of the output

      local_inv.ResizeTo ( other_block->local_rows, local_rows );

      axpy_row = other_block->local_firstrow;
      axpy_col = local_firstrow;

      for ( size_t lhs_row = 0; lhs_row < local_rows; ++lhs_row ) {

        for ( size_t rhs_row = 0; rhs_row < other_block->local_rows; ++rhs_row ) {

          lhs_off = psf.LocalEntryOffset ( lhs_row );
          rhs_off = other_block->local_row_offset [ rhs_row ];

          lhs_nnz = psf.NumConnections ( lhs_row );
          rhs_nnz = other_block->local_row_nnz [ rhs_row ];

          val = 0.0;

          k = 0;
          rhs_col = other_block->local_col[ rhs_off ];

          for ( j = 0; j < lhs_nnz; ++j ) {

            lhs_col = psf.Col ( lhs_off + j );
            
            while ( ( rhs_col < lhs_col ) && ( k < rhs_nnz - 1 ) ) {
              ++k;
              rhs_col = other_block->local_col[ rhs_off + k ];
            }

            if ( rhs_col == lhs_col ) {
              val += psf.Value ( lhs_off + j ) * other_block->data[ rhs_off + k ];
            }

          }

          local_inv.Set ( rhs_row, lhs_row, val );

        }
      }

      locglob.Axpy ( 1.0, local_inv, axpy_row, axpy_col );

    } else if ( ( np % 2 != 0 ) || ( shift != nshift -1 ) ) {
      // always compute non-transposed block if we have odd number of processes
      // or we have an even number of process and we are not on the last shift.

      local_inv.ResizeTo ( local_rows, other_block->local_rows );

      axpy_row = local_firstrow;
      axpy_col = other_block->local_firstrow;

      for ( size_t lhs_row = 0; lhs_row < local_rows; ++lhs_row ) {

        for ( size_t rhs_row = 0; rhs_row < other_block->local_rows; ++rhs_row ) {

          lhs_off = psf.LocalEntryOffset ( lhs_row );
          rhs_off = other_block->local_row_offset [ rhs_row ];

          lhs_nnz = psf.NumConnections ( lhs_row );
          rhs_nnz = other_block->local_row_nnz [ rhs_row ];

          val = 0.0;

          k = 0;
          rhs_col = other_block->local_col[ rhs_off ];

          for ( j = 0; j < lhs_nnz; ++j ) {

            lhs_col = psf.Col ( lhs_off + j );
            
            while ( ( rhs_col < lhs_col ) && ( k < rhs_nnz - 1 ) ) {
              ++k;
              rhs_col = other_block->local_col[ rhs_off + k ];
            }

            if ( rhs_col == lhs_col ) {
              val += psf.Value ( lhs_off + j ) * other_block->data[ rhs_off + k ];
            }

          }

          local_inv.Set ( lhs_row, rhs_row, val );

        }
      }

      locglob.Axpy ( 1.0, local_inv, axpy_row, axpy_col );

    }

    delete other_block;

    locglob.Detach();

    // free send buffer

    ret = MPI_Wait ( &send_size_request, &status );
    mpi_check ( psf.Comm(), ret );

    ret = MPI_Wait ( &send_request, &status );
    mpi_check ( psf.Comm(), ret );

    free ( sendbuf );

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

  out.ResizeTo ( W.Height(), W.Height() );

  // scale local portion of eigenvalues

  matrix_dist scaled ( D );

  size_t allocated = scaled.AllocatedMemory();
  double * raw = scaled.LocalBuffer();

  switch ( op ) {
    case EIG_SQRT:
      for ( size_t i = 0; i < allocated; ++i ) {
        raw[i] = sqrt ( raw[i] );
      }
      break;
    case EIG_INVSQRT:
      for ( size_t i = 0; i < allocated; ++i ) {
        raw[i] = 1.0 / sqrt ( raw[i] );
      }
      break;
    default:
      break;
  }

  // Get full local copy of eigenvalues

  matrix_local Dloc ( D.Height(), 1 );

  elem::AxpyInterface < double > globloc;
  globloc.Attach( elem::GLOBAL_TO_LOCAL, D );
  globloc.Axpy ( 1.0, Dloc, 0, 0 );
  globloc.Detach();

  // Compute temp = op(D) * W, by modifying our local elements

  matrix_dist temp ( W );

  matrix_local & local = temp.LocalMatrix();

  int hlocal = temp.LocalHeight();
  int wlocal = temp.LocalWidth();

  int rowoff = temp.ColShift();
  int rowstride = temp.ColStride();
  int row;

  double mval;

  for ( int i = 0; i < wlocal; ++i ) {
    for ( int j = 0; j < hlocal; ++j ) {
      row = rowoff + j * rowstride;
      mval = local.Get ( j, i );
      mval *= Dloc.Get ( row, 0 );
      local.Set ( j, i, mval );
    }
  }

  // compute out = W^T * ( op(D) * W )

  elem::Gemm ( elem::TRANSPOSE, elem::NORMAL, 1.0, W, temp, 0.0, out );

  return;
}


void harp::column_norm ( matrix_dist & mat, matrix_dist & S ) {

  S.ResizeTo ( mat.Height(), 1 );
  dist_matrix_zero ( S );

  // Sum local columns

  matrix_local Sloc ( S.Height(), 1 );

  matrix_local & local = mat.LocalMatrix();

  int hlocal = mat.LocalHeight();
  int wlocal = mat.LocalWidth();

  int rowoff = mat.ColShift();
  int rowstride = mat.ColStride();
  int row;

  double mval;

  for ( int i = 0; i < wlocal; ++i ) {
    for ( int j = 0; j < hlocal; ++j ) {
      row = rowoff + j * rowstride;
      mval = local.Get ( j, i ) + Sloc.Get ( row, 0 );
      Sloc.Set ( row, 0, mval );
    }
  }

  // Reduce

  elem::AxpyInterface < double > locglob;
  locglob.Attach( elem::LOCAL_TO_GLOBAL, S );
  locglob.Axpy ( 1.0, Sloc, 0, 0 );
  locglob.Detach();

  // Invert

  size_t allocated = S.AllocatedMemory();
  double * raw = S.LocalBuffer();

  for ( size_t i = 0; i < allocated; ++i ) {
    raw[i] = 1.0 / raw[i];
  }

  return;
}


void harp::apply_norm ( matrix_dist & S, matrix_dist & mat ) {

  // Get local copy of S

  matrix_local Sloc ( S.Height(), 1 );

  elem::AxpyInterface < double > globloc;
  globloc.Attach( elem::GLOBAL_TO_LOCAL, S );
  globloc.Axpy ( 1.0, Sloc, 0, 0 );
  globloc.Detach();

  // Apply to local matrix

  matrix_local & local = mat.LocalMatrix();

  int hlocal = mat.LocalHeight();
  int wlocal = mat.LocalWidth();

  int rowoff = mat.ColShift();
  int rowstride = mat.ColStride();
  int row;

  double mval;

  for ( int i = 0; i < wlocal; ++i ) {
    for ( int j = 0; j < hlocal; ++j ) {
      row = rowoff + j * rowstride;
      mval = local.Get ( j, i );
      mval *= Sloc.Get ( row, 0 );
      local.Set ( j, i, mval );
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


void harp::resolution ( matrix_dist & D, matrix_dist & W, matrix_dist & S, matrix_dist & R ) {

  R.ResizeTo ( W.Height(), W.Height() );
  dist_matrix_zero ( R );

  eigen_compose ( EIG_SQRT, D, W, R );

  column_norm ( R, S );

  apply_norm ( S, R );

  return;
}


void harp::extract ( matrix_dist & D, matrix_dist & W, matrix_dist & S, matrix_dist & z, matrix_dist & f ) {

  // compose ( W^T D^{-1/2} W )

  matrix_dist rtC ( W );
  dist_matrix_zero ( rtC );

  eigen_compose ( EIG_INVSQRT, D, W, rtC );

  // multiply C^{1/2} * z

  matrix_dist vtemp ( z );
  dist_matrix_zero ( vtemp );

  elem::Gemv ( elem::NORMAL, 1.0, rtC, z, 0.0, vtemp ); 

  // multiply S^-1 * vtemp

  f = vtemp;

  apply_norm ( S, f );

  return;
}




/*

  template < class P, class S, class R >
  void append_sky_comp ( boost::numeric::ublas::matrix_expression < P > const & psf, boost::numeric::ublas::matrix_expression < S > const & sky, size_t nspec, boost::numeric::ublas::matrix_expression < R > & output ) {

    typedef boost::numeric::ublas::matrix_range < R > R_view;


    // verify that length of sky components equals number of columns in design matrix

    size_t npix = psf().size1();
    size_t nbins = psf().size2();
    size_t ncomp = sky().size1();
    size_t nspecbins = (size_t) ( nbins / nspec );

    if ( sky().size2() != nspecbins ) {
      std::ostringstream o;
      o << "number of columns in sky component matrix (" << sky().size2() << ") does not match number of bins per spectrum in PSF (" << nspecbins << ")";
      MOAT_THROW( o.str().c_str() );
    }

    if ( ncomp > nspecbins ) {
      std::ostringstream o;
      o << "number of sky components (" << ncomp << ") exceeds the number of flux bins per spectrum (" << nbins << ")";
      MOAT_THROW( o.str().c_str() );
    }

    // resize output matrix and copy the original design matrix into the
    // first column block

    output().resize ( npix, nbins + ncomp * nspecbins );

    R_view original ( output(), mv_range( 0, npix ), mv_range( 0, nbins ) );

    original.assign ( psf() );

    // construct component matrix

    mat_compcol compmat ( nbins, ncomp * nspec, ncomp * nbins );

    for ( size_t i = 0; i < nspec; ++i ) {
      for ( size_t j = 0; j < ncomp; ++j ) {
        for ( size_t k = 0; k < nspecbins; ++k ) {
          compmat ( i * nspecbins + k, i * ncomp + j ) = sky() ( j, k );
        }
      }
    }

    // fill in sky portion of output design matrix

    R_view skyblock ( output(), mv_range ( 0, npix ), mv_range ( nbins, nbins + ncomp * nspecbins ) );

    boost::numeric::ublas::axpy_prod ( psf(), compmat, skyblock, true );

    return;
  }


*/


