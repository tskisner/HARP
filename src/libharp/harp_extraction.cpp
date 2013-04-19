// @COPYRIGHT@

#include <harp_internal.hpp>


using namespace std;
using namespace harp;



void harp::sub_spec ( matrix_dist & in, size_t total_nspec, size_t first_spec, size_t nspec, size_t first_lambda, size_t nlambda, matrix_dist & out ) {

  size_t in_nlambda = (size_t) ( in.Height() / total_nspec );

  if ( first_spec + nspec > total_nspec ) {
    std::ostringstream o;
    o << "input spec range (" << first_spec << " - " << first_spec + nspec - 1 << ") exceeds the number of input spectra (" << total_nspec << ")";
    HARP_THROW( o.str().c_str() );
  }

  if ( first_lambda + nlambda > in_nlambda ) {
    std::ostringstream o;
    o << "input lambda range (" << first_lambda << " - " << first_lambda + nlambda - 1 << ") exceeds the number of input points (" << in_nlambda << ")";
    HARP_THROW( o.str().c_str() );
  }

  if ( (size_t)out.Height() != nspec * nlambda ) {
    std::ostringstream o;
    o << "output matrix height (" << out.Height() << ") does not match parameters (" << nspec * nlambda << ")";
    HARP_THROW( o.str().c_str() );
  }

  dist_matrix_zero ( out );

  // FIXME: this should be changed to not store a full local copy of the input!

  matrix_local in_loc ( in.Height(), 1 );
  local_matrix_zero ( in_loc );

  elem::AxpyInterface < double > globloc;
  globloc.Attach( elem::GLOBAL_TO_LOCAL, in );
  globloc.Axpy ( 1.0, in_loc, 0, 0 );
  globloc.Detach();

  // Update local output matrix with proper slices from the input

  size_t hlocal = out.LocalHeight();
  size_t wlocal = out.LocalWidth();

  size_t rowoff = out.ColShift();
  size_t rowstride = out.ColStride();
  size_t row;

  double val;
  size_t out_spec;
  size_t out_lambda;

  if ( wlocal > 0 ) {
    for ( size_t j = 0; j < hlocal; ++j ) {
      row = rowoff + j * rowstride;
      out_spec = (size_t)( row / nlambda );
      out_lambda = row - out_spec * nlambda;
      out_spec += first_spec;
      out_lambda += first_lambda;
      val = in_loc.Get ( out_spec * in_nlambda + out_lambda, 0 );
      out.SetLocal ( j, 0, val );
    }
  }

  return;
}


void harp::accum_spec ( matrix_dist & full, size_t total_nspec, size_t first_spec, size_t nspec, size_t first_lambda, size_t nlambda, matrix_dist & chunk ) {

  size_t full_nlambda = (size_t) ( full.Height() / total_nspec );

  if ( first_spec + nspec > total_nspec ) {
    std::ostringstream o;
    o << "chunk spec range (" << first_spec << " - " << first_spec + nspec - 1 << ") exceeds the number of full spectra (" << total_nspec << ")";
    HARP_THROW( o.str().c_str() );
  }

  if ( first_lambda + nlambda > full_nlambda ) {
    std::ostringstream o;
    o << "chunk lambda range (" << first_lambda << " - " << first_lambda + nlambda - 1 << ") exceeds the number of full points (" << full_nlambda << ")";
    HARP_THROW( o.str().c_str() );
  }

  if ( (size_t)chunk.Height() != nspec * nlambda ) {
    std::ostringstream o;
    o << "chunk matrix height (" << chunk.Height() << ") does not match parameters (" << nspec * nlambda << ")";
    HARP_THROW( o.str().c_str() );
  }

  // FIXME: this should be changed to not store a full local copy of the matrix!

  matrix_local full_loc ( full.Height(), 1 );
  local_matrix_zero ( full_loc );

  // Copy our data into full local contribution

  size_t hlocal = chunk.LocalHeight();
  size_t wlocal = chunk.LocalWidth();

  size_t rowoff = chunk.ColShift();
  size_t rowstride = chunk.ColStride();
  size_t row;

  double val;
  size_t spec;
  size_t lambda;

  if ( wlocal > 0 ) {
    for ( size_t j = 0; j < hlocal; ++j ) {
      row = rowoff + j * rowstride;
      spec = (size_t)( row / nlambda );
      lambda = row - spec * nlambda;
      spec += first_spec;
      lambda += first_lambda;
      val = chunk.GetLocal ( j, 0 );
      full_loc.Set ( spec * nlambda + lambda, 0, val );
    }
  }

  // accumulate

  elem::AxpyInterface < double > locglob;
  locglob.Attach( elem::LOCAL_TO_GLOBAL, full );
  locglob.Axpy ( 1.0, full_loc, 0, 0 );
  locglob.Detach();

  return;
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

  size_t first_loc_row = psf.FirstLocalRow();
  size_t loc_height = psf.LocalHeight();
  size_t loc_entries = psf.NumLocalEntries();

  matrix_local local_in ( loc_height, 1 );
  local_matrix_zero ( local_in );

  elem::AxpyInterface < double > globloc;
  globloc.Attach( elem::GLOBAL_TO_LOCAL, in );
  globloc.Axpy ( 1.0, local_in, first_loc_row, 0 );
  globloc.Detach();

  // compute local output contribution

  matrix_local local_out ( npix, 1 );
  local_matrix_zero ( local_out );

  double val;
  double inval;
  double outval;
  size_t row;
  size_t col;

  for ( size_t loc = 0; loc < loc_entries; ++loc ) {
    row = psf.Row ( loc );
    col = psf.Col ( loc );
    val = psf.Value ( loc );
    inval = local_in.Get ( row - first_loc_row, 0 );
    outval = local_out.Get ( col, 0 );
    local_out.Set ( col, 0, outval + inval * val );
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
  size_t loc_entries = psf.NumLocalEntries();

  matrix_local local_z ( loc_height, 1 );
  local_matrix_zero ( local_z );

  double val;
  double zval;
  size_t row;
  size_t nnz;
  size_t col;

  for ( size_t loc = 0; loc < loc_entries; ++loc ) {
    row = psf.Row ( loc );
    col = psf.Col ( loc );
    val = psf.Value ( loc );
    zval = local_z.Get ( row - first_loc_row, 0 );
    local_z.Set ( row - first_loc_row, 0, zval + val * weight.Get( col, 0 ) );
  }

  elem::AxpyInterface < double > locglob;
  locglob.Attach( elem::LOCAL_TO_GLOBAL, z );
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

  elem::Matrix < double > local_inv ( local_rows, local_rows );

  for ( size_t i = 0; i < local_rows; ++i ) {
    for ( size_t j = 0; j < local_rows; ++j ) {
      local_inv.Set ( j, i, 0.0 );
    }
  }

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
        
        while ( ( rhs_col < lhs_col ) && ( k < rhs_nnz - 1 ) ) {
          ++k;
          rhs_col = psf.Col ( rhs_off + k );
        }

        if ( rhs_col == lhs_col ) {
          val += invnoise.Get( lhs_col, 0 ) * psf.Value ( lhs_off + j ) * psf.Value ( rhs_off + k );
        }

      }

      local_inv.Set ( lhs_row, rhs_row, val );

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

  for ( int shift = 0; shift < nshift; ++shift ) {

    if ( shift == 0 ) {
      // first shift, send our own data

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

    int ret = MPI_Isend ( (void*)(&sendbytes), 1, MPI_UNSIGNED_LONG, to_proc, send_size_key, psf.Comm(), &send_size_request );
    mpi_check ( psf.Comm(), ret );

    ret = MPI_Isend ( (void*)sendbuf, sendbytes, MPI_CHAR, to_proc, send_data_key, psf.Comm(), &send_request );
    mpi_check ( psf.Comm(), ret );

    // receive block from sender

    ret = MPI_Recv ( (void*)(&recvbytes), 1, MPI_UNSIGNED_LONG, from_proc, recv_size_key, psf.Comm(), &status );
    mpi_check ( psf.Comm(), ret );

    recvbuf = (char*)malloc ( recvbytes );
    if ( ! recvbuf ) {
      HARP_THROW( "cannot allocate receive buffer" );
    }

    ret = MPI_Recv ( (void*)recvbuf, recvbytes, MPI_CHAR, from_proc, recv_data_key, psf.Comm(), &status );
    mpi_check ( psf.Comm(), ret );

    // reconstruct sparse_block

    sparse_block * other_block = new sparse_block ( recvbuf, recvbytes );

    // compute block

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
              val += invnoise.Get( lhs_col, 0 ) * psf.Value ( lhs_off + j ) * other_block->data[ rhs_off + k ];
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
              val += invnoise.Get( lhs_col, 0 ) * psf.Value ( lhs_off + j ) * other_block->data[ rhs_off + k ];
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


void harp::resolution ( matrix_dist & D, matrix_dist & W, matrix_dist & S, matrix_dist & R ) {

  R = W;

  dist_matrix_zero ( R );

  eigen_compose ( EIG_SQRT, D, W, R );

  column_norm ( R, S );

  apply_norm ( S, R );

  return;
}


void harp::extract ( matrix_dist & D, matrix_dist & W, matrix_dist & S, matrix_dist & z, matrix_dist & Rf, matrix_dist & err, matrix_dist & f ) {

  dist_matrix_zero ( f );
  dist_matrix_zero ( Rf );
  dist_matrix_zero ( err );

  // compose ( W^T D^{-1/2} W )

  matrix_dist rtC ( W );
  dist_matrix_zero ( rtC );

  eigen_compose ( EIG_INVSQRT, D, W, rtC );

  // Store R * C for error calculation before calling Symv (which may
  // destroy upper triangle).

  matrix_dist RC ( rtC );

  apply_norm ( S, RC );

  // multiply rtC * z

  elem::Symv ( elem::LOWER, 1.0, rtC, z, 0.0, Rf );

  // normalize output spectra

  apply_norm ( S, Rf );


  // compute deconvolved spectra (numerically unstable, but useful for visualization).
  // R^-1 == ( W^T D^{-1/2} W ) S

  // first apply inverse norm

  matrix_dist srf ( Rf );

  apply_inverse_norm ( S, srf );

  // now apply rtC

  elem::Symv ( elem::LOWER, 1.0, rtC, srf, 0.0, f );


  // compute diagonal error on result.

  // set local err matrix copy

  matrix_local err_loc ( err.Height(), 1 );

  matrix_local & rc_loc = RC.Matrix();

  size_t hlocal = RC.LocalHeight();
  size_t wlocal = RC.LocalWidth();

  size_t rowoff = RC.ColShift();
  size_t rowstride = RC.ColStride();
  size_t row;

  size_t coloff = RC.RowShift();
  size_t colstride = RC.RowStride();
  size_t col;

  double mval;

  for ( size_t i = 0; i < wlocal; ++i ) {
    for ( size_t j = 0; j < hlocal; ++j ) {
      row = rowoff + j * rowstride;
      col = coloff + i * colstride;
      if ( row == col ) {
        err_loc.Set ( row, 0, sqrt ( RC.GetLocal ( j, i ) ) );
      }
    }
  }

  // reduce local error vector to global one.  There should be no overlap
  // of non-zero elements, since each diagonal value is owned by one process.

  dist_matrix_zero ( err );

  elem::AxpyInterface < double > locglob;
  locglob.Attach( elem::LOCAL_TO_GLOBAL, err );
  locglob.Axpy ( 1.0, err_loc, 0, 0 );
  locglob.Detach();

  return;
}


// Produce a design matrix suitable for simultaneous extraction
// and sky subtraction

void harp::sky_design ( matrix_sparse const & AT, std::vector < bool > const & sky, matrix_sparse & skyAT ) {

  size_t nbins = AT.Height();
  size_t npix = AT.Width();
  size_t nspec = sky.size();

  size_t nsky = 0;
  for ( size_t i = 0; i < nspec; ++i ) {
    if ( sky[i] ) {
      ++nsky;
    }
  }

  size_t nlambda = (size_t)( nbins / nspec );

  if ( nspec * nlambda != nbins ) {
    std::ostringstream o;
    o << "sky vector size (" << nspec << ") does not divide evenly into PSF spectral dimension (" << nbins << ")";
    HARP_THROW( o.str().c_str() );
  }

  size_t sky_nspec = nspec - nsky + 1;
  size_t sky_nbins = sky_nspec * nlambda;

  // build the full mapping of old bins to new bins

  vector < size_t > old_to_new ( nbins );

  size_t skystart = (nspec - nsky) * nlambda;
  size_t newoff = 0;
  size_t oldoff = 0;

  for ( size_t i = 0; i < nspec; ++i ) {
    if ( sky[i] ) {
      for ( size_t j = 0; j < nlambda; ++j ) {
        old_to_new[ oldoff ] = skystart + j;
        ++oldoff;
      }
    } else {
      for ( size_t j = 0; j < nlambda; ++j ) {
        old_to_new[ oldoff ] = newoff;
        ++oldoff;
        ++newoff;
      }
    }
  }

  // set up new design matrix and get our local range

  skyAT.ResizeTo ( sky_nbins, npix );

  size_t sky_firstrow = skyAT.FirstLocalRow();
  size_t sky_rows = skyAT.LocalHeight();

  size_t firstrow = AT.FirstLocalRow();
  size_t rows = AT.LocalHeight();

  // accumulate our local block of the original design matrix.  Since we do not
  // know a priori how many non-zeroes we will have

  skyAT.StartAssembly();

  // compute number of updates and reserve more space if needed

  size_t updates = 0;
  for ( size_t i = 0; i < rows; ++i ) {

  }

  /*

  // In order to build up the new matrix, send our data to the previous rank process
  // and receive from the next rank process.  Repeat this (number of process times)
  // to guarantee that all rows have been filled in.  

  int nshift = np - 1;

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

  for ( int shift = 0; shift < nshift; ++shift ) {

    if ( shift == 0 ) {
      // first shift, send our own data

      sparse_block * myblock = new sparse_block ( AT );

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

    int temp = MPI_Barrier ( AT.Comm() );

    int send_size_key = (shift * 2 * np) + 2 * myp;
    int send_data_key = (shift * 2 * np) + 2 * myp + 1;
    int recv_size_key = (shift * 2 * np) + 2 * from_proc;
    int recv_data_key = (shift * 2 * np) + 2 * from_proc + 1;

    int ret = MPI_Isend ( (void*)(&sendbytes), 1, MPI_UNSIGNED_LONG, to_proc, send_size_key, AT.Comm(), &send_size_request );
    mpi_check ( AT.Comm(), ret );

    ret = MPI_Isend ( (void*)sendbuf, sendbytes, MPI_CHAR, to_proc, send_data_key, AT.Comm(), &send_request );
    mpi_check ( AT.Comm(), ret );

    // receive block from sender

    ret = MPI_Recv ( (void*)(&recvbytes), 1, MPI_UNSIGNED_LONG, from_proc, recv_size_key, AT.Comm(), &status );
    mpi_check ( AT.Comm(), ret );

    recvbuf = (char*)malloc ( recvbytes );
    if ( ! recvbuf ) {
      HARP_THROW( "cannot allocate receive buffer" );
    }

    ret = MPI_Recv ( (void*)recvbuf, recvbytes, MPI_CHAR, from_proc, recv_data_key, AT.Comm(), &status );
    mpi_check ( AT.Comm(), ret );

    // reconstruct sparse_block

    sparse_block * other_block = new sparse_block ( recvbuf, recvbytes );

    // compute block

    locglob.Attach( elem::LOCAL_TO_GLOBAL, invcov );

    size_t axpy_row;
    size_t axpy_col;

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
            val += invnoise.Get( lhs_col, 0 ) * psf.Value ( lhs_off + j ) * other_block->data[ rhs_off + k ];
          }

        }

        local_inv.Set ( lhs_row, rhs_row, val );

      }
    }

    locglob.Axpy ( 1.0, local_inv, axpy_row, axpy_col );



    delete other_block;

    locglob.Detach();

    // free send buffer

    ret = MPI_Wait ( &send_size_request, &status );
    mpi_check ( psf.Comm(), ret );

    ret = MPI_Wait ( &send_request, &status );
    mpi_check ( psf.Comm(), ret );

    free ( sendbuf );

  }

  */


  return;
}



