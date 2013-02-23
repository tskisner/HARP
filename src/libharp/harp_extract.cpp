// @COPYRIGHT@

#include <harp_internal.hpp>


using namespace std;
using namespace harp;


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

  R.ResizeTo ( W.Height(), W.Height() );
  dist_matrix_zero ( R );

  eigen_compose ( EIG_SQRT, D, W, R );

  column_norm ( R, S );

  apply_norm ( S, R );

  return;
}


void harp::extract ( matrix_dist & D, matrix_dist & W, matrix_dist & S, matrix_dist & z, matrix_dist & f ) {

  dist_matrix_zero ( f );

  // compose ( W^T D^{-1/2} W )

  matrix_dist rtC ( W );
  dist_matrix_zero ( rtC );

  eigen_compose ( EIG_INVSQRT, D, W, rtC );

  // multiply C^{1/2} * z

  elem::Symv ( elem::LOWER, 1.0, rtC, z, 0.0, f ); 

  // multiply S^-1 * vtemp

  apply_norm ( S, f );

  return;
}


// Append columns to the design matrix to solve for a coefficient for the sky component

/*
void harp::sky_append ( matrix_sparse const & psf, size_t nspec, size_t specsize, matrix_local & skyspec, matrix_sparse & fullpsf ) {

  size_t nbins_orig = psf.Height();
  size_t npix = psf.Width();

  size_t nbins = nbins_orig + nspec;

  fullpsf.ResizeTo ( nbins, npix );

  



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
*/

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


