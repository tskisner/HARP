// @COPYRIGHT@

#include <harp_mpi_internal.hpp>


using namespace std;
using namespace harp;



void harp::mpi_sub_spec ( spec_slice_region const & full_region, spec_slice_region const & sub_region, mpi_matrix const & full_data, bool use_good_sub, mpi_matrix & sub_data ) {

  // verify that data dimensions match the sizes of the regions.  Select whether we
  // are using the full or good extent of the output. 

  if ( (full_region.n_spec * full_region.n_lambda) != full_data.Height() ) {
    std::ostringstream o;
    o << "input region total bins (" << (full_region.n_spec * full_region.n_lambda) << ") does not match input data size (" << full_data.Height() << ")";
    HARP_THROW( o.str().c_str() );
  }

  if ( (sub_region.n_spec * sub_region.n_lambda) != sub_data.Height() ) {
    std::ostringstream o;
    o << "output region total bins (" << (sub_region.n_spec * sub_region.n_lambda) << ") does not match output data size (" << sub_data.Height() << ")";
    HARP_THROW( o.str().c_str() );
  }

  size_t nsub;

  size_t sub_nlambda;
  size_t sub_nspec;
  size_t sub_firstspec;
  size_t sub_firstlambda;

  if ( use_good_sub ) {
    nsub = sub_region.n_good_spec * sub_region.n_good_lambda;
    sub_nlambda = sub_region.n_good_lambda;
    sub_nspec = sub_region.n_good_spec;
    sub_firstspec = sub_region.first_good_spec;
    sub_firstlambda = sub_region.first_good_lambda;
  } else {
    nsub = sub_region.n_spec * sub_region.n_lambda;
    sub_nlambda = sub_region.n_lambda;
    sub_nspec = sub_region.n_spec;
    sub_firstspec = sub_region.first_spec;
    sub_firstlambda = sub_region.first_lambda;
  }

  if ( sub_firstspec < full_region.first_spec ) {
    std::ostringstream o;
    o << "sub region first spec (" << sub_firstspec << ") is before first spec of full region (" << full_region.first_spec << ")";
    HARP_THROW( o.str().c_str() );
  }

  if ( sub_firstspec + sub_nspec > full_region.first_spec + full_region.n_spec ) {
    std::ostringstream o;
    o << "sub region last spec (" << (sub_firstspec + sub_nspec - 1) << ") is beyond last spec of full region (" << (full_region.first_spec + full_region.n_spec - 1) << ")";
    HARP_THROW( o.str().c_str() );
  }

  if ( sub_firstlambda < full_region.first_lambda ) {
    std::ostringstream o;
    o << "sub region first lambda (" << sub_firstlambda << ") is before first lambda of full region (" << full_region.first_lambda << ")";
    HARP_THROW( o.str().c_str() );
  }

  if ( sub_firstlambda + sub_nlambda > full_region.first_lambda + full_region.n_lambda ) {
    std::ostringstream o;
    o << "sub region last lambda (" << (sub_firstlambda + sub_nlambda - 1) << ") is beyond last lambda of full region (" << (full_region.first_lambda + full_region.n_lambda - 1) << ")";
    HARP_THROW( o.str().c_str() );
  }

  /*
  cerr << "DBG:  sub_spec" << endl;
  cerr << "DBG:     n_spec = " << sub_region.n_spec << endl;
  cerr << "DBG:     n_good_spec = " << sub_region.n_good_spec << endl;
  cerr << "DBG:     first_spec = " << sub_region.first_spec << endl;
  cerr << "DBG:     first_good_spec = " << sub_region.first_good_spec << endl;
  cerr << "DBG:     overlap_spec = " << sub_region.overlap_spec << endl;
  cerr << "DBG:     n_lambda = " << sub_region.n_lambda << endl;
  cerr << "DBG:     n_good_lambda = " << sub_region.n_good_lambda << endl;
  cerr << "DBG:     first_lambda = " << sub_region.first_lambda << endl;
  cerr << "DBG:     first_good_lambda = " << sub_region.first_good_lambda << endl;
  cerr << "DBG:     overlap_lambda = " << sub_region.overlap_lambda << endl;
  cerr << "DBG:     select_nspec = " << sub_nspec << endl;
  cerr << "DBG:     select_nlambda = " << sub_nlambda << endl;
  cerr << "DBG:     select_firstspec = " << sub_firstspec << endl;
  cerr << "DBG:     select_firstlambda = " << sub_firstlambda << endl;
  */

  // clear output

  mpi_matrix_zero ( sub_data );

  // FIXME: can this be changed to not store a full local copy of the input?

  elem_matrix_local full_loc ( full_data.Height(), 1 );
  local_matrix_zero ( full_loc );

  elem::AxpyInterface < double > globloc;
  globloc.Attach( elem::GLOBAL_TO_LOCAL, full_data );
  globloc.Axpy ( 1.0, full_loc, 0, 0 );
  globloc.Detach();

  // Update local output matrix with proper slices from the input

  size_t hlocal = sub_data.LocalHeight();
  size_t wlocal = sub_data.LocalWidth();

  size_t rowoff = sub_data.ColShift();
  size_t rowstride = sub_data.ColStride();
  size_t row;

  size_t global_spec;
  size_t global_lambda;

  size_t out_spec;
  size_t out_lambda;

  size_t in_spec;
  size_t in_lambda;

  double val;

  cerr << "mpi_sub proc " << sub_data.Grid().Rank() << " has local storage " << hlocal << " x " << wlocal << endl;

  if ( wlocal > 0 ) {
    for ( size_t j = 0; j < hlocal; ++j ) {
      // the global element of sub_data
      row = rowoff + j * rowstride;

      // compute the location in the full sub region

      out_spec = (size_t)( row / sub_region.n_lambda );
      out_lambda = row - ( out_spec * sub_region.n_lambda );

      global_spec = sub_region.first_spec + out_spec;
      global_lambda = sub_region.first_lambda + out_lambda;

      // are we within the selected region?

      if ( ( ( global_spec >= sub_firstspec ) && ( global_spec < sub_firstspec + sub_nspec ) ) && ( ( global_lambda >= sub_firstlambda ) && ( global_lambda < sub_firstlambda + sub_nlambda ) ) ) {

        // compute the location in the full region

        in_lambda = global_lambda - full_region.first_lambda;
        in_spec = global_spec - full_region.first_spec;

        // copy

        val = full_loc.Get ( in_spec * full_region.n_lambda + in_lambda, 0 );

        cerr << "proc " << sub_data.Grid().Rank() << " setting sub element " << row << endl;

        sub_data.SetLocal ( j, 0, val );

      }

    }
  }

  return;
}


void harp::mpi_accum_spec ( spec_slice_region const & sub_region, spec_slice_region const & full_region, mpi_matrix const & sub_data, bool use_good_sub, mpi_matrix & full_data ) {

  // verify that data dimensions match the sizes of the regions.  Select whether we
  // are using the full or good extent of the input. 

  if ( (full_region.n_spec * full_region.n_lambda) != full_data.Height() ) {
    std::ostringstream o;
    o << "output region total bins (" << (full_region.n_spec * full_region.n_lambda) << ") does not match output data size (" << full_data.Height() << ")";
    HARP_THROW( o.str().c_str() );
  }

  if ( (sub_region.n_spec * sub_region.n_lambda) != sub_data.Height() ) {
    std::ostringstream o;
    o << "input region total bins (" << (sub_region.n_spec * sub_region.n_lambda) << ") does not match input data size (" << sub_data.Height() << ")";
    HARP_THROW( o.str().c_str() );
  }

  size_t nsub;

  size_t sub_nlambda;
  size_t sub_nspec;
  size_t sub_firstspec;
  size_t sub_firstlambda;

  if ( use_good_sub ) {
    nsub = sub_region.n_good_spec * sub_region.n_good_lambda;
    sub_nlambda = sub_region.n_good_lambda;
    sub_nspec = sub_region.n_good_spec;
    sub_firstspec = sub_region.first_good_spec;
    sub_firstlambda = sub_region.first_good_lambda;
  } else {
    nsub = sub_region.n_spec * sub_region.n_lambda;
    sub_nlambda = sub_region.n_lambda;
    sub_nspec = sub_region.n_spec;
    sub_firstspec = sub_region.first_spec;
    sub_firstlambda = sub_region.first_lambda;
  }

  if ( sub_firstspec < full_region.first_spec ) {
    std::ostringstream o;
    o << "sub region first spec (" << sub_firstspec << ") is before first spec of full region (" << full_region.first_spec << ")";
    HARP_THROW( o.str().c_str() );
  }

  if ( sub_firstspec + sub_nspec > full_region.first_spec + full_region.n_spec ) {
    std::ostringstream o;
    o << "sub region last spec (" << (sub_firstspec + sub_nspec - 1) << ") is beyond last spec of full region (" << (full_region.first_spec + full_region.n_spec - 1) << ")";
    HARP_THROW( o.str().c_str() );
  }

  if ( sub_firstlambda < full_region.first_lambda ) {
    std::ostringstream o;
    o << "sub region first lambda (" << sub_firstlambda << ") is before first lambda of full region (" << full_region.first_lambda << ")";
    HARP_THROW( o.str().c_str() );
  }

  if ( sub_firstlambda + sub_nlambda > full_region.first_lambda + full_region.n_lambda ) {
    std::ostringstream o;
    o << "sub region last lambda (" << (sub_firstlambda + sub_nlambda - 1) << ") is beyond last lambda of full region (" << (full_region.first_lambda + full_region.n_lambda - 1) << ")";
    HARP_THROW( o.str().c_str() );
  }

  // FIXME: can this be changed to not store a full local copy of the input?

  elem_matrix_local full_loc ( full_data.Height(), 1 );
  local_matrix_zero ( full_loc );

  // Update output data with proper slices from the input

  size_t global_spec;
  size_t global_lambda;

  size_t out_spec;
  size_t out_lambda;

  size_t in_spec;
  size_t in_lambda;

  size_t hlocal = sub_data.LocalHeight();
  size_t wlocal = sub_data.LocalWidth();

  size_t rowoff = sub_data.ColShift();
  size_t rowstride = sub_data.ColStride();
  size_t row;

  double val;

  if ( wlocal > 0 ) {
    for ( size_t j = 0; j < hlocal; ++j ) {
      // the global element of sub_data
      row = rowoff + j * rowstride;

      // compute the location in the full sub region

      in_spec = (size_t)( row / sub_region.n_lambda );
      in_lambda = row - ( in_spec * sub_region.n_lambda );

      global_spec = sub_region.first_spec + in_spec;
      global_lambda = sub_region.first_lambda + in_lambda;

      // are we within the selected region?

      if ( ( ( global_spec >= sub_firstspec ) && ( global_spec < sub_firstspec + sub_nspec ) ) && ( ( global_lambda >= sub_firstlambda ) && ( global_lambda < sub_firstlambda + sub_nlambda ) ) ) {

        // compute the location in the full region

        out_lambda = global_lambda - full_region.first_lambda;
        out_spec = global_spec - full_region.first_spec;

        // copy

        val = sub_data.GetLocal ( j, 0 );

        full_loc.Set ( out_spec * full_region.n_lambda + out_lambda, 0, val );

      }

    }
  }

  // accumulate

  elem::AxpyInterface < double > locglob;
  locglob.Attach( elem::LOCAL_TO_GLOBAL, full_data );
  locglob.Axpy ( 1.0, full_loc, 0, 0 );
  locglob.Detach();

  return;
}




/*

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

  matrix_local local_inv ( local_rows, local_rows );
  local_matrix_zero ( local_inv );

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

  elem::AxpyInterface < double > locglob;
  locglob.Attach( elem::LOCAL_TO_GLOBAL, invcov );
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

    int send_size_key = (shift * 2 * np) + 2 * myp;
    int send_data_key = (shift * 2 * np) + 2 * myp + 1;
    int recv_size_key = (shift * 2 * np) + 2 * from_proc;
    int recv_data_key = (shift * 2 * np) + 2 * from_proc + 1;

    // send our data, then wait receive the next block

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

    size_t axpy_row = 0;
    size_t axpy_col = 0;
    bool participate;

    if ( from_proc > myp ) {
      // we are computing the transposed block of the output

      participate = true;

      local_inv.ResizeTo ( other_block->local_rows, local_rows );
      local_matrix_zero ( local_inv );

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

    } else if ( ( np % 2 != 0 ) || ( shift != nshift -1 ) ) {
      // always compute non-transposed block if we have odd number of processes
      // or we have an even number of process and we are not on the last shift.

      participate = true;

      local_inv.ResizeTo ( local_rows, other_block->local_rows );
      local_matrix_zero ( local_inv );

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

    } else {
      // for even numbers of processes, on the last pass, half the processors have
      // duplicate blocks and sit idle.

      participate = false;

    }

    delete other_block;

    // Wait for everyone to finish sending and calculating their blocks before
    // entering into Axpy communication code.

    ret = MPI_Wait ( &send_size_request, &status );
    mpi_check ( psf.Comm(), ret );

    ret = MPI_Wait ( &send_request, &status );
    mpi_check ( psf.Comm(), ret );

    free ( sendbuf );

    ret = MPI_Barrier ( psf.Comm() );

    // accumulate to global matrix

    locglob.Attach( elem::LOCAL_TO_GLOBAL, invcov );
    if ( participate ) {
      locglob.Axpy ( 1.0, local_inv, axpy_row, axpy_col );
    }
    locglob.Detach();

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


void harp::extract ( matrix_dist & D, matrix_dist & W, matrix_dist & S, matrix_dist & z, matrix_dist & Rf, matrix_dist & f ) {

  dist_matrix_zero ( f );
  dist_matrix_zero ( Rf );

  // compose ( W^T D^{-1/2} W )

  matrix_dist invrtC ( W );
  dist_matrix_zero ( invrtC );

  eigen_compose ( EIG_INVSQRT, D, W, invrtC );

  // Compute R * C

  matrix_dist RC ( invrtC );

  apply_norm ( S, RC );

  // Compute R * f.

  elem::Gemv ( elem::NORMAL, 1.0, RC, z, 0.0, Rf );

  // compute deconvolved spectra (numerically unstable, but useful for visualization).
  // R^-1 == ( W^T D^{-1/2} W ) S

  matrix_dist temp ( Rf );

  apply_inverse_norm ( S, temp );

  // now apply invrtC.  This destroys upper triangle of invrtC!

  elem::Symv ( elem::LOWER, 1.0, invrtC, temp, 0.0, f );

  return;
}

*/

