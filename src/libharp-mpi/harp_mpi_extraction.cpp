// @COPYRIGHT@

#include <harp_mpi_internal.hpp>


using namespace std;
using namespace harp;



void harp::mpi_sub_spec ( spec_slice_region const & full_region, spec_slice_region const & sub_region, mpi_matrix const & full_data, bool use_good_sub, mpi_matrix & sub_data ) {

  // verify that the process grid matches between input and output

  int np;
  int myp;

  MPI_Comm_size ( full_data.Grid().Comm(), &np );
  MPI_Comm_rank ( full_data.Grid().Comm(), &myp );

  if ( full_data.Grid() != sub_data.Grid() ) {
    std::ostringstream o;
    o << "input and output matrices have different process grids!";
    HARP_MPI_ABORT( myp, o.str().c_str() );
  }

  // verify that data dimensions match the sizes of the regions.  Select whether we
  // are using the full or good extent of the output. 

  if ( (full_region.n_spec * full_region.n_lambda) != full_data.Height() ) {
    std::ostringstream o;
    o << "input region total bins (" << (full_region.n_spec * full_region.n_lambda) << ") does not match input data size (" << full_data.Height() << ")";
    HARP_MPI_ABORT( myp, o.str().c_str() );
  }

  if ( (sub_region.n_spec * sub_region.n_lambda) != sub_data.Height() ) {
    std::ostringstream o;
    o << "output region total bins (" << (sub_region.n_spec * sub_region.n_lambda) << ") does not match output data size (" << sub_data.Height() << ")";
    HARP_MPI_ABORT( myp, o.str().c_str() );
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
    HARP_MPI_ABORT( myp, o.str().c_str() );
  }

  if ( sub_firstspec + sub_nspec > full_region.first_spec + full_region.n_spec ) {
    std::ostringstream o;
    o << "sub region last spec (" << (sub_firstspec + sub_nspec - 1) << ") is beyond last spec of full region (" << (full_region.first_spec + full_region.n_spec - 1) << ")";
    HARP_MPI_ABORT( myp, o.str().c_str() );
  }

  if ( sub_firstlambda < full_region.first_lambda ) {
    std::ostringstream o;
    o << "sub region first lambda (" << sub_firstlambda << ") is before first lambda of full region (" << full_region.first_lambda << ")";
    HARP_MPI_ABORT( myp, o.str().c_str() );
  }

  if ( sub_firstlambda + sub_nlambda > full_region.first_lambda + full_region.n_lambda ) {
    std::ostringstream o;
    o << "sub region last lambda (" << (sub_firstlambda + sub_nlambda - 1) << ") is beyond last lambda of full region (" << (full_region.first_lambda + full_region.n_lambda - 1) << ")";
    HARP_MPI_ABORT( myp, o.str().c_str() );
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

  //cerr << "mpi_sub proc " << sub_data.Grid().Rank() << " has local storage " << hlocal << " x " << wlocal << endl;

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

        //cerr << "proc " << sub_data.Grid().Rank() << " setting sub element " << row << endl;

        sub_data.SetLocal ( j, 0, val );

      }

    }
  }

  return;
}


void harp::mpi_accum_spec ( spec_slice_region const & sub_region, spec_slice_region const & full_region, mpi_matrix const & sub_data, bool use_good_sub, mpi_matrix & full_data ) {

  // verify that the process grid matches between input and output

  int np;
  int myp;

  MPI_Comm_size ( full_data.Grid().Comm(), &np );
  MPI_Comm_rank ( full_data.Grid().Comm(), &myp );

  if ( full_data.Grid() != sub_data.Grid() ) {
    std::ostringstream o;
    o << "input and output matrices have different process grids!";
    HARP_MPI_ABORT( myp, o.str().c_str() );
  }

  // verify that data dimensions match the sizes of the regions.  Select whether we
  // are using the full or good extent of the input. 

  if ( (full_region.n_spec * full_region.n_lambda) != full_data.Height() ) {
    std::ostringstream o;
    o << "output region total bins (" << (full_region.n_spec * full_region.n_lambda) << ") does not match output data size (" << full_data.Height() << ")";
    HARP_MPI_ABORT( myp, o.str().c_str() );
  }

  if ( (sub_region.n_spec * sub_region.n_lambda) != sub_data.Height() ) {
    std::ostringstream o;
    o << "input region total bins (" << (sub_region.n_spec * sub_region.n_lambda) << ") does not match input data size (" << sub_data.Height() << ")";
    HARP_MPI_ABORT( myp, o.str().c_str() );
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
    HARP_MPI_ABORT( myp, o.str().c_str() );
  }

  if ( sub_firstspec + sub_nspec > full_region.first_spec + full_region.n_spec ) {
    std::ostringstream o;
    o << "sub region last spec (" << (sub_firstspec + sub_nspec - 1) << ") is beyond last spec of full region (" << (full_region.first_spec + full_region.n_spec - 1) << ")";
    HARP_MPI_ABORT( myp, o.str().c_str() );
  }

  if ( sub_firstlambda < full_region.first_lambda ) {
    std::ostringstream o;
    o << "sub region first lambda (" << sub_firstlambda << ") is before first lambda of full region (" << full_region.first_lambda << ")";
    HARP_MPI_ABORT( myp, o.str().c_str() );
  }

  if ( sub_firstlambda + sub_nlambda > full_region.first_lambda + full_region.n_lambda ) {
    std::ostringstream o;
    o << "sub region last lambda (" << (sub_firstlambda + sub_nlambda - 1) << ") is beyond last lambda of full region (" << (full_region.first_lambda + full_region.n_lambda - 1) << ")";
    HARP_MPI_ABORT( myp, o.str().c_str() );
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


void harp::mpi_noise_weighted_spec ( mpi_matrix_sparse const & AT, elem_matrix_local const & invnoise, vector_mask const & mask, elem_matrix_local const & img, mpi_matrix & z ) {

  int myp = AT.comm().rank();
  int np = AT.comm().size();

  size_t nbins = AT.rows();
  size_t npix = AT.cols();

  if ( invnoise.Height() != npix ) {
    std::ostringstream o;
    o << "number of rows in inverse noise covariance (" << invnoise.Height() << ") does not match number of pixels in PSF (" << npix << ")";
    HARP_MPI_ABORT( myp, o.str().c_str() );
  }

  if ( img.Height() != npix ) {
    std::ostringstream o;
    o << "number of elements in image vector (" << img.Height() << ") does not match number of pixels in PSF (" << npix << ")";
    HARP_MPI_ABORT( myp, o.str().c_str() );
  }

  z.Resize ( nbins, 1 );
  mpi_matrix_zero ( z );

  // apply noise covariance to image

  elem_matrix_local weight ( npix, 1 );
  for ( size_t i = 0; i < npix; ++i ) {
    weight.Set ( i, 0, invnoise.Get(i,0) * img.Get(i,0) * (double)mask[i] );
  }

  // accumulate local pieces of z

  size_t local_firstrow = AT.block().firstrow;
  size_t local_nrows = AT.block().rows;
  size_t local_nvals = AT.block().vals;

  elem_matrix_local local_z ( local_nrows, 1 );
  local_matrix_zero ( local_z );

  double val;
  double zval;
  size_t row;
  size_t col;

  for ( size_t loc = 0; loc < local_nvals; ++loc ) {
    row = AT.block().row[ loc ];
    col = AT.block().col[ loc ];
    val = AT.block().data[ loc ];
    zval = local_z.Get ( row - local_firstrow, 0 );
    local_z.Set ( row - local_firstrow, 0, zval + val * weight.Get ( col, 0 ) );
  }

  // accumulate to output

  elem::AxpyInterface < double > locglob;
  locglob.Attach( elem::LOCAL_TO_GLOBAL, z );
  locglob.Axpy ( 1.0, local_z, local_firstrow, 0 );
  locglob.Detach();

  return;
}


void harp::mpi_inverse_covariance ( mpi_matrix_sparse const & AT, elem_matrix_local const & invnoise, vector_mask const & mask, mpi_matrix & invcov ) {

  int myp = AT.comm().rank();
  int np = AT.comm().size();

  size_t nbins = AT.rows();
  size_t npix = AT.cols();

  // check consistent sizes

  if ( invnoise.Height() != npix ) {
    std::ostringstream o;
    o << "number of rows in inverse noise covariance (" << invnoise.Height() << ") does not match number of pixels in PSF (" << npix << ")";
    HARP_MPI_ABORT( myp, o.str().c_str() );
  }

  invcov.Resize ( nbins, nbins );
  mpi_matrix_zero ( invcov );

  // First, accumulate our own diagonal subblock

  size_t local_firstrow = AT.block().firstrow;
  size_t local_nrows = AT.block().rows;
  size_t local_nvals = AT.block().vals;

  elem_matrix_local local_inv ( local_nrows, local_nrows );
  local_matrix_zero ( local_inv );

  double val;
  size_t lhs_off;
  size_t rhs_off;
  size_t lhs_nnz;
  size_t rhs_nnz;
  size_t lhs_col;
  size_t rhs_col;
  size_t j, k;

  for ( size_t lhs_row = 0; lhs_row < local_nrows; ++lhs_row ) {

    for ( size_t rhs_row = 0; rhs_row <= lhs_row; ++rhs_row ) {

      lhs_off = AT.block().row_offset[ lhs_row ];
      rhs_off = AT.block().row_offset[ rhs_row ];

      lhs_nnz = AT.block().row_nnz[ lhs_row ];
      rhs_nnz = AT.block().row_nnz[ rhs_row ];

      val = 0.0;

      k = 0;

      rhs_col = AT.block().col[ rhs_off ];

      for ( j = 0; j < lhs_nnz; ++j ) {

        lhs_col = AT.block().col[ lhs_off + j ];
        
        while ( ( rhs_col < lhs_col ) && ( k < rhs_nnz - 1 ) ) {
          ++k;
          rhs_col = AT.block().col[ rhs_off + k ];
        }

        if ( rhs_col == lhs_col ) {
          val += invnoise.Get( lhs_col, 0 ) * (double)mask[ lhs_col ] * AT.block().data[ lhs_off + j ] * AT.block().data[ rhs_off + k ];
        }

      }

      //cerr << myp << " INVC (" << lhs_row << "," << rhs_row << ") = " << val << endl;

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

  mpi_matrix_sparse_block send_block;
  mpi_matrix_sparse_block recv_block;

  for ( int shift = 0; shift < nshift; ++shift ) {

    int send_key = (shift * np) + myp;
    int recv_key = (shift * np) + from_proc;

    // send our data, then wait receive the next block

    boost::mpi::request send_req;

    if ( shift == 0 ) {
      // first shift, send our own data

      send_req = AT.comm().isend ( to_proc, send_key, AT.block() );

    } else {
      // pass along the buffer

      send_block = recv_block;

      send_req = AT.comm().isend ( to_proc, send_key, send_block );
      
    }

    // receive block from sender

    AT.comm().recv ( from_proc, recv_key, recv_block );

    // compute output block

    size_t axpy_row = 0;
    size_t axpy_col = 0;
    bool participate;

    if ( from_proc > myp ) {
      // we are computing the transposed block of the output

      participate = true;

      local_inv.Resize ( recv_block.rows, local_nrows );
      local_matrix_zero ( local_inv );

      axpy_row = recv_block.firstrow;
      axpy_col = local_firstrow;

      for ( size_t lhs_row = 0; lhs_row < local_nrows; ++lhs_row ) {

        for ( size_t rhs_row = 0; rhs_row < recv_block.rows; ++rhs_row ) {

          lhs_off = AT.block().row_offset[ lhs_row ];
          rhs_off = recv_block.row_offset[ rhs_row ];

          lhs_nnz = AT.block().row_nnz[ lhs_row ];
          rhs_nnz = recv_block.row_nnz[ rhs_row ];

          val = 0.0;

          k = 0;
          rhs_col = recv_block.col[ rhs_off ];

          for ( j = 0; j < lhs_nnz; ++j ) {

            lhs_col = AT.block().col[ lhs_off + j ];
            
            while ( ( rhs_col < lhs_col ) && ( k < rhs_nnz - 1 ) ) {
              ++k;
              rhs_col = recv_block.col[ rhs_off + k ];
            }

            if ( rhs_col == lhs_col ) {
              val += invnoise.Get( lhs_col, 0 ) * (double)mask[ lhs_col ] * AT.block().data[ lhs_off + j ] * recv_block.data[ rhs_off + k ];
            }

          }

          local_inv.Set ( rhs_row, lhs_row, val );

        }
      }

    } else if ( ( np % 2 != 0 ) || ( shift != nshift -1 ) ) {
      // always compute non-transposed block if we have odd number of processes
      // or we have an even number of process and we are not on the last shift.

      participate = true;

      local_inv.Resize ( local_nrows, recv_block.rows );
      local_matrix_zero ( local_inv );

      axpy_row = local_firstrow;
      axpy_col = recv_block.firstrow;

      for ( size_t lhs_row = 0; lhs_row < local_nrows; ++lhs_row ) {

        for ( size_t rhs_row = 0; rhs_row < recv_block.rows; ++rhs_row ) {

          lhs_off = AT.block().row_offset[ lhs_row ];
          rhs_off = recv_block.row_offset[ rhs_row ];

          lhs_nnz = AT.block().row_nnz[ lhs_row ];
          rhs_nnz = recv_block.row_nnz[ rhs_row ];

          val = 0.0;

          k = 0;
          rhs_col = recv_block.col[ rhs_off ];

          for ( j = 0; j < lhs_nnz; ++j ) {

            lhs_col = AT.block().col[ lhs_off + j ];
            
            while ( ( rhs_col < lhs_col ) && ( k < rhs_nnz - 1 ) ) {
              ++k;
              rhs_col = recv_block.col[ rhs_off + k ];
            }

            if ( rhs_col == lhs_col ) {
              val += invnoise.Get( lhs_col, 0 ) * (double)mask[ lhs_col ] * AT.block().data[ lhs_off + j ] * recv_block.data[ rhs_off + k ];
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

    // Wait for everyone to finish sending and calculating their blocks before
    // entering into Axpy communication code.

    send_req.wait();

    AT.comm().barrier();

    // accumulate to global matrix

    locglob.Attach( elem::LOCAL_TO_GLOBAL, invcov );
    if ( participate ) {
      locglob.Axpy ( 1.0, local_inv, axpy_row, axpy_col );
    }
    locglob.Detach();

  }

  return;
}


void harp::mpi_resolution ( mpi_matrix & D, mpi_matrix & W, mpi_matrix & S, mpi_matrix & R ) {

  R = W;

  mpi_matrix_zero ( R );

  mpi_eigen_compose ( EIG_SQRT, D, W, R );

  mpi_column_norm ( R, S );

  mpi_apply_norm ( S, R );

  return;
}


void harp::mpi_extract ( mpi_matrix & D, mpi_matrix & W, mpi_matrix & S, mpi_matrix & z, mpi_matrix & Rf, mpi_matrix & f ) {

  mpi_matrix_zero ( f );
  mpi_matrix_zero ( Rf );

  // compose ( W^T D^{-1/2} W )

  mpi_matrix composed ( W );
  mpi_matrix_zero ( composed );

  mpi_eigen_compose ( EIG_INVSQRT, D, W, composed );

  // Compute R * C

  mpi_apply_norm ( S, composed );

  // Compute R * f.

  elem::Gemv ( elem::NORMAL, 1.0, composed, z, 0.0, Rf );

  // compute deconvolved spectra (numerically unstable, but useful for visualization).

  mpi_eigen_compose ( EIG_INV, D, W, composed );

  elem::Gemv ( elem::NORMAL, 1.0, composed, z, 0.0, f );

  return;
}


void harp::mpi_extract_slices ( mpi_spec_slice_p slice, mpi_psf_p design, elem_matrix_local const & img, elem_matrix_local const & img_inv_var, mpi_matrix const & truth, mpi_matrix & Rf, mpi_matrix & f, mpi_matrix & err, mpi_matrix & Rtruth, std::map < std::string, double > & profile, bool lambda_mask, std::string const & status_prefix ) {

  // This function assumes that the slice has been distributed over the gangs, and that
  // the design matrix and spectral domain quantities have been distributed over the
  // global communicator (we redistribute these inside this function).

  // check dimensions and clear output

  if ( img.Height() != img_inv_var.Height() ) {
    HARP_MPI_ABORT( slice->gang_comm().rank(), "image and inverse pixel variance must have the same size" );
  }

  size_t img_rows = design->img_rows();
  size_t img_cols = design->img_cols();

  size_t psf_imgsize = img_rows * img_cols;

  if ( img.Height() != psf_imgsize ) {
    HARP_MPI_ABORT( slice->gang_comm().rank(), "image size must match PSF image dimensions" );
  }

  size_t psf_specsize = design->n_spec() * design->n_lambda();

  if ( ( truth.Height() > 0 ) && ( truth.Height() != psf_specsize ) ) {
    HARP_MPI_ABORT( slice->gang_comm().rank(), "input truth size must either be zero or match PSF spectral dimensions" );
  }

  Rf.Resize ( psf_specsize, 1 );
  mpi_matrix_zero ( Rf );

  f.Resize ( psf_specsize, 1 );
  mpi_matrix_zero ( f );

  err.Resize ( psf_specsize, 1 );
  mpi_matrix_zero ( err );

  Rtruth.Resize ( truth.Height(), 1 );
  mpi_matrix_zero ( Rtruth );

  // timing variables

  profile.clear();

  vector < string > counters;
  counters.push_back ( "design" );
  counters.push_back ( "inverse" );
  counters.push_back ( "eigen" );
  counters.push_back ( "norm" );
  counters.push_back ( "nsespec" );
  counters.push_back ( "extract" );

  for ( vector < string > :: const_iterator itcounter = counters.begin(); itcounter != counters.end(); ++itcounter ) {
    profile[ (*itcounter) ] = 0.0;
  }

  // create a region description that contains all spectral bins

  spec_slice_region full_region = slice->full_region();

  // gang-distributed quantities and redistribution

  elem::Grid gang_grid ( slice->gang_comm() );

  mpi_matrix gang_truth ( truth.Height(), 1, gang_grid );
  mpi_matrix_zero ( gang_truth );

  mpi_matrix gang_Rtruth ( truth.Height(), 1, gang_grid );
  mpi_matrix_zero ( gang_Rtruth );

  mpi_matrix gang_Rf ( psf_specsize, 1, gang_grid );
  mpi_matrix_zero ( gang_Rf );
  
  mpi_matrix gang_f ( psf_specsize, 1, gang_grid );
  mpi_matrix_zero ( gang_f );
  
  mpi_matrix gang_err ( psf_specsize, 1, gang_grid );
  mpi_matrix_zero ( gang_err );

  mpi_gang_distribute ( truth, gang_truth );

  mpi_psf_p gang_design ( design->redistribute ( slice->gang_comm() ) );

  // In order to keep the printed output clean, we use one-sided MPI calls
  // to implement a locking mechanism used by the rank zero process in all gangs.
  // unfortunately, I don't think boost::mpi supports this, so we have to call
  // the raw C API...

  int * root_lockbuf;

  MPI_Win printlock;

  if ( slice->gang_comm().rank() == 0 ) {
    int ret;
    if ( slice->rank_comm().rank() == 0 ) {
      ret = MPI_Alloc_mem ( sizeof(int), MPI_INFO_NULL, (void*)&root_lockbuf );
      ret = MPI_Win_create ( (void*)(root_lockbuf), sizeof(int), sizeof(int), MPI_INFO_NULL, slice->rank_comm(), &printlock );
    } else {
      ret = MPI_Win_create ( NULL, 0, sizeof(int), MPI_INFO_NULL, slice->rank_comm(), &printlock );
    }
  }

  // Process all spectral slices assigned to our gang

  size_t region_index = 0;

  vector < spec_slice_region > regions = slice->regions();

  for ( vector < spec_slice_region > :: const_iterator regit = regions.begin(); regit != regions.end(); ++regit ) {

    double tstart = MPI_Wtime();

    size_t nbins = regit->n_spec * regit->n_lambda;

    // slice-specific quantities distributed over the gang

    mpi_matrix slice_Rf ( nbins, 1, gang_grid );
    mpi_matrix_zero ( slice_Rf );

    mpi_matrix slice_f ( nbins, 1, gang_grid );
    mpi_matrix_zero ( slice_f );

    mpi_matrix slice_err ( nbins, 1, gang_grid );
    mpi_matrix_zero ( slice_err );

    mpi_matrix slice_truth ( 1, 1, gang_grid );
    mpi_matrix slice_Rtruth ( 1, 1, gang_grid );

    // extract sub-data for this slice out of the global input and output
    // spectral domain products.  the sub_spec command below includes the
    // overlap for each region.  later when accumulating, we only accumulate
    // the "good" portion of each region.

    if ( truth.Height() > 0 ) {
      slice_truth.Resize ( nbins, 1 );
      mpi_matrix_zero ( slice_truth );

      slice_Rtruth.Resize ( nbins, 1 );
      mpi_matrix_zero ( slice_Rtruth );

      mpi_sub_spec ( full_region, (*regit), gang_truth, false, slice_truth );
    }

    // build the list of spectral points we want for the projection.  also
    // build the pixel mask.  We are going to mask all pixels in the wavelength
    // direction which extend beyond the bin centers of the extreme bins.  In the
    // spec direction, we assume that there are either bundle gaps or that the 
    // extraction is being done across all spectra.

    double tsubstart = MPI_Wtime();

    vector_mask mask ( img_inv_var.Height() );

    // start by masking all pixels...
    mask.clear();

    size_t xoff;
    size_t yoff;
    size_t nx;
    size_t ny;

    // select our spectral bins

    map < size_t, set < size_t > > speclambda;

    for ( size_t s = 0; s < regit->n_spec; ++s ) {
      for ( size_t l = 0; l < regit->n_lambda; ++l ) {
        speclambda[ s + regit->first_spec ].insert ( l + regit->first_lambda );
      }
    }

    if ( lambda_mask ) {

      // we un-mask all pixels that touch only the good bins we are solving for

      for ( size_t s = 0; s < regit->n_good_spec; ++s ) {
        for ( size_t l = 0; l < regit->n_good_lambda; ++l ) {

          gang_design->extent ( s + regit->first_good_spec, l + regit->first_good_lambda, xoff, yoff, nx, ny );

          for ( size_t j = 0; j < nx; ++j ) {
            for ( size_t i = 0; i < ny; ++i ) {
              mask[ ( xoff + j ) * img_rows + yoff + i ] = 1;
            }
          }

        }
      }

    } else {

      // we un-mask all pixels that touch any of the bins we are solving for

      for ( size_t s = 0; s < regit->n_spec; ++s ) {
        for ( size_t l = 0; l < regit->n_lambda; ++l ) {

          gang_design->extent ( s + regit->first_spec, l + regit->first_lambda, xoff, yoff, nx, ny );

          for ( size_t j = 0; j < nx; ++j ) {
            for ( size_t i = 0; i < ny; ++i ) {
              mask[ ( xoff + j ) * img_rows + yoff + i ] = 1;
            }
          }

        }
      }

    }

    // get the projection matrix for this slice

    mpi_matrix_sparse AT ( slice->gang_comm(), nbins, img.Height() );

    gang_design->project_transpose ( speclambda, AT );

    double tsubstop = MPI_Wtime();
    double time_design = ( tsubstop - tsubstart );

    // build the inverse spectral covariance for this slice

    tsubstart = MPI_Wtime();

    mpi_matrix invC ( nbins, nbins, gang_grid );
    mpi_matrix_zero ( invC );

    mpi_inverse_covariance ( AT, img_inv_var, mask, invC );

    tsubstop = MPI_Wtime();
    double time_inverse = ( tsubstop - tsubstart );

    // eigendecompose

    tsubstart = MPI_Wtime();

    mpi_matrix eig_vals ( nbins, 1, gang_grid );
    mpi_matrix eig_vecs ( nbins, nbins, gang_grid );

    bool regularize = lambda_mask;
    mpi_eigen_decompose ( invC, eig_vals, eig_vecs, regularize );

    tsubstop = MPI_Wtime();
    double time_eigen = ( tsubstop - tsubstart );

    // if we are convolving input truth spectra, then we need to explicitly compute the resolution
    // matrix, and so we do that while computing the column norm that we need.  If we are not
    // processing truth spectra, we can just compute the norm.

    tsubstart = MPI_Wtime();

    if ( truth.Height() > 0 ) {

      mpi_matrix res ( nbins, nbins, gang_grid );

      mpi_resolution ( eig_vals, eig_vecs, slice_err, res );

      elem::Gemv ( elem::NORMAL, 1.0, res, slice_truth, 0.0, slice_Rtruth );

      mpi_accum_spec ( (*regit), full_region, slice_Rtruth, true, gang_Rtruth );

    } else {

      mpi_norm ( eig_vals, eig_vecs, slice_err );

    }

    tsubstop = MPI_Wtime();
    double time_norm = ( tsubstop - tsubstart );

    // compute the "noise weighted spectra", which is the RHS of the extraction
    // equation, A^T N^-1 p

    tsubstart = MPI_Wtime();

    mpi_matrix z_spec ( nbins, 1, gang_grid );

    mpi_noise_weighted_spec ( AT, img_inv_var, mask, img, z_spec );

    tsubstop = MPI_Wtime();
    double time_nsespec = ( tsubstop - tsubstart );

    // extract the spectra for this region

    tsubstart = MPI_Wtime();

    mpi_extract ( eig_vals, eig_vecs, slice_err, z_spec, slice_Rf, slice_f );

    tsubstop = MPI_Wtime();
    double time_extract = ( tsubstop - tsubstart );

    // accumulate results to gang solution

    mpi_accum_spec ( (*regit), full_region, slice_Rf, true, gang_Rf );
    mpi_accum_spec ( (*regit), full_region, slice_f, true, gang_f );
    mpi_accum_spec ( (*regit), full_region, slice_err, true, gang_err );

    // accumulate timing for this region

    double tstop = MPI_Wtime();
    double time_chunk = tstop - tstart;

    profile [ "design" ] += time_design;
    profile [ "inverse" ] += time_inverse;
    profile [ "eigen" ] += time_eigen;
    profile [ "norm" ] += time_norm;
    profile [ "nsespec" ] += time_nsespec;
    profile [ "extract" ] += time_extract;

    // optionally write progress

    if ( ( slice->gang_comm().rank() == 0 ) && ( status_prefix != "" ) ) {

      int ret = MPI_Win_lock ( MPI_LOCK_EXCLUSIVE, 0, 0, printlock );

      cout << status_prefix << " finished chunk " << region_index << endl;
      cout << status_prefix << "   computing A^T = " << time_design << " seconds" << endl;
      cout << status_prefix << "   building inverse covariance = " << time_inverse << " seconds" << endl;
      cout << status_prefix << "   eigendecompose inverse covariance = " << time_eigen << " seconds" << endl;
      if ( truth.Height() > 0 ) {
        cout << status_prefix << "   compute column norm and resolution convolved truth = " << time_norm << " seconds" << endl;
      } else {
        cout << status_prefix << "   compute column norm = " << time_norm << " seconds" << endl;
      }
      cout << status_prefix << "   compute noise weighted spec = " << time_nsespec << " seconds" << endl;
      cout << status_prefix << "   extraction = " << time_extract << " seconds" << endl;
      cout << status_prefix << "   total chunk time = " << time_chunk << " seconds" << endl;

      // free the lock
      ret = MPI_Win_unlock ( 0, printlock );

    }

    ++region_index;

  }

  // free the one-sided printer window

  if ( slice->gang_comm().rank() == 0 ) {
    int ret;
    MPI_Win_free ( &printlock );
    if ( slice->rank_comm().rank() == 0 ) {
      ret = MPI_Free_mem ( (void*)root_lockbuf );
    }
  }

  return;
}











