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

  // check that communicators are the same

  if ( (MPI_Comm)(AT.comm()) != (MPI_Comm)(z.Grid().Comm()) ) {
    std::ostringstream o;
    o << "design matrix and noise weighted spec have different communicators";
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


/*
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

