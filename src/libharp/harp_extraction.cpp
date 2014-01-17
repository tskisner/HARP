// @COPYRIGHT@

#include <harp_internal.hpp>


using namespace std;
using namespace harp;


void harp::sub_spec ( spec_slice_region const & full_region, bool use_good_full, spec_slice_region const & sub_region, bool use_good_sub, vector_double const & full_data, vector_double & sub_data ) {

  // verify that data input dimensions match the size of the full region.  Select whether we
  // are using the full or good extent of the input and output slices.

  size_t nfull;
  size_t nsub;

  size_t full_nlambda;
  size_t full_nspec;
  size_t full_firstspec;
  size_t full_firstlambda;

  size_t sub_nlambda;
  size_t sub_nspec;
  size_t sub_firstspec;
  size_t sub_firstlambda;

  if ( use_good_full ) {
    nfull = full_region.n_good_spec * full_region.n_good_lambda;
    full_nlambda = full_region.n_good_lambda;
    full_nspec = full_region.n_good_spec;
    full_firstspec = full_region.first_good_spec;
    full_firstlambda = full_region.first_good_lambda;
  } else {
    nfull = full_region.n_spec * full_region.n_lambda;
    full_nlambda = full_region.n_lambda;
    full_nspec = full_region.n_spec;
    full_firstspec = full_region.first_spec;
    full_firstlambda = full_region.first_lambda;
  }

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

  if ( nfull != full_data.size() ) {
    std::ostringstream o;
    o << "input region total bins (" << nfull << ") does not match input data size (" << full_data.size() << ")";
    HARP_THROW( o.str().c_str() );
  }

  if ( sub_firstspec < full_firstspec ) {
    std::ostringstream o;
    o << "sub region first spec (" << sub_firstspec << ") is before first spec of full region (" << full_firstspec << ")";
    HARP_THROW( o.str().c_str() );
  }

  if ( sub_firstspec + sub_nspec > full_firstspec + full_nspec ) {
    std::ostringstream o;
    o << "sub region last spec (" << (sub_firstspec + sub_nspec - 1) << ") is beyond last spec of full region (" << (full_firstspec + full_nspec - 1) << ")";
    HARP_THROW( o.str().c_str() );
  }

  if ( sub_firstlambda < full_firstlambda ) {
    std::ostringstream o;
    o << "sub region first lambda (" << sub_firstlambda << ") is before first lambda of full region (" << full_firstlambda << ")";
    HARP_THROW( o.str().c_str() );
  }

  if ( sub_firstlambda + sub_nlambda > full_firstlambda + full_nlambda ) {
    std::ostringstream o;
    o << "sub region last lambda (" << (sub_firstlambda + sub_nlambda - 1) << ") is beyond last lambda of full region (" << (full_firstlambda + full_nlambda - 1) << ")";
    HARP_THROW( o.str().c_str() );
  }

  sub_data.resize ( nsub );
  sub_data.clear();

  // Update output data with proper slices from the input

  size_t out_spec;
  size_t out_lambda;
  size_t in_spec;
  size_t in_lambda;

  for ( size_t i = 0; i < nsub; ++i ) {
    out_spec = sub_firstspec + (size_t)( i / sub_nlambda );
    out_lambda = sub_firstlambda + ( i - out_spec * sub_nlambda );
    in_lambda = out_lambda - full_firstlambda;
    in_spec = out_spec - full_firstspec;
    sub_data[i] = full_data [ in_spec * full_nlambda + in_lambda ];
  }

  return;
}


void harp::accum_spec ( spec_slice_region const & sub_region, bool use_good_sub, spec_slice_region const & full_region, bool use_good_full, vector_double const & sub_data, vector_double & full_data ) {

  // verify that data input dimensions match the size of the sub region.  Select whether we
  // are using the full or good extent of the input and output slices.

  size_t nfull;
  size_t nsub;

  size_t full_nlambda;
  size_t full_nspec;
  size_t full_firstspec;
  size_t full_firstlambda;

  size_t sub_nlambda;
  size_t sub_nspec;
  size_t sub_firstspec;
  size_t sub_firstlambda;

  if ( use_good_full ) {
    nfull = full_region.n_good_spec * full_region.n_good_lambda;
    full_nlambda = full_region.n_good_lambda;
    full_nspec = full_region.n_good_spec;
    full_firstspec = full_region.first_good_spec;
    full_firstlambda = full_region.first_good_lambda;
  } else {
    nfull = full_region.n_spec * full_region.n_lambda;
    full_nlambda = full_region.n_lambda;
    full_nspec = full_region.n_spec;
    full_firstspec = full_region.first_spec;
    full_firstlambda = full_region.first_lambda;
  }

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

  if ( nsub != sub_data.size() ) {
    std::ostringstream o;
    o << "input region total bins (" << nsub << ") does not match input data size (" << sub_data.size() << ")";
    HARP_THROW( o.str().c_str() );
  }

  if ( nfull != full_data.size() ) {
    std::ostringstream o;
    o << "output region total bins (" << nfull << ") does not match output data size (" << full_data.size() << ")";
    HARP_THROW( o.str().c_str() );
  }

  if ( sub_firstspec < full_firstspec ) {
    std::ostringstream o;
    o << "sub region first spec (" << sub_firstspec << ") is before first spec of full region (" << full_firstspec << ")";
    HARP_THROW( o.str().c_str() );
  }

  if ( sub_firstspec + sub_nspec > full_firstspec + full_nspec ) {
    std::ostringstream o;
    o << "sub region last spec (" << (sub_firstspec + sub_nspec - 1) << ") is beyond last spec of full region (" << (full_firstspec + full_nspec - 1) << ")";
    HARP_THROW( o.str().c_str() );
  }

  if ( sub_firstlambda < full_firstlambda ) {
    std::ostringstream o;
    o << "sub region first lambda (" << sub_firstlambda << ") is before first lambda of full region (" << full_firstlambda << ")";
    HARP_THROW( o.str().c_str() );
  }

  if ( sub_firstlambda + sub_nlambda > full_firstlambda + full_nlambda ) {
    std::ostringstream o;
    o << "sub region last lambda (" << (sub_firstlambda + sub_nlambda - 1) << ") is beyond last lambda of full region (" << (full_firstlambda + full_nlambda - 1) << ")";
    HARP_THROW( o.str().c_str() );
  }

  // Update output data with proper slices from the input

  size_t out_spec;
  size_t out_lambda;
  size_t in_spec;
  size_t in_lambda;

  for ( size_t i = 0; i < nsub; ++i ) {
    in_spec = sub_firstspec + (size_t)( i / sub_nlambda );
    in_lambda = sub_firstlambda + ( i - out_spec * sub_nlambda );
    out_lambda = out_lambda - full_firstlambda;
    out_spec = out_spec - full_firstspec;
    full_data [ out_spec * full_nlambda + out_lambda ] += sub_data[i];
  }

  return;
}


void harp::spec_project ( matrix_double_sparse const & AT, vector_double const & spectra, vector_double & img ) {

  // FIXME:  for now, we just use the (unthreaded) boost sparse matrix-vector product.  If this
  // operation dominates the cost in any way, we can add a threaded implementation here.

  size_t nbins = AT.size1();
  size_t npix = AT.size2();

  if ( spectra.size() != nbins ) {
    std::ostringstream o;
    o << "number of spectral input bins (" << spectra.size() << ") does not match number of bins in PSF (" << nbins << ")";
    HARP_THROW( o.str().c_str() );
  }

  img.resize ( npix );

  boost::numeric::ublas::axpy_prod ( spectra, AT, img, true );

  return;
}


void harp::noise_weighted_spec ( matrix_double_sparse const & AT, vector_double const & invnoise, vector_double const & img, vector_double & z ) {

  size_t nbins = AT.size1();
  size_t npix = AT.size2();

  if ( invnoise.size() != npix ) {
    std::ostringstream o;
    o << "number of rows in inverse noise covariance (" << invnoise.size() << ") does not match number of pixels in A^T (" << npix << ")";
    HARP_THROW( o.str().c_str() );
  }

  if ( img.size() != npix ) {
    std::ostringstream o;
    o << "number of elements in image vector (" << img.size() << ") does not match number of pixels in A^T (" << npix << ")";
    HARP_THROW( o.str().c_str() );
  }

  z.resize ( nbins );

  // apply noise covariance to image

  vector_double imgnse ( npix );
  for ( size_t i = 0; i < npix; ++i ) {
    imgnse[i] = invnoise[i] * img[i];
  }

  // accumulate z

  // FIXME:  for now, we just use the (unthreaded) boost sparse matrix-vector product.  If this
  // operation dominates the cost in any way, we can add a threaded implementation here.

  boost::numeric::ublas::axpy_prod ( AT, imgnse, z, true );

  return;
}


void harp::inverse_covariance ( matrix_double_sparse const & AT, vector_double const & invnoise, matrix_double & invC ) {

  // check consistent sizes

  size_t nbins = AT.size1();
  size_t npix = AT.size2();

  if ( invnoise.size() != npix ) {
    std::ostringstream o;
    o << "number of rows in inverse noise covariance (" << invnoise.size() << ") does not match number of pixels in A^T (" << npix << ")";
    HARP_THROW( o.str().c_str() );
  }

  invC.resize ( nbins, nbins );
  invC.clear();

  // Accumulate output C^-1 = ( A^T N^-1 A ).  We have two sets of iterators over the sparse elements of A^T,
  // called "left" and "right" in analogy with the order of the terms in the equation.
  // This allows us to accumulate the (dense) elements of invC without ever storing a dense copy of A^T.
  // With some extra logic, we could accumulate only the upper or lower triangle and then copy the symmetric
  // elements to the other triangle.  This is not done in the serial code here.  One could also do thread-level
  // parallel accumulation of the output, in analogy with the MPI code.  This is not done here.

  matrix_double_sparse :: const_iterator1 left_rowit = AT.begin1();
  matrix_double_sparse :: const_iterator1 right_rowit;

  matrix_double_sparse :: const_iterator2 left_colit;
  matrix_double_sparse :: const_iterator2 right_colit;

  for ( left_rowit = AT.begin1(); left_rowit != AT.end1(); ++left_rowit ) {

    for ( right_rowit = AT.begin1(); right_rowit != AT.end1(); ++right_rowit ) {

      right_colit = right_rowit.begin();

      double val = 0.0;

      for ( left_colit = left_rowit.begin(); left_colit != left_rowit.end(); ++left_colit ) {

        while ( ( right_colit != right_rowit.end() ) && ( right_colit.index2() < left_colit.index2() ) ) {
          ++right_colit;
        }

        if ( right_colit.index2() == left_colit.index2() ) {
          val += invnoise[ left_colit.index2() ] * (*right_colit) * (*left_colit);
        }

      }

      invC ( left_rowit.index1(), right_rowit.index1() ) = val;

    }

  }

  return;
}


void harp::resolution ( vector_double const & D, matrix_double const & W, vector_double & S, matrix_double & R ) {

  R.resize ( W.size1(), W.size2() );
  R.clear();

  eigen_compose ( EIG_SQRT, D, W, R );

  column_norm ( R, S );

  apply_norm ( S, R );

  return;
}


void harp::extract ( vector_double const & D, matrix_double const & W, vector_double const & S, vector_double const & z, vector_double & Rf, vector_double & f ) {

  Rf.resize ( z.size() );
  f.resize ( z.size() );

  Rf.clear();
  f.clear();

  // compose ( W^T D^{-1/2} W )

  matrix_double invrtC ( W.size1(), W.size2() );
  invrtC.clear();

  eigen_compose ( EIG_INVSQRT, D, W, invrtC );

  // Compute R * C

  matrix_double RC ( invrtC );

  apply_norm ( S, RC );

  // Compute R * f.

  boost::numeric::bindings::blas::gemv ( 1.0, RC, z, 0.0, Rf );

  // compute deconvolved spectra (numerically unstable, but useful for visualization).
  // R^-1 == ( W^T D^{-1/2} W ) S

  vector_double temp ( Rf );

  apply_inverse_norm ( S, temp );

  // now apply invrtC.  This destroys upper triangle of invrtC!

  boost::numeric::bindings::blas::symv ( 1.0, boost::numeric::bindings::lower( invrtC ), temp, 0.0, f );

  return;
}


void harp::extract_slices ( spec_slice_p slice, psf_p design, vector_double const & img, vector_double const & img_inv_var, vector_double const & truth, vector_double & Rf, vector_double & f, vector_double & err, vector_double & Rtruth, map < string, double > & profile, bool region_threads, string const & status_prefix ) {

  // check dimensions and clear output

  if ( img.size() != img_inv_var.size() ) {
    HARP_THROW( "image and inverse pixel variance must have the same size" );
  }

  size_t psf_imgsize = design->img_rows() * design->img_cols();

  if ( img.size() != psf_imgsize ) {
    HARP_THROW( "image size must match PSF image dimensions" );
  }

  size_t psf_specsize = design->n_spec() * design->n_lambda();

  if ( ( truth.size() > 0 ) && ( truth.size() != psf_specsize ) ) {
    HARP_THROW( "input truth size must either be zero or match PSF spectral dimensions" );
  }

  Rf.resize ( psf_specsize );
  Rf.clear();

  f.resize ( psf_specsize );
  f.clear();

  err.resize ( psf_specsize );
  err.clear();

  Rtruth.resize ( psf_specsize );
  Rtruth.clear();

  profile.clear();

  // timing variables

  profile [ "design" ] = 0.0;
  profile [ "inverse" ] = 0.0;
  profile [ "eigen" ] = 0.0;
  profile [ "norm" ] = 0.0;
  profile [ "nsespec" ] = 0.0;
  profile [ "extract" ] = 0.0;

  // create a region description that contains all spectral bins

  spec_slice_region full_region = slice->full_region();

  // Matrices and vectors that we use for all slices.  These are resized for each slice,
  // but if the slices have the same number of bins then the resize does nothing and
  // we save doing the allocation / free each time through the loop.

  vector_double slice_truth;
  vector_double slice_Rtruth;
  vector_double slice_Rf;
  vector_double slice_f;
  vector_double slice_err;

  matrix_double invC;
  matrix_double res;
  vector_double eig_vals;
  matrix_double eig_vecs;
  vector_double z_spec;

  // Process all spectral slices

  // FIXME: implement option to thread this loop.

  size_t region_index = 0;

  for ( vector < spec_slice_region > :: iterator regit = slice->regions(0).begin(); regit != slice->regions(0).end(); ++regit ) {

    double tstart = wtime();

    size_t nbins = regit->n_spec * regit->n_lambda;

    // extract sub-data for this slice out of the global input and output
    // spectral domain products.  the sub_spec command below includes the
    // overlap for each region.  later when accumulating, we only accumulate
    // the "good" portion of each region.

    if ( truth.size() > 0 ) {
      slice_truth.resize ( nbins );
      slice_Rtruth.resize ( nbins );

      sub_spec ( full_region, false, (*regit), false, truth, slice_truth );
    }
    slice_Rf.resize ( nbins );
    slice_f.resize ( nbins );
    slice_err.resize ( nbins );

    // build the list of spectral points we want for the projection

    double tsubstart = wtime();

    map < size_t, set < size_t > > speclambda;

    for ( size_t s = 0; s < regit->n_spec; ++s ) {
      for ( size_t l = 0; l < regit->n_lambda; ++l ) {
        speclambda[ s + regit->first_spec ].insert ( l + regit->first_lambda );
      }
    }

    // get the projection matrix for this slice

    matrix_double_sparse AT;

    design->project_transpose ( speclambda, AT );

    double tsubstop = wtime();
    double time_design = ( tsubstop - tsubstart );

    // build the inverse spectral covariance for this slice

    tsubstart = wtime();

    invC.resize ( nbins, nbins );

    inverse_covariance ( AT, img_inv_var, invC );

    tsubstop = wtime();
    double time_inverse = ( tsubstop - tsubstart );

    // eigendecompose

    tsubstart = wtime();

    eig_vals.resize ( nbins );
    eig_vecs.resize ( nbins, nbins );

    eigen_decompose ( invC, eig_vals, eig_vecs );

    tsubstop = wtime();
    double time_eigen = ( tsubstop - tsubstart );

    // if we are convolving input truth spectra, then we need to explicitly compute the resolution
    // matrix, and so we do that while computing the column norm that we need.  If we are not
    // processing truth spectra, we can just compute the norm.

    tsubstart = wtime();

    if ( truth.size() > 0 ) {

      res.resize ( nbins, nbins );

      resolution ( eig_vals, eig_vecs, slice_err, res );

      boost::numeric::bindings::blas::gemv ( 1.0, res, slice_truth, 0.0, slice_Rtruth );

      accum_spec ( (*regit), true, full_region, false, slice_Rtruth, Rtruth );

    } else {

      norm ( eig_vals, eig_vecs, slice_err );

    }

    tsubstop = wtime();
    double time_norm = ( tsubstop - tsubstart );

    // compute the "noise weighted spectra", which is the RHS of the extraction
    // equation, A^T N^-1 p

    tsubstart = wtime();

    z_spec.resize ( nbins );

    noise_weighted_spec ( AT, img_inv_var, img, z_spec );

    tsubstop = wtime();
    double time_nsespec = ( tsubstop - tsubstart );

    // extract the spectra for this region

    tsubstart = wtime();

    extract ( eig_vals, eig_vecs, slice_err, z_spec, slice_Rf, slice_f );

    tsubstop = wtime();
    double time_extract = ( tsubstop - tsubstart );

    // accumulate results to global solution

    accum_spec ( (*regit), true, full_region, false, slice_Rf, Rf );
    accum_spec ( (*regit), true, full_region, false, slice_f, f );
    accum_spec ( (*regit), true, full_region, false, slice_err, err );

    // accumulate timing for this region

    double tstop = wtime();
    double time_chunk = tstop - tstart;

    profile [ "design" ] += time_design;
    profile [ "inverse" ] += time_inverse;
    profile [ "eigen" ] += time_eigen;
    profile [ "norm" ] += time_norm;
    profile [ "nsespec" ] += time_nsespec;
    profile [ "extract" ] += time_extract;

    // optionally write progress

    if ( status_prefix != "" ) {

      cout << status_prefix << " finished chunk " << region_index << endl;
      cout << status_prefix << "   computing A^T = " << time_design << " seconds" << endl;
      cout << status_prefix << "   building inverse covariance = " << time_inverse << " seconds" << endl;
      cout << status_prefix << "   eigendecompose inverse covariance = " << time_eigen << " seconds" << endl;
      if ( truth.size() > 0 ) {
        cout << status_prefix << "   compute column norm and resolution convolved truth = " << time_norm << " seconds" << endl;
      } else {
        cout << status_prefix << "   compute column norm = " << time_norm << " seconds" << endl;
      }
      cout << status_prefix << "   compute noise weighted spec = " << time_nsespec << " seconds" << endl;
      cout << status_prefix << "   extraction = " << time_extract << " seconds" << endl;
      cout << status_prefix << "   total chunk time = " << time_chunk << " seconds" << endl;

    }

    ++region_index;

  }

  return;
}





