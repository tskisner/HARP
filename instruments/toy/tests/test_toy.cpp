#include <iostream>

#include <fstream>

#include <harp_test.hpp>

#include <boost/random.hpp>

using namespace std;
using namespace harp;


void toy_pcgmle_prec ( data_vec_view & in, data_vec_view & out, int_vec_view & flags, void * data ) {
  data_vec * prec = (data_vec *) data;
  
  data_vec_view :: const_iterator vit;
  for ( vit = in.begin(); vit != in.end(); ++vit ) {
    size_t pos = vit.index();
    if ( flags[ pos ] == 0 ) {
      out[ pos ] = (*prec)[pos] * (*vit);
    } else {
      out[ pos ] = 0.0;
    }
  }
  
  return;
}


void toy_pcgmle_report ( double const & norm, double const & deltazero, int const & iter, double const & alpha, double const & beta, double const & delta, double const & epsilon ) {
  
  double relerr = sqrt ( delta / deltazero );
  
  cerr << "  PCG iter " << iter << ": alpha = " << alpha << " beta = " << beta << " delta = " << delta << " epsilon = " << epsilon << " rel err = " << relerr << endl;
  
  return;
}


void toy_pcgmle_profile ( string const & name, string const & desc, double & totaltime, double & opencltime, map < string, long long int > & papi ) {
  
  cerr << "Profiling:   " << desc << ":  " << totaltime << " seconds" << endl;
  
  return;
}


void harp::test_toy ( string const & datadir ) {
  
  cerr << "Testing toy image construction..." << endl;
  
  std::map < std::string, std::string > params;
  
  params[ "path" ] = datadir + "/test_image.fits";
  params[ "hdu" ] = "1";
  
  image_p testpix ( image::create ( string("toy"), params ) );
  
  params[ "hdu" ] = "2";
  
  image_p testnoise ( image::create ( string("toy"), params ) );
  
  cerr << "  (PASSED)" << endl;
  
  
  cerr << "Testing toy image read..." << endl;
  
  dense_mat pix ( testpix->rows(), testpix->cols() );
  dense_mat noise ( testpix->rows(), testpix->cols() );
  
  dense_mat_view pixview ( pix, mv_range ( 0, pix.size1() ), mv_range ( 0, pix.size2() ) );
  dense_mat_view noiseview ( noise, mv_range ( 0, noise.size1() ), mv_range ( 0, noise.size2() ) );
  
  testpix->read ( 0, 0, pixview );
  testnoise->read ( 0, 0, noiseview );
  
  cerr << "  (PASSED)" << endl;
  
  
  cerr << "Testing toy spectrum construction..." << endl;
  
  params.clear();
  
  params[ "path" ] = datadir + "/test_spectra.fits";
  params[ "hdu" ] = "1";
  params[ "pos" ] = "3";
  
  spectrum_p testspec ( spectrum::create ( string("toy"), params ) );
  
  cerr << "  (PASSED)" << endl;
  
  
  cerr << "Testing toy spectrum read..." << endl;
  
  data_vec spec ( testspec->size() );
  data_vec_view specview ( spec, mv_range ( 0, spec.size() ) );
  
  testspec->read ( specview );

  //for ( size_t i = 0; i < spec.size(); ++i ) {
  //  cout << spec[i] << endl;
  //}
  
  cerr << "  (PASSED)" << endl;
  
  
  cerr << "Testing toy PSF construction..." << endl;
  
  params.clear();
  
  params[ "path" ] = datadir + "/test_psf.fits";
  
  params[ "corr" ] = "20";
  
  psf_p testpsf ( psf::create ( string("toy"), params ) );
  
  cerr << "  found " << testpsf->nspec() << " spectra" << endl;
  
  cerr << "  each with " << testpsf->specsize(0) << " flux bins" << endl;
  
  cerr << "  (PASSED)" << endl;
  
  
  cerr << "Testing toy PSF lambda read..." << endl;
  
  data_vec lambda;
  
  for ( size_t i = 0; i < testpsf->nspec(); ++i ) {
    testpsf->lambda ( i, lambda );
  }
  
  cerr << "  (PASSED)" << endl;
  
  
  cerr << "Testing toy PSF calculation of sparse projection matrix..." << endl;
  
  moat::profile * prof = moat::profile::get ( );

  prof->reg ( "PCG_PSF", "constructing projection matrix" );
  
  size_t nbins = testpsf->nspec() * testpsf->specsize(0);
  size_t npix = testpix->rows() * testpix->cols();
  
  sparse_mat projmat ( npix, nbins );
  
  sparse_mat_view projview ( projmat, mv_range ( 0, npix ), mv_range ( 0, nbins ) );
  
  prof->start ( "PCG_PSF" );
  
  testpsf->projection ( 0, testpsf->nspec() - 1, 0, testpsf->specsize(0) - 1, (size_t)0, testpix->cols() - 1, (size_t)0, testpix->rows() - 1, projview );
  
  prof->stop ( "PCG_PSF" );
  
  cerr << "  (PASSED)" << endl;

  
  cerr << "Testing toy PCG-MLE solve..." << endl;
  
  data_vec inspec ( nbins );
  data_vec outspec ( nbins );
  int_vec flags ( nbins );
  
  data_vec measured ( npix );
  
  for ( size_t b = 0; b < nbins; ++b ) {
    flags[b] = 0;
    //if ( ( b % 60 < 5 ) || ( b % 60 > 55 ) ) {
    //  flags[b] = 1;
    //}
    
    if ( ( b % 10 == 0 ) && ( b % 60 != 0 ) ) {
      inspec[b] = 2000.0;
    } else {
      inspec[b] = 0.0;
    }
    outspec[b] = 0.0;
  }
  
  // write input spectra
  
  
  boost::numeric::ublas::axpy_prod ( projview, inspec, measured, true );
  
  dense_mat outmat ( testpix->rows(), testpix->cols() );
  
  dense_mat_view outview ( outmat, mv_range ( 0, testpix->rows() ), mv_range ( 0, testpix->cols() ) );
  
  size_t pixoff = 0;
  for ( size_t i = 0; i < testpix->rows(); ++i ) {
    for ( size_t j = 0; j < testpix->cols(); ++j ) {
      outmat( i, j ) = measured[pixoff];
      ++pixoff;
    }
  }
  
  
  params.clear();
  
  params[ "hdu" ] = "1";
  
  std::ostringstream o;
  o << testpix->rows();
  params[ "rows" ] = o.str();
  
  o.str("");
  o << testpix->cols();
  params[ "cols" ] = o.str();
  
  image_p outsigimage ( image::create ( string("toy"), params ) );
  
  outsigimage->write ( "!" + datadir + "/toy_MLE_inputs.fits.out", 0, 0, outview );
  
  
  typedef boost::ecuyer1988 base_generator_type;
  base_generator_type generator(42u);
  
  data_vec rms ( npix );
  
  for ( size_t i = 0; i < npix; ++i ) {
    rms[i] = sqrt( 16.0 + measured[i] );
    
    boost::normal_distribution < double > dist ( 0.0, rms[i] );
    
    boost::variate_generator < base_generator_type&, boost::normal_distribution < double > > gauss ( generator, dist );
    
    measured[i] += gauss();
  }
  
  
  pixoff = 0;
  for ( size_t i = 0; i < testpix->rows(); ++i ) {
    for ( size_t j = 0; j < testpix->cols(); ++j ) {
      outmat( i, j ) = measured[pixoff];
      ++pixoff;
    }
  }
  
  
  params.clear();
  
  params[ "hdu" ] = "2";
  
  o.str("");
  o << testpix->rows();
  params[ "rows" ] = o.str();
  
  o.str("");
  o << testpix->cols();
  params[ "cols" ] = o.str();
  
  image_p outsnimage ( image::create ( string("toy"), params ) );
  
  outsnimage->write ( datadir + "/toy_MLE_inputs.fits.out", 0, 0, outview );
  
  // construct inverse noise covariance and preconditioner
  
  sparse_mat invnoise ( npix, npix );
  data_vec precdata ( nbins );
  
  for ( size_t i = 0; i < npix; ++i ) {
    invnoise( i, i ) = 1.0 / ( rms[i] * rms[i] );
  }
  
  for ( size_t i = 0; i < nbins; ++i ) {
    precdata[i] = 0.0;
    for ( size_t j = 0; j < npix; ++j ) {
      precdata[i] += projmat( j, i ) * projmat( j, i ) * invnoise( j, j );
    }
  }
  
  for ( size_t i = 0; i < nbins; ++i ) {
    if ( precdata[i] > 0.0 ) {
      precdata[i] = 1.0 / precdata[i];
    } else {
      precdata[i] = 1.0;
    }
  }
  
   
  prof->reg ( "PCG_PREC", "applying preconditioner" );
  prof->reg ( "PCG_PMV", "projection matrix-vector multiply" );
  prof->reg ( "PCG_NMV", "N^-1 matrix-vector multiply" );
  prof->reg ( "PCG_VEC", "vector ops time" );
  prof->reg ( "PCG_TOT", "Total PCG time" );
  
  
  data_vec rhs ( nbins );
  data_vec q ( nbins );
  data_vec r ( nbins );
  data_vec s ( nbins );
  data_vec d ( nbins );
  
  double err = moat::la::pcg_mle < sparse_mat, sparse_mat, data_vec, int_vec > ( true, true, projmat, invnoise, measured, outspec, q, r, s, d, flags, rhs, 20, 1.0e-12, toy_pcgmle_prec, (void*)&precdata, toy_pcgmle_report, "PCG_TOT", "PCG_VEC", "PCG_PMV", "PCG_NMV", "PCG_PREC" );
  
  prof->stop_all();
  
  prof->query ( toy_pcgmle_profile );
  
  prof->unreg ( "PCG_PREC" );
  prof->unreg ( "PCG_PMV" );
  prof->unreg ( "PCG_NMV" );
  prof->unreg ( "PCG_VEC" );
  prof->unreg ( "PCG_TOT" );
  prof->unreg ( "PCG_PSF" );
  
  string outdata = datadir + "/toy_MLE_spectra.out";
  
  fstream out;
  out.open ( outdata.c_str(), ios::out );
  
  for ( size_t i = 0; i < nbins; ++i ) {
    if ( flags[i] == 0 ) {
      out << i << " " << inspec[i] << " " << outspec[i] << " " << sqrt(precdata[i]) << endl;
    }
  }
  
  out.close();
  
  cerr << "  (PASSED)" << endl;
  
  return;
}



