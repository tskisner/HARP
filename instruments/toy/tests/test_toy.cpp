#include <iostream>

#include <fstream>

#include <harp_test.hpp>

#include <boost/random.hpp>

using namespace std;
using namespace harp;




void toy_profile ( string const & name, string const & desc, double & totaltime, double & opencltime, map < string, long long int > & papi ) {
  
  cerr << "Profiling:   " << desc << ":  " << totaltime << " seconds" << endl;
  
  return;
}


void harp::test_toy ( string const & datadir ) {
  
  cerr << "Testing toy image construction..." << endl;
  
  std::map < std::string, std::string > params;
  
  params[ "path" ] = datadir + "/test_image.fits";
  params[ "signal" ] = "1";
  params[ "noise" ] = "2";
  
  image_p testimg ( image::create ( string("toy"), params ) );
  
  cerr << "  (PASSED)" << endl;
  
  
  cerr << "Testing toy image read..." << endl;
  
  mat_denserow pix ( testimg->rows(), testimg->cols() );
  mat_denserow noise ( testimg->rows(), testimg->cols() );
  
  testimg->read ( 0, 0, pix );
  testimg->read_noise ( 0, 0, noise );
  
  cerr << "  (PASSED)" << endl;
  
  
  cerr << "Testing toy spectrum construction..." << endl;
  
  params.clear();
  
  params[ "path" ] = datadir + "/test_spectra.fits";
  params[ "hdu" ] = "1";
  
  spec_p testspec ( spec::create ( string("toy"), params ) );
  
  cerr << "  (PASSED)" << endl;
  
  
  cerr << "Testing toy spectrum read..." << endl;
  
  vec_dense spec ( testspec->spectrum_size( 3 ) );
  
  testspec->read_spectrum ( 3, spec );

  //for ( size_t i = 0; i < spec.size(); ++i ) {
  //  cout << spec[i] << endl;
  //}
  
  cerr << "  (PASSED)" << endl;
  
  
  cerr << "Testing toy PSF construction..." << endl;
  
  params.clear();
  
  params[ "path" ] = datadir + "/test_psf.fits";
  
  params[ "corr" ] = "20";
  
  params[ "binning" ] = "4";
  
  psf_p testpsf ( psf::create ( string("toy"), params ) );
  
  cerr << "  found " << testpsf->nspec() << " spectra" << endl;
  
  cerr << "  each with " << testpsf->specsize(0) << " flux bins" << endl;
  
  cerr << "  (PASSED)" << endl;
  
  
  cerr << "Testing toy PSF lambda read..." << endl;
  
  vec_dense lambda;
  
  for ( size_t i = 0; i < testpsf->nspec(); ++i ) {
    testpsf->lambda ( i, lambda );
  }
  
  cerr << "  (PASSED)" << endl;
  
  
  cerr << "Testing toy PSF calculation of sparse projection matrix..." << endl;
  
  moat::profile * prof = moat::profile::get ( );

  prof->reg ( "PCG_PSF", "compute projection matrix" );
  prof->reg ( "PCG_REMAP", "remap projection matrix" );
  
  size_t nbins = testpsf->nspec() * testpsf->specsize(0);
  size_t npix = testimg->rows() * testimg->cols();
  
  mat_sparserow projmat ( npix, nbins );

  testpsf->projection ( string("PCG_PSF"), string("PCG_REMAP"), 0, testpsf->nspec() - 1, 0, testpsf->specsize(0) - 1, (size_t)0, testimg->cols() - 1, (size_t)0, testimg->rows() - 1, projmat );
  
  cerr << "  (PASSED)" << endl;

  exit(0);

  
  cerr << "Testing toy solve..." << endl;
  
  vec_dense inspec ( nbins );
  vec_dense outspec ( nbins );
  vec_flag flags ( nbins );
  
  vec_dense measured ( npix );
  
  for ( size_t b = 0; b < nbins; ++b ) {
    flags[b] = 0;
    //if ( ( b % 60 < 5 ) || ( b % 60 > 55 ) ) {
    //  flags[b] = 1;
    //}
    
    if ( ( b % 5 == 0 ) && ( b % testpsf->specsize(0) != 0 ) ) {
      inspec[b] = 1000.0;
    } else {
      inspec[b] = 0.0;
    }
    outspec[b] = 0.0;
  }
  
  // write input spectra
  
  /*
  boost::numeric::ublas::axpy_prod ( projmat, inspec, measured, true );
  
  dense_rowmat outmat ( testimg->rows(), testimg->cols() );
  
  dense_rowmat_view outview ( outmat, mv_range ( 0, testimg->rows() ), mv_range ( 0, testimg->cols() ) );
  
  size_t pixoff = 0;
  for ( size_t i = 0; i < testimg->rows(); ++i ) {
    for ( size_t j = 0; j < testimg->cols(); ++j ) {
      outmat( i, j ) = measured[pixoff];
      ++pixoff;
    }
  }
  */
  
  /*
  params.clear();
  
  params[ "signal" ] = "1";
  params[ "noise" ] = "2";
  
  std::ostringstream o;
  o << testimg->rows();
  params[ "rows" ] = o.str();
  
  o.str("");
  o << testimg->cols();
  params[ "cols" ] = o.str();
  
  image_p outsigimage ( image::create ( string("toy"), params ) );
  
  outsigimage->write ( "!" + datadir + "/toy_MLE_inputs.fits.out", 0, 0, outview );
  */
  
  /*
  typedef boost::ecuyer1988 base_generator_type;
  base_generator_type generator(42u);
  
  vec_dense rms ( npix );
  
  for ( size_t i = 0; i < npix; ++i ) {
    rms[i] = sqrt( 16.0 + measured[i] );
    
    boost::normal_distribution < double > dist ( 0.0, rms[i] );
    
    boost::variate_generator < base_generator_type&, boost::normal_distribution < double > > gauss ( generator, dist );
    
    measured[i] += gauss();
  }
  
  pixoff = 0;
  for ( size_t i = 0; i < testimg->rows(); ++i ) {
    for ( size_t j = 0; j < testimg->cols(); ++j ) {
      outmat( i, j ) = measured[pixoff];
      ++pixoff;
    }
  }
  */
  
  /*
  params.clear();
  
  params[ "signal" ] = "4";
  
  o.str("");
  o << testimg->rows();
  params[ "rows" ] = o.str();
  
  o.str("");
  o << testimg->cols();
  params[ "cols" ] = o.str();
  
  image_p outsnimage ( image::create ( string("toy"), params ) );
  
  outsnimage->write ( datadir + "/toy_MLE_inputs.fits.out", 0, 0, outview );
  */

  
  // construct inverse noise covariance and preconditioner
  
  /*
  comp_rowmat invnoise ( npix, npix );
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
  */

  /*
  
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

  */
  
  return;
}



