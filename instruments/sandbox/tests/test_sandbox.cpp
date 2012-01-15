#include <iostream>
#include <fstream>

#include <harp_test.hpp>

#include <boost/random.hpp>

#include <moat.hpp>

extern "C" {
#include <unistd.h>
}

using namespace std;
using namespace harp;





void sandbox_profile ( string const & name, string const & desc, double & totaltime, double & opencltime, map < string, long long int > & papi ) {
  
  cerr << "Profiling:   " << desc << ":  " << totaltime << " seconds" << endl;
  
  return;
}


void harp::test_sandbox ( string const & datadir ) {
  
  cerr << "Testing sandbox image construction..." << endl;
  
  boost::property_tree::ptree props;
  props.put ( "path", datadir + "/test_medium_input_image.fits" );
  props.put ( "signal", 1 );
  props.put ( "noise", 2 );
  image_p testimg ( image::create ( string("sandbox"), props ) );
  
  cerr << "  (PASSED)" << endl;
  
  
  cerr << "Testing sandbox image read..." << endl;
  
  mat_denserow pix ( testimg->rows(), testimg->cols() );
  mat_denserow noise ( testimg->rows(), testimg->cols() );
  
  testimg->read ( 0, 0, pix );
  testimg->read_noise ( 0, 0, noise );
  
  cerr << "  (PASSED)" << endl;
  
  
  cerr << "Testing sandbox spectrum construction..." << endl;
  
  props.clear();
  props.put ( "path", datadir + "/test_medium_input_spectra.fits" );
  props.put ( "hdu", 1 );
  spec_p testspec ( spec::create ( string("sandbox"), props ) );
  
  cerr << "  (PASSED)" << endl;
  
  
  cerr << "Testing sandbox spectrum read..." << endl;
  
  vec_dense spec ( testspec->spectrum_size( 3 ) );
  
  testspec->read_spectrum ( 3, spec );

  //for ( size_t i = 0; i < spec.size(); ++i ) {
  //  cout << spec[i] << endl;
  //}
  
  cerr << "  (PASSED)" << endl;
  
  
  cerr << "Testing sandbox PSF construction..." << endl;
  
  props.clear();
  props.put ( "path", datadir + "/test_medium_psf.fits" );
  props.put ( "corr", 10 );
  psf_p testpsf ( psf::create ( string("sandbox"), props ) );
  
  cerr << "  found " << testpsf->nspec() << " spectra" << endl;
  
  cerr << "  each with " << testpsf->specsize(0) << " flux bins" << endl;
  
  cerr << "  (PASSED)" << endl;
  
  
  cerr << "Testing sandbox PSF lambda read..." << endl;
  
  vec_dense lambda;
  
  for ( size_t i = 0; i < testpsf->nspec(); ++i ) {
    testpsf->lambda ( i, lambda );
  }
  
  cerr << "  (PASSED)" << endl;
  
  
  cerr << "Testing sandbox PSF calculation of sparse projection matrix..." << endl;
  
  moat::profile * prof = moat::profile::get ( );

  prof->reg ( "SANDBOX_PSF", "compute projection matrix" );
  prof->reg ( "SANDBOX_REMAP", "remap projection matrix" );
  
  size_t nbins = testpsf->nspec() * testpsf->specsize(0);
  size_t npix = testimg->rows() * testimg->cols();
  
  mat_compcol projmat ( npix, nbins );

  testpsf->projection ( string("SANDBOX_PSF"), string("SANDBOX_REMAP"), 0, testpsf->nspec() - 1, 0, testpsf->specsize(0) - 1, (size_t)0, testimg->cols() - 1, (size_t)0, testimg->rows() - 1, projmat );
  
  cerr << "  (PASSED)" << endl;

  
  cerr << "Testing sandbox solve..." << endl << endl;
  
  vec_dense inspec ( nbins );
  vec_dense outspec ( nbins );
  vec_flag flags ( nbins );
  
  vec_dense measured ( npix );
  vec_dense imgnoise ( npix );

  
  // read input spectrum

  testspec->read ( inspec );

  // read image data and noise

  testimg->read ( measured );
  testimg->read_noise ( imgnoise );
  
  for ( size_t b = 0; b < nbins; ++b ) {
    flags[b] = 0;
    /*
    //if ( ( b % 60 < 5 ) || ( b % 60 > 55 ) ) {
    //  flags[b] = 1;
    //}
    
    if ( ( b % 5 == 0 ) && ( b % testpsf->specsize(0) != 0 ) ) {
      inspec[b] = 2000.0;
    } else {
      inspec[b] = 0.0;
    }
    */
    outspec[b] = 0.0;
  }

  // write input spectra

  cerr << "Write input check spectra..." << endl;

  std::ostringstream o;

  string delpath = datadir + "/test_medium_input_spectra_check.fits.out";
  int ret = unlink ( delpath.c_str() );

  props.clear();
  props.put ( "hdu", 1 );
  props.put ( "nspec", testpsf->nspec() );
  props.put ( "specsize", testpsf->specsize(0) );

  spec_p specinput ( spec::create ( string("sandbox"), props ) );

  specinput->write ( datadir + "/test_medium_input_spectra_check.fits.out", inspec );

  cerr << "  (PASSED)" << endl;

  // construct signal image

  vec_dense imgcheck ( npix );

  cerr << "Check input noisy image..." << endl;

  boost::numeric::ublas::axpy_prod ( projmat, inspec, imgcheck, true );

  // add noise to image check

  typedef boost::ecuyer1988 base_generator_type;
  base_generator_type generator(42u);
  
  vec_dense rms ( npix );
  
  for ( size_t i = 0; i < npix; ++i ) {
    //rms[i] = sqrt( 16.0 + measured[i] );
    rms[i] = 1.0 / sqrt ( imgnoise[i] );
    
    boost::normal_distribution < double > dist ( 0.0, rms[i] );
    
    boost::variate_generator < base_generator_type&, boost::normal_distribution < double > > gauss ( generator, dist );
    
    imgcheck[i] += gauss();
    
    //imgnoise[i] = 1.0 / (rms[i] * rms[i]);
  }

  cerr << "  (PASSED)" << endl;


  // write out image

  cerr << "Write input check image and noise variance..." << endl;

  delpath = datadir + "/test_medium_input_image_check.fits.out";
  ret = unlink ( delpath.c_str() );

  props.clear();
  props.put ( "signal", 1 );
  props.put ( "noise", 2 );
  props.put ( "rows", testimg->rows() );
  props.put ( "cols", testimg->cols() );

  image_p outimage ( image::create ( string("sandbox"), props ) );
  
  outimage->write ( datadir + "/test_medium_input_image_check.fits.out", imgcheck );
  outimage->write_noise ( datadir + "/test_medium_input_image_check.fits.out", imgnoise );

  cerr << "  (PASSED)" << endl;

  
  // construct inverse pixel noise covariance

  cerr << "Construct pixel noise covariance..." << endl;
  
  mat_comprow invnoise ( npix, npix, npix );

  for ( size_t i = 0; i < npix; ++i ) {
    invnoise ( i, i ) = imgnoise[i];
  }
  
  cerr << "  (PASSED)" << endl;

  // noise weighted spectra

  cerr << "Compute Z" << endl;

  vec_dense z ( nbins );

  noise_weighted_spec < mat_compcol, mat_comprow, vec_dense > ( projmat, invnoise, measured, z );

  //vec_dense ztemp ( npix );
  //boost::numeric::ublas::axpy_prod ( invnoise, measured, ztemp, true );
  //boost::numeric::ublas::axpy_prod ( ztemp, projmat, z, true ); 
  
  delpath = datadir + "/test_medium_Z_spectra.fits.out";
  ret = unlink ( delpath.c_str() );

  props.clear();
  props.put ( "hdu", 1 );
  props.put ( "nspec", testpsf->nspec() );
  props.put ( "specsize", testpsf->specsize(0) );

  spec_p zspec ( spec::create ( string("sandbox"), props ) );
  zspec->write ( datadir + "/test_medium_Z_spectra.fits.out", z );

  cerr << "  (PASSED)" << endl;


  // construct the inverse spectral covariance matrix

  mat_comprow invcov ( nbins, nbins );

  mat_dynrow builder ( nbins, nbins );

  mat_compcol temp ( npix, nbins );

  cerr << "Multiply N^-1 A..." << endl;

  moat::la::multiply_mm < mat_comprow, mat_compcol, mat_compcol > ( invnoise, projmat, temp, false, true, true, "" );

  //boost::numeric::ublas::axpy_prod ( invnoise, projmat, temp, true );

  cerr << "  (PASSED)" << endl;


  cerr << "Multiply A^T (N^-1 A)..." << endl;

  moat::la::multiply_mm < mat_compcol, mat_compcol, mat_dynrow > ( projmat, temp, builder, true, true, true, std::string("") );
  
  //boost::numeric::ublas::axpy_prod ( boost::numeric::ublas::trans ( projmat ), temp, builder, true );

  mat_dynrow::iterator2 itcol;
  mat_dynrow::iterator1 itrow;

  for ( itcol = builder.begin2(); itcol != builder.end2(); ++itcol ) {
    for ( itrow = itcol.begin(); itrow != itcol.end(); ++itrow ) {
      invcov ( itrow.index1(), itrow.index2() ) = (*itrow);
    }
  }

  delpath = datadir + "/test_medium_invC.fits.out";
  ret = unlink ( delpath.c_str() );

  props.clear();
  props.put ( "signal", 1 );
  props.put ( "noise", 2 );
  props.put ( "rows", nbins );
  props.put ( "cols", nbins );

  image_p invcimage ( image::create ( string("sandbox"), props ) );
  
  invcimage->write ( datadir + "/test_medium_invC.fits.out", imgcheck );

  cerr << "  (PASSED)" << endl;

  
  // extraction

  cerr << "Doing extraction..." << endl;

  extract_dense < mat_comprow > ( invcov, z, outspec );

  cerr << "  (PASSED)" << endl;


  // write output  

  string outdata = datadir + "/test_medium_output.out";

  fstream out;
  out.open ( outdata.c_str(), ios::out );
  
  for ( size_t i = 0; i < nbins; ++i ) {
    out << i << " " << inspec[i] << " " << outspec[i] << " " << 1.0 / sqrt( invcov(i,i) )<< endl;
  }
  
  out.close();


  cerr << "  (PASSED)" << endl;
     
  return;
}



