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



extern "C" void dgesdd_ ( char *, int *, int *, double *, int *, double *, double *, int *, double *, int *, double *, int *, int *, int * );




void lapack_svd ( int dim, double * mat, double * eigen ) {
  char jobz = 'O';
  int lda = dim;
  int ldu = dim;
  double * U = NULL;
  double * VT = moat::double_alloc ( dim * dim );
  int lwork;
  double * work;
  int * iwork = moat::int_alloc ( 8 * dim );
  int info;

  lwork = dim * ( 5 * dim + 7 );
  work = moat::double_alloc ( lwork );

  dgesdd_ ( &jobz, &dim, &dim, mat, &lda, eigen, U, &ldu, VT, &dim, work, &lwork, iwork, &info );

  free ( iwork );
  free ( work );
  free ( VT );

  return;
}







void toy_profile ( string const & name, string const & desc, double & totaltime, double & opencltime, map < string, long long int > & papi ) {
  
  cerr << "Profiling:   " << desc << ":  " << totaltime << " seconds" << endl;
  
  return;
}


void harp::test_toy ( string const & datadir ) {
  
  cerr << "Testing toy image construction..." << endl;
  
  std::map < std::string, std::string > params;
  params[ "path" ] = datadir + "/test_medium_input_image.fits";
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
  params[ "path" ] = datadir + "/test_medium_input_spectra.fits";
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
  params[ "path" ] = datadir + "/test_medium_psf.fits";
  params[ "corr" ] = "10";
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
  
  mat_compcol projmat ( npix, nbins );

  testpsf->projection ( string("PCG_PSF"), string("PCG_REMAP"), 0, testpsf->nspec() - 1, 0, testpsf->specsize(0) - 1, (size_t)0, testimg->cols() - 1, (size_t)0, testimg->rows() - 1, projmat );
  
  cerr << "  (PASSED)" << endl;

  
  cerr << "Testing toy solve..." << endl << endl;
  
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
  /*
  cerr << "Write input spectra..." << endl;

  std::ostringstream o;

  string delpath = datadir + "/test_input_spectra.fits.out";
  int ret = unlink ( delpath.c_str() );

  params.clear();
  params[ "hdu" ] = "1";
  o.str("");
  o << testpsf->nspec();
  params[ "nspec" ] = o.str();
  o.str("");
  o << testpsf->specsize(0);
  params[ "specsize" ] = o.str();
  spec_p specinput ( spec::create ( string("toy"), params ) );

  specinput->write ( datadir + "/test_input_spectra.fits.out", inspec );

  cerr << "  (PASSED)" << endl;

  // construct signal image

  cerr << "Construct input noisy image..." << endl;

  boost::numeric::ublas::axpy_prod ( projmat, inspec, measured, true );

  // add noise to image

  typedef boost::ecuyer1988 base_generator_type;
  base_generator_type generator(42u);
  
  vec_dense rms ( npix );
  
  for ( size_t i = 0; i < npix; ++i ) {
    //rms[i] = sqrt( 16.0 + measured[i] );
    rms[i] = sqrt( measured[i] / 16.0 );
    
    boost::normal_distribution < double > dist ( 0.0, rms[i] );
    
    boost::variate_generator < base_generator_type&, boost::normal_distribution < double > > gauss ( generator, dist );
    
    measured[i] += gauss();
    imgnoise[i] = 1.0 / (rms[i] * rms[i]);
  }

  cerr << "  (PASSED)" << endl;

  // write out image

  cerr << "Write input image and noise variance..." << endl;

  delpath = datadir + "/test_input_image.fits.out";
  ret = unlink ( delpath.c_str() );

  params.clear();

  params[ "signal" ] = "1";
  params[ "noise" ] = "2";  
  o.str("");
  o << testimg->rows();
  params[ "rows" ] = o.str();
  o.str("");
  o << testimg->cols();
  params[ "cols" ] = o.str();
  
  image_p outimage ( image::create ( string("toy"), params ) );
  
  outimage->write ( datadir + "/test_input_image.fits.out", measured );
  outimage->write_noise ( datadir + "/test_input_image.fits.out", imgnoise );

  cerr << "  (PASSED)" << endl;
  */

  
  // construct inverse pixel noise covariance

  cerr << "Construct pixel noise covariance..." << endl;
  
  mat_comprow invnoise ( npix, npix, npix );

  for ( size_t i = 0; i < npix; ++i ) {
    invnoise ( i, i ) = imgnoise[i];
  }
  
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

  boost::numeric::ublas::axpy_prod ( boost::numeric::ublas::trans ( projmat ), temp, builder, true );

  mat_dynrow::iterator2 itcol;
  mat_dynrow::iterator1 itrow;

  for ( itcol = builder.begin2(); itcol != builder.end2(); ++itcol ) {
    for ( itrow = itcol.begin(); itrow != itcol.end(); ++itrow ) {
      invcov ( itrow.index1(), itrow.index2() ) = (*itrow);
    }
  }

  cerr << "  (PASSED)" << endl;

  // noise weighted spectra

  cerr << "compute Z" << endl;

  vec_dense z ( nbins );
  vec_dense ztemp ( npix );

  boost::numeric::ublas::axpy_prod ( invnoise, measured, ztemp, true );

  boost::numeric::ublas::axpy_prod ( ztemp, projmat, z, true ); 
  
  cerr << "  (PASSED)" << endl;


  // remap to row-major?

  size_t i, j, k;



  
  // extraction
  
  vec_dense errors ( nbins );
  vector < vec_dense > vrinspec ( 1 );
  vector < vec_dense > vinspec ( 1 );
  vinspec[0].resize ( nbins );
  vrinspec[0].resize ( nbins );

  vinspec[0] = inspec;

  extract < mat_comprow, vec_dense > ( invcov, z, outspec, errors, vinspec, vrinspec );
  

  cerr << "final spectra" << endl;


  string outdata = datadir + "/test_medium_output.out";

  fstream out;
  out.open ( outdata.c_str(), ios::out );
  
  for ( size_t i = 0; i < nbins; ++i ) {
    out << i << " " << (vrinspec[0])[i] << " " << outspec[i] << " " << errors[i] << " " << ((vrinspec[0])[i] - outspec[i]) / errors[i] << endl;
  }
  
  out.close();


  cerr << "  (PASSED)" << endl;
     
  return;
}



