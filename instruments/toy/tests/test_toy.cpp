#include <iostream>

#include <harp_test.hpp>

using namespace std;
using namespace harp;


void harp::test_toy ( string const & datadir ) {
  
  cerr << "Testing toy image creation..." << endl;
  
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
  
  //for ( size_t i = 0; i < testpix->cols(); ++i ) {
  //  for ( size_t j = 0; j < testpix->rows(); ++j ) {
  //    cout << pix ( j, i ) << " ";
  //  }
  //  cout << endl;
  //}
  
  cerr << "  (PASSED)" << endl;
  
  
  cerr << "Testing toy spectrum creation..." << endl;
  
  params.clear();
  
  params[ "path" ] = datadir + "/test_spectra.fits";
  params[ "hdu" ] = "1";
  params[ "pos" ] = "3";
  
  spectrum_p testspec ( spectrum::create ( string("toy"), params ) );
  
  cerr << "  (PASSED)" << endl;
  
  
  cerr << "Testing toy spectrum read..." << endl;
  
  data_vec spec ( testspec->size() );
  
  testspec->read ( spec );

  //for ( size_t i = 0; i < spec.size(); ++i ) {
  //  cout << spec[i] << endl;
  //}
  
  cerr << "  (PASSED)" << endl;
  
  
  cerr << "Testing toy PSF creation..." << endl;
  
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
    //cerr << "( " << i << " )" << endl;
    //for ( size_t j = 0; j < lambda.size(); ++j ) {
    //  cerr << lambda[j] << " ";
    //}
    //cerr << endl << endl;
  }
  
  cerr << "  (PASSED)" << endl;
  
  
  cerr << "Testing toy PSF read of sparse projection matrix..." << endl;
  
  size_t nbins = testpsf->nspec() * testpsf->specsize(0);
  size_t npix = testpix->rows() * testpix->cols();
  
  sparse_mat projmat ( npix, nbins );
  
  sparse_mat_view projview ( projmat, mv_range ( 0, npix ), mv_range ( 0, nbins ) );
  
  testpsf->projection ( 0, testpsf->nspec() - 1, 0, testpsf->specsize(0) - 1, (size_t)0, testpix->cols() - 1, (size_t)0, testpix->rows() - 1, projview );
  
  data_vec inspec ( nbins );
  data_vec outvec ( testpix->rows() * testpix->cols() );
  
  for ( size_t b = 0; b < nbins; ++b ) {
    if ( (b+5) % 10 == 0 ) {
      inspec[b] = 5.0;
    } else {
      inspec[b] = 0.0;
    }
  }
  
  //cerr << "start axpy_prod" << endl;
  boost::numeric::ublas::axpy_prod ( projview, inspec, outvec, true );
  //cerr << "finished axpy_prod" << endl;
  
  dense_mat outmat ( testpix->rows(), testpix->cols() );
  dense_mat_view outview ( outmat, mv_range ( 0, testpix->rows() ), mv_range ( 0, testpix->cols() ) );
  
  size_t pixoff = 0;
  for ( size_t i = 0; i < testpix->rows(); ++i ) {
    for ( size_t j = 0; j < testpix->cols(); ++j ) {
      outmat( i, j ) = outvec[pixoff];
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
  
  image_p outimage ( image::create ( string("toy"), params ) );
  
  outimage->write ( "!" + datadir + "/test_psf_prod.fits", 0, 0, outview );
  
  
  
  cerr << "  (PASSED)" << endl;
  
  
  
  return;
}