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
  
  
  
  return;
}