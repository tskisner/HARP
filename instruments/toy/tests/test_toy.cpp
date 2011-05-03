#include <iostream>

#include <harp_test.hpp>

using namespace std;
using namespace harp;


void harp::test_toy ( string const & datadir ) {
  
  cerr << "Testing toy image creation..." << endl; 
  
  std::map < std::string, std::string > params;
  
  params[ "path" ] = datadir + "/test_image.fits";
  
  image_p testimg ( image::create ( string("toy"), params ) );
  
  cerr << "  (PASSED)" << endl;
  
  
  
  
  return;
}