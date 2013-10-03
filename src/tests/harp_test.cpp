
#include <iostream>

#include <harp_test.hpp>

using namespace std;
using namespace harp;


int main ( int argc, char *argv[] ) {

  string datadir = "testdata";
  
  fits::test ( datadir );

  test_spec_specter ( datadir );

  test_spec_sim ( datadir );

  //test_image_sim ( datadir );

  //test_image_fits ( datadir );

  //test_psf_gauss ( datadir );

  //test_psf_gauss_sim ( datadir );
  
  return 0;
}
