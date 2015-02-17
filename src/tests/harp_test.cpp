
#include <iostream>

#include <harp_test.hpp>

using namespace std;
using namespace harp;


int main ( int argc, char *argv[] ) {

  string datadir = "testdata";
  
  fits::test ( datadir );

  test_serialize ( datadir );

  test_linalg ( datadir );

  test_spec_simspecter ( datadir );

  test_psf_gauss_sim ( datadir );

  test_psf_gauss ( datadir );

  test_image_simfits ( datadir );

  test_desi ( datadir );

  test_specslice ( datadir );

  test_extract ( datadir );

  return 0;
}
