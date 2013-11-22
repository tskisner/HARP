
#include <iostream>

#include <harp_test.hpp>

using namespace std;
using namespace harp;


// serialization class exports (should only be placed in front of a main() !!! )

#include <harp_serialization.hpp>


int main ( int argc, char *argv[] ) {

  string datadir = "testdata";
  
  fits::test ( datadir );

  test_serialize ( datadir );

  test_lapack ( datadir );

  test_spec_simspecter ( datadir );

  test_psf_gauss_sim ( datadir );

  test_psf_gauss ( datadir );

  test_image_simfits ( datadir );

  return 0;
}
