
#include <iostream>

#include <harp_test.hpp>

using namespace std;
using namespace harp;


int main ( int argc, char *argv[] ) {

  cliq::Initialize( argc, argv );

  int np;
  int myp;

  MPI_Comm_size ( MPI_COMM_WORLD, &np );
  MPI_Comm_rank ( MPI_COMM_WORLD, &myp );

  string datadir = "testdata";

  if ( myp == 0 ) {
    cerr << endl;
  }
  
  fits::test ( datadir );

  test_elemental ( datadir );

  test_invcov ( datadir );

  //test_tinyKLT ( datadir );

  test_spec_specter ( datadir );

  test_spec_sim ( datadir );

  test_psf_gauss ( datadir );

  test_sim_extract ( datadir );

  test_small_extract ( datadir );
  
  if ( myp == 0 ) {
    cerr << endl;
  }

  cliq::Finalize();
  
  return 0;
}
