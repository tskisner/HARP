
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

  // run built-in tests

  if ( myp == 0 ) {
    cerr << endl;
  }
  
  fits::test ( datadir );

  test_elemental ( datadir );

  test_invcov ( datadir );

  //test_tinyKLT ( datadir );
  
  if ( myp == 0 ) {
    cerr << endl;
  }
  

  // run format tests
  
# include "harp_testcommands.cpp"

  //cerr << endl;

  cliq::Finalize();
  //MPI_Finalize ();
  
  return 0;
}
