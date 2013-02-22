#include <iostream>
#include <fstream>
#include <sstream>

#include <harp_test.hpp>

#include <boost/random.hpp>


extern "C" {
#include <unistd.h>
}

using namespace std;
using namespace harp;


void harp::test_boss ( string const & datadir ) {

  int np;
  int myp;

  MPI_Comm_size ( MPI_COMM_WORLD, &np );
  MPI_Comm_rank ( MPI_COMM_WORLD, &myp );
  
  if ( myp == 0 ) {
    cerr << "Testing boss_specter operations..." << endl;
  }

  if ( myp == 0 ) {
    cerr << "  (PASSED)" << endl;
  }
     
  return;
}



