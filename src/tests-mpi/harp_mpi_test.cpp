
#include <iostream>

#include <harp_mpi_test.hpp>

using namespace std;
using namespace harp;


int main ( int argc, char *argv[] ) {

  El::Initialize ( argc, argv );

  boost::mpi::environment env;
  boost::mpi::communicator comm;

  int np = comm.size();
  int myp = comm.rank();


  string datadir = "../tests/testdata";

  if ( myp == 0 ) {
    cout << endl;
  }
  
  mpi_test_linalg ( datadir );

  mpi_test_specslice ( datadir );

  mpi_test_spec ( datadir );

  mpi_test_extract ( datadir );

  
  if ( myp == 0 ) {
    cout << endl;
  }

  El::Finalize();
  
  return 0;
}
