// @COPYRIGHT@

#include <harp_mpi_math_internal.hpp>

extern "C" {
  #include <sys/time.h>
}


using namespace std;
using namespace harp;


void harp::mpi_dist_1D ( boost::mpi::communicator const & comm, size_t n, size_t & myn, size_t & offset ) {

  int rank = comm.rank();
  int nproc = comm.size();

  myn = (size_t) ( n / (size_t)(nproc) );
  
  size_t leftover = n % (size_t)(nproc);
  
  if ( (size_t)rank < leftover ) {
    ++myn;
    offset = myn * (size_t)rank;
  } else {
    offset = ( ( myn + 1 ) * leftover ) + ( myn * ( (size_t)rank - leftover ) );
  }
  
  return;
}
