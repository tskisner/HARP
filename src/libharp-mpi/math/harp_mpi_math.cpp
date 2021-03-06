/*
  High Performance Astrophysical Reconstruction and Processing (HARP)
  (c) 2014-2015, The Regents of the University of California, 
  through Lawrence Berkeley National Laboratory.  See top
  level LICENSE file for details.
*/

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

  if ( myn == 0 ) {
    offset = 0;
  }
  
  return;
}
