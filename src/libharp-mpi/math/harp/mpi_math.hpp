// @COPYRIGHT@

#ifndef HARP_MPI_MATH_HPP
#define HARP_MPI_MATH_HPP

#include <mpi.h>

#include <boost/mpi.hpp>

#include <harp/math.hpp>


namespace harp {


  #define HARP_MPI_ABORT(proc, msg) \
  std::cerr << "File " << __FILE__ << ", Line " << __LINE__ << ", Proc " << proc << " : " << msg << std::endl; \
  MPI_Abort ( MPI_COMM_WORLD, 1 ) \


  void mpi_dist_1D ( boost::mpi::communicator const & comm, size_t n, size_t & myn, size_t & offset );


}

#include <harp/mpi_linalg.hpp>
#include <harp/mpi_comm.hpp>

#endif
