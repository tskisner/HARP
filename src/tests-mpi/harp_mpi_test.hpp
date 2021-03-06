/*
  High Performance Astrophysical Reconstruction and Processing (HARP)
  (c) 2014-2015, The Regents of the University of California, 
  through Lawrence Berkeley National Laboratory.  See top
  level LICENSE file for details.
*/

#ifndef HARP_MPI_TEST_HPP
#define HARP_MPI_TEST_HPP

#include <harp_mpi_internal.hpp>


namespace harp {

  void mpi_test_linalg ( std::string const & datadir );

  void mpi_test_specslice ( std::string const & datadir );

  void mpi_test_spec ( std::string const & datadir );

  void mpi_test_extract ( std::string const & datadir );

  
}

  
#endif