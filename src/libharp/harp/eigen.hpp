// @COPYRIGHT@

#ifndef HARP_EIGEN_HPP
#define HARP_EIGEN_HPP

#include <iostream>
#include <fstream>
#include <climits>

#include <moat.hpp>


namespace harp {

  int lapack_ev ( int dim, double * mat, double * eigen );

  template < typename M, typename V >
  void extract ( M & inout, V & z, V & spectra, V & errors, std::vector < V > & intestspec, std::vector < V > & outtestspec ) {



    return;
  }


}

#endif
