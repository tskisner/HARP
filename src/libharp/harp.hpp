// @COPYRIGHT@

#ifndef HARP_HPP
#define HARP_HPP

#include <string>
#include <memory>
#include <vector>
#include <map>
#include <sstream>

#include <boost/smart_ptr/shared_ptr.hpp>
#include <boost/smart_ptr/weak_ptr.hpp>
#include <boost/enable_shared_from_this.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>

extern "C" {
  #include <fitsio.h>
}

#include <moat.hpp>

namespace harp {
  
  // typedef for sparse row major matrix
  typedef boost::numeric::ublas::compressed_matrix < double, boost::numeric::ublas::row_major > sparse_mat;
  typedef boost::numeric::ublas::matrix_range < sparse_mat > sparse_mat_view;
  
  // typedef for dense row major matrix
  typedef boost::numeric::ublas::matrix < double, boost::numeric::ublas::row_major > dense_mat;
  typedef boost::numeric::ublas::matrix_range < dense_mat > dense_mat_view;
  
}

#include <harp/image.hpp>
#include <harp/spectra.hpp>
#include <harp/psf.hpp>

#endif



