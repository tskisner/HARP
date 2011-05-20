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

#include <moat.hpp>

namespace harp {
  
  // typedef for sparse row major matrix
  typedef boost::numeric::ublas::mapped_matrix < double, boost::numeric::ublas::row_major > sparse_mat;
  typedef boost::numeric::ublas::matrix_range < sparse_mat > sparse_mat_view;
  
  // typedef for dense row major matrix
  typedef boost::numeric::ublas::matrix < double, boost::numeric::ublas::row_major > dense_mat;
  typedef boost::numeric::ublas::matrix_range < dense_mat > dense_mat_view;
  
  // typedef for dense vector
  typedef boost::numeric::ublas::vector < double > data_vec;
  typedef boost::numeric::ublas::vector_range < data_vec > data_vec_view;
  
  // typedef for int vector
  typedef boost::numeric::ublas::vector < int > int_vec;
  typedef boost::numeric::ublas::vector_range < int_vec > int_vec_view;
  
  // typedef for range
  typedef boost::numeric::ublas::range mv_range;
  
}


#include <harp/fits.hpp>

#include <harp/image.hpp>
#include <harp/spectrum.hpp>
#include <harp/psf.hpp>

#endif



