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
  
  // typedef for sparse mapped row major matrix
  typedef boost::numeric::ublas::mapped_matrix < double, boost::numeric::ublas::row_major > sparse_rowmat;
  typedef boost::numeric::ublas::matrix_range < sparse_rowmat > sparse_rowmat_view;
  
  // typedef for sparse mapped column major matrix
  typedef boost::numeric::ublas::mapped_matrix < double, boost::numeric::ublas::column_major > sparse_colmat;
  typedef boost::numeric::ublas::matrix_range < sparse_colmat > sparse_colmat_view;
  
  // typedef for compressed row major matrix
  typedef boost::numeric::ublas::compressed_matrix < double, boost::numeric::ublas::row_major > comp_rowmat;
  typedef boost::numeric::ublas::matrix_range < comp_rowmat > comp_rowmat_view;
  
  // typedef for compressed column major matrix
  typedef boost::numeric::ublas::compressed_matrix < double, boost::numeric::ublas::column_major > comp_colmat;
  typedef boost::numeric::ublas::matrix_range < comp_colmat > comp_colmat_view;
  
  // typedef for dense row major matrix
  typedef boost::numeric::ublas::matrix < double, boost::numeric::ublas::row_major > dense_rowmat;
  typedef boost::numeric::ublas::matrix_range < dense_rowmat > dense_rowmat_view;
  
  // typedef for dense column major matrix
  typedef boost::numeric::ublas::matrix < double, boost::numeric::ublas::column_major > dense_colmat;
  typedef boost::numeric::ublas::matrix_range < dense_colmat > dense_colmat_view;
  
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
#include <harp/spec.hpp>
#include <harp/psf.hpp>

#endif



