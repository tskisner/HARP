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
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/banded.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>

#include <boost/property_tree/ptree.hpp>


#include <moat.hpp>



namespace harp {
  
  typedef boost::numeric::ublas::compressed_matrix < double, boost::numeric::ublas::row_major > mat_comprow;
  typedef boost::numeric::ublas::compressed_matrix < double, boost::numeric::ublas::column_major > mat_compcol;

  typedef boost::numeric::ublas::mapped_matrix < double, boost::numeric::ublas::row_major > mat_dynrow;
  typedef boost::numeric::ublas::mapped_matrix < double, boost::numeric::ublas::column_major > mat_dyncol;

  typedef boost::numeric::ublas::diagonal_matrix < double > mat_diag;
  
  typedef boost::numeric::ublas::matrix < double, boost::numeric::ublas::row_major > mat_denserow;
  typedef boost::numeric::ublas::matrix < double, boost::numeric::ublas::column_major > mat_densecol;

  typedef boost::numeric::ublas::mapped_vector < double > vec_sparse;
  typedef boost::numeric::ublas::vector < double > vec_dense;
  typedef boost::numeric::ublas::mapped_vector < uint8_t > vec_flag;
  typedef boost::numeric::ublas::vector < int > vec_denseint;
  
}


#include <harp/fits.hpp>

#include <harp/image.hpp>
#include <harp/spec.hpp>
#include <harp/psf.hpp>

#include <harp/metadata.hpp>

#endif



