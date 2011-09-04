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

#include <moat.hpp>

#include <fftw3.h>

#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>

namespace harp {
  
  typedef Eigen::SparseMatrix < double, Eigen::RowMajor > mat_comprow;
  typedef Eigen::DynamicSparseMatrix < double, Eigen::RowMajor > mat_dynrow;
  typedef Eigen::DiagonalMatrix < double, Eigen::Dynamic > mat_diag;
  typedef Eigen::Matrix < double, Eigen::Dynamic, Eigen::Dynamic > mat_denserow;
  
  typedef Eigen::SparseVector < double > vec_sparse;
  typedef Eigen::Matrix < double, Eigen::Dynamic, 1 > vec_dense;
  typedef Eigen::Matrix < uint8_t, Eigen::Dynamic, 1 > vec_flag;
  typedef Eigen::Matrix < int, Eigen::Dynamic, 1 > vec_int;
  
}


#include <harp/fits.hpp>

#include <harp/image.hpp>
#include <harp/spec.hpp>
#include <harp/psf.hpp>

#endif



