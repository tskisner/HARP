// @COPYRIGHT@

#ifndef HARP_LINALG_HPP
#define HARP_LINALG_HPP


#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/storage.hpp>

#include <boost/numeric/bindings/ublas.hpp>
#include <boost/numeric/bindings/blas.hpp>
#include <boost/numeric/bindings/lapack.hpp>
#include <boost/numeric/bindings/views.hpp>


namespace harp {

  typedef enum {
    EIG_NONE,
    EIG_SQRT,
    EIG_INVSQRT
  } eigen_op;

  typedef boost::numeric::ublas::matrix < double, boost::numeric::ublas::column_major > matrix_double;

  typedef boost::numeric::ublas::compressed_matrix < double, boost::numeric::ublas::row_major > matrix_double_sparse;

  typedef boost::numeric::ublas::vector < double > vector_double;

  typedef boost::numeric::ublas::matrix < float, boost::numeric::ublas::column_major > matrix_float;

  typedef boost::numeric::ublas::compressed_matrix < float, boost::numeric::ublas::row_major > matrix_float_sparse;

  typedef boost::numeric::ublas::vector < float > vector_float;


  void check_column_major ( boost::numeric::ublas::column_major_tag );


  void check_column_major ( boost::numeric::ublas::row_major_tag );


  template < class M >
  void check_column_major ( boost::numeric::ublas::matrix_expression < M > const & matrix ) {
    typedef typename M::orientation_category orientation_category;
    check_column_major ( orientation_category() );
    return;
  }


  void eigen_decompose ( matrix_double const & invcov, vector_double & D, matrix_double & W );

  void eigen_compose ( eigen_op op, vector_double const & D, matrix_double const & W, matrix_double & out );

  void column_norm ( matrix_double const & mat, vector_double & S );

  void apply_norm ( vector_double const & S, matrix_double & mat );

  void apply_inverse_norm ( vector_double const & S, matrix_double & mat );

  void apply_norm ( vector_double const & S, vector_double & vec );

  void apply_inverse_norm ( vector_double const & S, vector_double & vec );

  void norm ( vector_double const & D, matrix_double const & W, vector_double & S );


}


#endif
