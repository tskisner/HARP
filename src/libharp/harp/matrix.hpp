// @COPYRIGHT@

#ifndef HARP_MATRIX_HPP
#define HARP_MATRIX_HPP


namespace harp {

  typedef boost::numeric::ublas::matrix < double, boost::numeric::ublas::column_major > matrix_double;

  typedef boost::numeric::ublas::compressed_matrix < double, boost::numeric::ublas::row_major > matrix_double_sparse;

  typedef boost::numeric::ublas::vector < double > vector_double;


  void check_column_major ( boost::numeric::ublas::column_major_tag );


  void check_column_major ( boost::numeric::ublas::row_major_tag );


  template < class M >
  void check_column_major ( boost::numeric::ublas::matrix_expression < M > const & matrix ) {
    typedef typename M::orientation_category orientation_category;
    check_column_major ( orientation_category() );
    return;
  }


  /*
  void eigen_decompose ( matrix_serial const & invcov, matrix_dist & D, matrix_dist & W );

  void eigen_compose ( eigen_op op, matrix_dist const & D, matrix_dist const & W, matrix_dist & out );

  void column_norm ( matrix_dist const & mat, matrix_dist & S );

  void apply_norm ( matrix_dist const & S, matrix_dist & mat );

  void apply_inverse_norm ( matrix_dist const & S, matrix_dist & mat );

  void norm ( matrix_dist const & D, matrix_dist const & W, matrix_dist & S );

  void gang_distribute ( matrix_dist const & mat, matrix_dist & gmat );

  void gang_accum ( matrix_dist const & gmat, matrix_dist & mat );

  */

}

#endif

