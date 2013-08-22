// @COPYRIGHT@

#ifndef HARP_MATRIX_HPP
#define HARP_MATRIX_HPP


namespace harp {

  class sparse_block {

    public :
      sparse_block ( matrix_sparse const & orig );
      sparse_block ( char * packed, size_t nbytes );
      ~sparse_block () { }

      char * pack ( size_t & nbytes );

      int global_rows;
      int local_firstrow;
      int local_rows;
      int local_vals;
      std::vector < int > local_row;
      std::vector < int > local_col;
      std::vector < int > local_row_offset;
      std::vector < int > local_row_nnz;
      std::vector < double > data;

  };

  void dist_matrix_zero ( matrix_dist & mat );

  void local_matrix_zero ( matrix_local & mat );

  void eigen_decompose ( matrix_dist const & invcov, matrix_dist & D, matrix_dist & W );

  void eigen_compose ( eigen_op op, matrix_dist const & D, matrix_dist const & W, matrix_dist & out );

  void column_norm ( matrix_dist const & mat, matrix_dist & S );

  void apply_norm ( matrix_dist const & S, matrix_dist & mat );

  void apply_inverse_norm ( matrix_dist const & S, matrix_dist & mat );

  void norm ( matrix_dist const & D, matrix_dist const & W, matrix_dist & S );

  void gang_distribute ( matrix_dist const & mat, matrix_dist & gmat );

  void gang_accum ( matrix_dist const & gmat, matrix_dist & mat );

}

#endif
