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

  void eigen_decompose ( matrix_dist & invcov, matrix_dist & D, matrix_dist & W );

  void eigen_compose ( eigen_op op, matrix_dist & D, matrix_dist & W, matrix_dist & out );

  void column_norm ( matrix_dist & mat, matrix_dist & S );

  void apply_norm ( matrix_dist & S, matrix_dist & mat );

  void norm ( matrix_dist & D, matrix_dist & W, matrix_dist & S );

  void sub_block ( matrix_dist & in, int firstrow, int firstcol, int nrow, int ncol, matrix_dist & out );

}

#endif

