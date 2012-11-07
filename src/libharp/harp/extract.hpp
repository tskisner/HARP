// @COPYRIGHT@

#ifndef HARP_EXTRACT_HPP
#define HARP_EXTRACT_HPP


namespace harp {

  // helpers for communicating sparse matrix blocks

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

  void spec_project ( matrix_sparse const & m, matrix_dist const & in, matrix_local & out );

  void noise_weighted_spec ( matrix_sparse const & psf, matrix_local const & invnoise, matrix_local const & img, matrix_dist & z );

  void inverse_covariance ( matrix_sparse const & psf, matrix_local const & invnoise, matrix_dist & invcov );

  void eigenpairs ( matrix_dist & invcov, matrix_dist & D, matrix_dist & W );

  void norm ( matrix_dist & D, matrix_dist & W, matrix_dist & S );

  void resolution ( matrix_dist & D, matrix_dist & W, matrix_dist & S, matrix_dist & R );

  void extract ( matrix_dist & D, matrix_dist & W, matrix_dist & S, matrix_dist & z, matrix_dist & f );


}

#endif

