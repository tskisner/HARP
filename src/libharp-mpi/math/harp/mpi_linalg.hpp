// @COPYRIGHT@

#ifndef HARP_MPI_LINALG_HPP
#define HARP_MPI_LINALG_HPP

#include <elemental.hpp>


namespace harp {

  // convenience typedefs

  typedef elem::Matrix < double > elem_matrix_local;

  typedef elem::DistMatrix < double, elem::MC, elem::MR > elem_matrix;

  typedef elem_matrix mpi_matrix;


  // distributed sparse matrix

  class mpi_matrix_sparse_block {

    friend class boost::serialization::access;

    public :

      mpi_matrix_sparse_block ( );

      ~mpi_matrix_sparse_block () { }

      size_t firstrow;
      size_t rows;
      size_t vals;
      std::vector < size_t > row;
      std::vector < size_t > col;
      std::vector < size_t > row_offset;
      std::vector < size_t > row_nnz;
      std::vector < double > data;

    private :

      template < class Archive >
      void serialize ( Archive & ar, const unsigned int version ) {
        ar & BOOST_SERIALIZATION_NVP(firstrow);
        ar & BOOST_SERIALIZATION_NVP(rows);
        ar & BOOST_SERIALIZATION_NVP(vals);
        ar & BOOST_SERIALIZATION_NVP(row);
        ar & BOOST_SERIALIZATION_NVP(col);
        ar & BOOST_SERIALIZATION_NVP(row_offset);
        ar & BOOST_SERIALIZATION_NVP(row_nnz);
        ar & BOOST_SERIALIZATION_NVP(data);
        return;
      }

  };


  class mpi_matrix_sparse {

    public :

      mpi_matrix_sparse ( boost::mpi::communicator const & comm, size_t rows, size_t cols );

      ~mpi_matrix_sparse () { }

      boost::mpi::communicator const & comm ( ) { return comm_; }

      mpi_matrix_sparse_block & block ( ) { return block_; }

      size_t rows () { return rows_; }
      size_t cols () { return cols_; }

    private :

      boost::mpi::communicator comm_;
      mpi_matrix_sparse_block block_;
      size_t rows_;
      size_t cols_;

  };


  void mpi_matrix_zero ( mpi_matrix & mat );

  void local_matrix_zero ( elem_matrix_local & mat );

  void mpi_eigen_decompose ( mpi_matrix const & invcov, mpi_matrix & D, mpi_matrix & W );

  void mpi_eigen_compose ( eigen_op op, mpi_matrix const & D, mpi_matrix const & W, mpi_matrix & out );

  void mpi_column_norm ( mpi_matrix const & mat, mpi_matrix & S );

  void mpi_apply_norm ( mpi_matrix const & S, mpi_matrix & mat );

  void mpi_apply_inverse_norm ( mpi_matrix const & S, mpi_matrix & mat );

  void mpi_norm ( mpi_matrix const & D, mpi_matrix const & W, mpi_matrix & S );

  void mpi_gang_distribute ( mpi_matrix const & mat, mpi_matrix & gmat );

  void mpi_gang_accum ( mpi_matrix const & gmat, mpi_matrix & mat );


}


#endif
