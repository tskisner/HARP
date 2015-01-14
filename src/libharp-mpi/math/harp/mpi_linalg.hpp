/*
  High Performance Astrophysical Reconstruction and Processing (HARP)
  (c) 2014-2015, The Regents of the University of California, 
  through Lawrence Berkeley National Laboratory.  See top
  level LICENSE file for details.
*/

#ifndef HARP_MPI_LINALG_HPP
#define HARP_MPI_LINALG_HPP

#include <El.hpp>


namespace harp {

  // convenience typedefs

  typedef El::Matrix < double > elem_matrix_local;

  typedef El::DistMatrix < double, El::MC, El::MR > elem_matrix;

  typedef elem_matrix mpi_matrix;


  // conversion routines between ublas and elemental column vectors

  template < class V >
  void ublas_to_elem ( boost::numeric::ublas::vector_expression < V > const & in, elem_matrix_local & out ) {
    typedef V vector_type;

    out.Resize ( in().size(), 1 );
    for ( size_t i = 0; i < in().size(); ++i ) {
      out.Set ( i, 0, in()[i] );
    }

    return;
  }

  template < class V >
  void ublas_to_elem ( boost::numeric::ublas::vector_expression < V > const & in, elem_matrix & out ) {
    typedef V vector_type;

    out.Resize ( in().size(), 1 );

    // populate local elements

    size_t hlocal = out.LocalHeight();
    size_t wlocal = out.LocalWidth();

    size_t rowoff = out.ColShift();
    size_t rowstride = out.ColStride();
    size_t row;

    for ( size_t i = 0; i < wlocal; ++i ) {
      for ( size_t j = 0; j < hlocal; ++j ) {
        row = rowoff + j * rowstride;
        out.SetLocal ( j, i, in()[ row ] );
      }
    }

    return;
  }

  template < class V >
  void elem_to_ublas ( elem_matrix_local const & in, boost::numeric::ublas::vector_expression < V > & out ) {
    typedef V vector_type;

    if ( in.Width() != 1 ) {
      HARP_MPI_ABORT( 0, "elem_to_ublas only works with single-column matrices" );
    }

    out().resize ( in.Height() );

    for ( size_t i = 0; i < in.Height(); ++i ) {
      out()[i] = in.Get ( i, 0 );
    }

    return;
  }

  template < class V >
  void elem_to_ublas ( elem_matrix const & in, boost::numeric::ublas::vector_expression < V > & out ) {
    typedef V vector_type;

    if ( in.Width() != 1 ) {
      HARP_MPI_ABORT( 0, "elem_to_ublas only works with single-column matrices" );
    }

    // get a local copy

    elem_matrix_local local ( in.Height(), 1 );

    El::AxpyInterface < double > globloc;
    globloc.Attach( El::GLOBAL_TO_LOCAL, in );
    globloc.Axpy ( 1.0, local, 0, 0 );
    globloc.Detach();

    // copy to output

    elem_to_ublas ( local, out() );

    return;
  }


  // distributed sparse matrix

  class mpi_matrix_sparse_block {

    friend class boost::serialization::access;

    public :

      mpi_matrix_sparse_block ( ) { }

      ~mpi_matrix_sparse_block ( ) { }

      // default copy and assignment operators are sufficient...

      size_t firstrow;
      size_t rows;
      size_t vals;
      std::vector < size_t > row;
      std::vector < size_t > col;
      std::vector < size_t > row_offset;
      std::vector < size_t > row_nnz;
      vector_double data;

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

      void clear ();

      boost::mpi::communicator const & comm ( ) const { return comm_; }

      mpi_matrix_sparse_block & block ( ) { return block_; }

      mpi_matrix_sparse_block const & block ( ) const { return block_; }

      size_t rows () const { return rows_; }
      size_t cols () const { return cols_; }

    private :

      boost::mpi::communicator comm_;
      mpi_matrix_sparse_block block_;
      size_t rows_;
      size_t cols_;

  };


  void mpi_matrix_zero ( mpi_matrix & mat );

  void local_matrix_zero ( elem_matrix_local & mat );

  void mpi_eigen_decompose ( mpi_matrix const & invcov, mpi_matrix & D, mpi_matrix & W, bool regularize );

  void mpi_eigen_compose ( eigen_op op, mpi_matrix const & D, mpi_matrix const & W, mpi_matrix & out );

  void mpi_column_norm ( mpi_matrix const & mat, mpi_matrix & S );

  void mpi_apply_norm ( mpi_matrix const & S, mpi_matrix & mat );

  void mpi_apply_inverse_norm ( mpi_matrix const & S, mpi_matrix & mat );

  void mpi_norm ( mpi_matrix const & D, mpi_matrix const & W, mpi_matrix & S );

  void mpi_sparse_mv_trans ( mpi_matrix_sparse const & AT, mpi_matrix const & in, elem_matrix_local & out );

  void mpi_gang_distribute ( mpi_matrix const & mat, mpi_matrix & gmat );

  void mpi_gang_accum ( mpi_matrix const & gmat, mpi_matrix & mat );


}


#endif
