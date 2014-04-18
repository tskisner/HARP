// @COPYRIGHT@

#ifndef HARP_MPI_PSF_HPP
#define HARP_MPI_PSF_HPP


namespace harp {

  class mpi_psf {
    
    public :

      mpi_psf ( boost::mpi::communicator const & comm, std::string const & type, boost::property_tree::ptree const & props );

      mpi_psf ( ) { }

      ~mpi_psf ( ) { }

      mpi_psf * redistribute ( boost::mpi::communicator const & comm );

      boost::mpi::communicator comm ( ) { return comm_; }

      size_t n_spec ( ) const;
      
      size_t n_lambda ( ) const;

      size_t img_rows ( ) const;

      size_t img_cols ( ) const;

      vector_double lambda ( ) const;

      void extent ( size_t spec, size_t lambda, size_t & x_offset, size_t & y_offset, size_t & n_x, size_t & n_y ) const;

      // The base class provides a default implementation of these 2 functions, so that derived classes
      // only need to implement the response() method.
      
      void project_transpose ( std::map < size_t, std::set < size_t > > const & speclambda, mpi_matrix_sparse & AT ) const;

      // These are convenience functions if you want the whole projection matrix

      void project_transpose ( mpi_matrix_sparse & AT ) const;
      
      std::string type ( ) const;
      
    private :
    
      boost::mpi::communicator comm_;
      psf_p local_;      
      
  };
  
  typedef boost::shared_ptr < harp::mpi_psf > mpi_psf_p;
  typedef boost::weak_ptr < harp::mpi_psf > mpi_psf_wp;

}


#endif

