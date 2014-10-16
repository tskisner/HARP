/*
  High Performance Astrophysical Reconstruction and Processing (HARP)
  (c) 2014-2015, The Regents of the University of California, 
  through Lawrence Berkeley National Laboratory.  See top
  level LICENSE file for details.
*/

#ifndef HARP_MPI_PSF_HPP
#define HARP_MPI_PSF_HPP


namespace harp {

  class mpi_psf {
    
    public :

      mpi_psf ( boost::mpi::communicator const & comm, std::string const & type, boost::property_tree::ptree const & props );

      mpi_psf ( ) { }

      ~mpi_psf ( ) { }

      mpi_psf * redistribute ( boost::mpi::communicator const & newcomm );

      boost::mpi::communicator comm ( ) { return comm_; }

      size_t n_spec ( ) const;
      
      size_t n_lambda ( ) const;

      size_t img_rows ( ) const;

      size_t img_cols ( ) const;

      vector_double lambda ( ) const;

      void extent ( size_t spec, size_t lambda, size_t & x_offset, size_t & y_offset, size_t & n_x, size_t & n_y ) const;

      void extent_multi ( std::map < size_t, std::set < size_t > > const & speclambda, std::vector < size_t > & x_offset, std::vector < size_t > & y_offset, std::vector < size_t > & n_x, std::vector < size_t > & n_y ) const;

      void project_transpose ( std::map < size_t, std::set < size_t > > const & speclambda, mpi_matrix_sparse & AT ) const;

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

