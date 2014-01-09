// @COPYRIGHT@

#ifndef HARP_MPI_IMAGE_HPP
#define HARP_MPI_IMAGE_HPP


namespace harp {

  class mpi_image {

    public :

      mpi_image ( boost::mpi::communicator const & comm, boost::property_tree::ptree const & props );
      
      ~mpi_image ( ) { }
      
      size_t n_rows ( ) const;
      
      size_t n_cols ( ) const;
      
      void values ( vector_double & data ) const;

      void inv_variance ( vector_double & invvar ) const;

      boost::property_tree::ptree metadata ( ) const;
      
      void values ( matrix_double & data ) const;

      void inv_variance ( matrix_double & invvar ) const;

      std::string format ( ) const;
      
    private :
    
      boost::mpi::communicator comm_;
      image_p local_;
      
  };
  
  typedef boost::shared_ptr < harp::mpi_image > mpi_image_p;
  typedef boost::weak_ptr < harp::mpi_image > mpi_image_wp;

}


#endif

