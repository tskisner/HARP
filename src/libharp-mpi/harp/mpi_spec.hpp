// @COPYRIGHT@

#ifndef HARP_MPI_SPEC_HPP
#define HARP_MPI_SPEC_HPP


namespace harp {

  class mpi_spec {
    
    public :

      mpi_spec ( boost::mpi::communicator const & comm, boost::property_tree::ptree const & props );
      
      ~mpi_spec ( ) { }

      size_t n_spec ( ) const;

      size_t n_lambda ( ) const;

      void values ( vector_double & data ) const;

      void lambda ( vector_double & lambda ) const;

      void targets ( std::vector < target > & target_list ) const;

      std::string format ( ) const;
      
    private :

      boost::mpi::communicator comm_;
      spec_p local_;
      
  };

  typedef boost::shared_ptr < harp::mpi_spec > mpi_spec_p;
  typedef boost::weak_ptr < harp::mpi_spec > mpi_spec_wp;

}


#endif
