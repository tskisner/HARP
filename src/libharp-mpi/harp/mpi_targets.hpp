// @COPYRIGHT@

#ifndef HARP_MPI_TARGETS_HPP
#define HARP_MPI_TARGETS_HPP


namespace harp {

  class mpi_targets {
    
    public :

      mpi_targets ( boost::mpi::communicator const & comm, std::string const & type, boost::property_tree::ptree const & props );
      
      ~mpi_targets ( ) { }

      boost::mpi::communicator comm ( ) { return comm_; }

      size_t n_objects ( ) const;

      std::vector < object_p > objects ( ) const;

      std::string type ( ) const;
      
    private :

      boost::mpi::communicator comm_;
      targets_p local_;
      
  };

  typedef boost::shared_ptr < harp::mpi_targets > mpi_targets_p;
  typedef boost::weak_ptr < harp::mpi_targets > mpi_targets_wp;

}


#endif

