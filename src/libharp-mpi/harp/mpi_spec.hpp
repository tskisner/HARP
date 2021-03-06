/*
  High Performance Astrophysical Reconstruction and Processing (HARP)
  (c) 2014-2015, The Regents of the University of California, 
  through Lawrence Berkeley National Laboratory.  See top
  level LICENSE file for details.
*/

#ifndef HARP_MPI_SPEC_HPP
#define HARP_MPI_SPEC_HPP


namespace harp {

  class mpi_spec {
    
    public :

      mpi_spec ( boost::mpi::communicator const & comm, std::string const & type, boost::property_tree::ptree const & props );
      
      ~mpi_spec ( ) { }

      boost::mpi::communicator comm ( ) { return comm_; }

      size_t n_spec ( ) const;

      size_t n_lambda ( ) const;

      void values ( vector_double & data ) const;

      void values ( mpi_matrix & data ) const;

      void lambda ( vector_double & lambda ) const;

      std::string type ( ) const;

      spec_p local ();
      
    private :

      boost::mpi::communicator comm_;
      spec_p local_;
      
  };

  typedef boost::shared_ptr < harp::mpi_spec > mpi_spec_p;
  typedef boost::weak_ptr < harp::mpi_spec > mpi_spec_wp;

}


#endif

