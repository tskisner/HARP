// @COPYRIGHT@

#ifndef HARP_MPI_SPEC_SLICE_HPP
#define HARP_MPI_SPEC_SLICE_HPP


namespace harp {

  class mpi_spec_slice {
    
    public :

      mpi_spec_slice ( boost::mpi::communicator const & comm, size_t nspec, size_t nlambda, size_t chunk_nspec, size_t chunk_nlambda, size_t overlap_spec, size_t overlap_lambda );
      
      ~mpi_spec_slice ( ) { }

      spec_slice_region const & full_region ( ) const { return local_->full_region(); }

      std::vector < spec_slice_region > regions ( ) const;
      
    private :

      boost::mpi::communicator comm_;
      spec_slice_p local_;

  };
  
  typedef boost::shared_ptr < harp::mpi_spec_slice > mpi_spec_slice_p;

}


#endif

