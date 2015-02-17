/*
  High Performance Astrophysical Reconstruction and Processing (HARP)
  (c) 2014-2015, The Regents of the University of California, 
  through Lawrence Berkeley National Laboratory.  See top
  level LICENSE file for details.
*/

#include <harp_mpi_internal.hpp>


using namespace std;
using namespace harp;


harp::mpi_targets::mpi_targets ( boost::mpi::communicator const & comm, std::string const & type, boost::property_tree::ptree const & props ) {

  comm_ = comm;

  int rank = comm.rank();

  plugin_registry & reg = plugin_registry::get();

  if ( rank == 0 ) {
    // instantiate
    local_.reset ( reg.create_targets ( type, props ) );
  }

  // broadcast to all processes

  mpi_comm_bcast ( comm_, local_, 0 );
  
}


size_t harp::mpi_targets::n_objects ( ) const {
  return local_->n_objects();
}


vector < object_p > harp::mpi_targets::objects ( ) const {
  return local_->objects();
}


std::string harp::mpi_targets::type ( ) const {
  return local_->type();
}


targets_p harp::mpi_targets::local ( ) {
  return local_;
}



