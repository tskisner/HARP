// @COPYRIGHT@

#include <harp_mpi_internal.hpp>


using namespace std;
using namespace harp;


harp::mpi_spec::mpi_spec ( boost::mpi::communicator const & comm, std::string const & type, boost::property_tree::ptree const & props ) {

  comm_ = comm;

  int rank = comm.rank();

  if ( rank == 0 ) {
    // ONLY rank zero should use the plugin registry!
    plugin_registry & reg = plugin_registry::get();

    // instantiate
    local_.reset ( reg.create_spec ( type, props ) );
  }
    
  // broadcast to all processes

  boost::mpi::broadcast ( comm_, local_, 0 );

}


size_t harp::mpi_spec::n_spec ( ) const {
  return local_->n_spec();
}


size_t harp::mpi_spec::n_lambda ( ) const {
  return local_->n_lambda();
}


void harp::mpi_spec::values ( vector_double & data ) const {
  local_->values ( data );
  return;
}


void harp::mpi_spec::lambda ( vector_double & lambda_vals ) const {
  local_->lambda ( lambda_vals );
  return;
}


void harp::mpi_spec::targets ( std::vector < target > & target_list ) const {
  local_->targets ( target_list );
  return;
}


std::string harp::mpi_spec::type ( ) const {
  return local_->type();
}




