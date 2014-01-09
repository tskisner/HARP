// @COPYRIGHT@

#include <harp_mpi_internal.hpp>


using namespace std;
using namespace harp;


harp::mpi_spec::mpi_spec ( boost::mpi::communicator const & comm, boost::property_tree::ptree const & props ) {

  comm_ = comm;

  int rank = comm.rank();

  if ( rank == 0 ) {
    // instantiate
    local_.reset ( spec::create ( props ) );
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


std::string harp::mpi_spec::format ( ) const {
  return local_->format();
}




