// @COPYRIGHT@

#include <harp_mpi_internal.hpp>


using namespace std;
using namespace harp;


harp::mpi_image::mpi_image ( boost::mpi::communicator const & comm, std::string const & type, boost::property_tree::ptree const & props ) {

  comm_ = comm;

  int rank = comm.rank();

  if ( rank == 0 ) {
    // ONLY rank zero should use the plugin registry!
    plugin_registry & reg = plugin_registry::get();

    // instantiate
    local_.reset ( reg.create_image ( type, props ) );
  }
    
  // broadcast to all processes

  boost::mpi::broadcast ( comm_, local_, 0 );

}

      
size_t harp::mpi_image::n_rows ( ) const {
  return local_->n_rows();
}

      
size_t harp::mpi_image::n_cols ( ) const {
  return local_->n_cols();
}


void harp::mpi_image::values ( vector_double & data ) const {
  local_->values ( data );
  return;
}


void harp::mpi_image::inv_variance ( vector_double & invvar ) const {
  local_->inv_variance ( invvar );
  return;
}


boost::property_tree::ptree harp::mpi_image::metadata ( ) const {
  return local_->metadata();
}


void harp::mpi_image::values ( matrix_double & data ) const {
  local_->values ( data );
  return;
}


void harp::mpi_image::inv_variance ( matrix_double & invvar ) const {
  local_->inv_variance ( invvar );
  return;
}


std::string harp::mpi_image::type ( ) const {
  return local_->type();
}






