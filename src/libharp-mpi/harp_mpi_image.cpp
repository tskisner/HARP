// @COPYRIGHT@

#include <harp_mpi_internal.hpp>


using namespace std;
using namespace harp;


harp::mpi_image::mpi_image ( boost::mpi::communicator const & comm, boost::property_tree::ptree const & props ) {

  comm_ = comm;

  int rank = comm.rank();

  if ( rank == 0 ) {
    // instantiate
    local_.reset ( image::create ( props ) );
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


std::string harp::mpi_image::format ( ) const {
  return local_->format();
}






