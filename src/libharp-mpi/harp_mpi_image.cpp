/*
  High Performance Astrophysical Reconstruction and Processing (HARP)
  (c) 2014-2015, The Regents of the University of California, 
  through Lawrence Berkeley National Laboratory.  See top
  level LICENSE file for details.
*/

#include <harp_mpi_internal.hpp>


using namespace std;
using namespace harp;


harp::mpi_image::mpi_image ( boost::mpi::communicator const & comm, std::string const & type, boost::property_tree::ptree const & props ) {

  comm_ = comm;

  int rank = comm.rank();

  plugin_registry & reg = plugin_registry::get();

  if ( rank == 0 ) {
    // instantiate
    local_.reset ( reg.create_image ( type, props ) );
  }
    
  // broadcast to all processes

  mpi_comm_bcast ( comm_, local_, 0 );

}

      
size_t harp::mpi_image::n_rows ( ) const {
  return local_->n_rows();
}

      
size_t harp::mpi_image::n_cols ( ) const {
  return local_->n_cols();
}


void harp::mpi_image::values ( vector_double & data ) const {
  int rank = comm_.rank();
  if ( rank == 0 ) {
    local_->values ( data );
  }
  boost::mpi::broadcast ( comm_, data, 0 );
  return;
}


void harp::mpi_image::inv_variance ( vector_double & invvar ) const {
  int rank = comm_.rank();
  if ( rank == 0 ) {
    local_->inv_variance ( invvar );
  }
  boost::mpi::broadcast ( comm_, invvar, 0 );
  return;
}


void harp::mpi_image::values ( matrix_double & data ) const {
  int rank = comm_.rank();
  if ( rank == 0 ) {
    local_->values ( data );
  }
  boost::mpi::broadcast ( comm_, data, 0 );
  return;
}


void harp::mpi_image::mask ( vector_mask & msk ) const {
  int rank = comm_.rank();
  if ( rank == 0 ) {
    local_->mask ( msk );
  }
  boost::mpi::broadcast ( comm_, msk, 0 );
  return;
}


void harp::mpi_image::mask ( matrix_mask & msk ) const {
  int rank = comm_.rank();
  if ( rank == 0 ) {
    local_->mask ( msk );
  }
  boost::mpi::broadcast ( comm_, msk, 0 );
  return;
}


void harp::mpi_image::inv_variance ( matrix_double & invvar ) const {
  int rank = comm_.rank();
  if ( rank == 0 ) {
    local_->inv_variance ( invvar );
  }
  boost::mpi::broadcast ( comm_, invvar, 0 );
  return;
}


void harp::mpi_image::values ( elem_matrix_local & data ) const {
  vector_double temp;
  values ( temp );
  ublas_to_elem ( temp, data );
  return;
}


void harp::mpi_image::inv_variance ( elem_matrix_local & invvar ) const {
  vector_double temp;
  inv_variance ( temp );
  ublas_to_elem ( temp, invvar );
  return;
}


std::string harp::mpi_image::type ( ) const {
  return local_->type();
}


image_p harp::mpi_image::local ( ) {
  return local_;
}




