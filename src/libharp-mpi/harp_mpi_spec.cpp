/*
  High Performance Astrophysical Reconstruction and Processing (HARP)
  (c) 2014-2015, The Regents of the University of California, 
  through Lawrence Berkeley National Laboratory.  See top
  level LICENSE file for details.
*/

#include <harp_mpi_internal.hpp>


using namespace std;
using namespace harp;


harp::mpi_spec::mpi_spec ( boost::mpi::communicator const & comm, std::string const & type, boost::property_tree::ptree const & props ) {

  comm_ = comm;

  int rank = comm.rank();

  plugin_registry & reg = plugin_registry::get();

  if ( rank == 0 ) {
    // instantiate
    local_.reset ( reg.create_spec ( type, props ) );
  }

  // broadcast to all processes

  mpi_comm_bcast ( comm_, local_, 0 );
  
}


size_t harp::mpi_spec::n_spec ( ) const {
  return local_->n_spec();
}


size_t harp::mpi_spec::n_lambda ( ) const {
  return local_->n_lambda();
}


void harp::mpi_spec::values ( vector_double & data ) const {
  int rank = comm_.rank();
  if ( rank == 0 ) {
    local_->values ( data );
  }
  boost::mpi::broadcast ( comm_, data, 0 );
  return;
}


void harp::mpi_spec::values ( mpi_matrix & data ) const {
  vector_double temp;
  values ( temp );
  ublas_to_elem ( temp, data );
  return;
}


void harp::mpi_spec::lambda ( vector_double & lambda_vals ) const {
  int rank = comm_.rank();
  if ( rank == 0 ) {
    local_->lambda ( lambda_vals );
  }
  boost::mpi::broadcast ( comm_, lambda_vals, 0 );
  return;
}


std::string harp::mpi_spec::type ( ) const {
  return local_->type();
}


spec_p harp::mpi_spec::local ( ) {
  return local_;
}


