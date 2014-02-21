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
    
  // broadcast to all processes.  Instead of using the boost::mpi packed archive types, we manually
  // serialize to a binary archive and broadcast the bytes.  This is due to some complicated details:
  // 
  // 1.  boost serialization of derived types using a base class pointer requires explicit "registration"
  // of the type using BOOST_CLASS_EXPORT.  This export command must only be found in one compiled source
  // file (the file where the plugin implementation is located in our case).
  // 
  // 2.  The BOOST_CLASS_EXPORT must come after ALL archive headers have been #include'd.  These are included
  // in harp/plugin.hpp, but only include the serial archive types.  The export command explicitly instantiates
  // template specialization of the class serialize function for each archive type.
  //
  // 3.  In order to use boost::mpi packed types, we would have to #include the mpi archive headers in the 
  // plugin implementation file before BOOST_CLASS_EXPORT.  However this would require compiling all
  // plugins with the MPI compiler!
  //
  // since we want to use a single DLL plugin file for both serial and MPI versions of HARP, we are left
  // serializing all types as a vector of bytes and sending that...

  std::ostringstream oss;
  boost::archive::text_oarchive oa(oss);
oa << value;


  boost::archive::binary_oarchive ba ( ofs );
    oa << BOOST_SERIALIZATION_NVP(one);

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




