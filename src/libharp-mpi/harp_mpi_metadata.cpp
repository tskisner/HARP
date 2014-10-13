// @COPYRIGHT@

#include <harp_mpi_internal.hpp>


using namespace std;
using namespace harp;


mpi_spec_p harp::mpi_load_spec ( boost::mpi::communicator const & comm, boost::property_tree::ptree const & tree ) {

  boost::property_tree::ptree props;
  string type = tree.get < string > ( "type" );
  props = tree.get_child ( "props" );

  return mpi_spec_p ( new mpi_spec ( comm, type, props ) );
}


mpi_image_p harp::mpi_load_image ( boost::mpi::communicator const & comm, boost::property_tree::ptree const & tree ) {

  boost::property_tree::ptree props;
  string type = tree.get < string > ( "type" );
  props = tree.get_child ( "props" );

  return mpi_image_p ( new mpi_image ( comm, type, props ) );
}


mpi_psf_p harp::mpi_load_psf ( boost::mpi::communicator const & comm, boost::property_tree::ptree const & tree ) {

  boost::property_tree::ptree props;
  string type = tree.get < string > ( "type" );
  props = tree.get_child ( "props" );

  return mpi_psf_p ( new mpi_psf ( comm, type, props ) );
}


mpi_targets_p harp::mpi_load_targets ( boost::mpi::communicator const & comm, boost::property_tree::ptree const & tree ) {

  boost::property_tree::ptree props;
  string type = tree.get < string > ( "type" );
  props = tree.get_child ( "props" );

  return mpi_targets_p ( new mpi_targets ( comm, type, props ) );
}


harp::mpi_group::mpi_group ( boost::mpi::communicator const & comm, boost::property_tree::ptree const & tree ) {
  load ( comm, tree );
}


void harp::mpi_group::load ( boost::mpi::communicator const & comm, boost::property_tree::ptree const & tree ) {

  handle_.reset();

  boost::property_tree::ptree psf_tree = tree.get_child ( "psf" );
  handle_ = mpi_load_psf ( comm, psf_tree );

  imgs_.clear();

  boost::property_tree::ptree::const_iterator v = tree.begin();

  while ( v != tree.end() ) {

    if ( v->first == "image" ) {
      imgs_.push_back ( mpi_load_image ( comm, v->second ) );
    }

    ++v;
  }

  return;
}


std::list < mpi_group_p > harp::mpi_load_groups ( boost::mpi::communicator const & comm, std::string const & path ) {
  // Read JSON into a property tree
  
  boost::property_tree::ptree tree;
  boost::property_tree::json_parser::read_json ( path, tree );

  std::list < mpi_group_p > ret;

  boost::property_tree::ptree::const_iterator v = tree.begin();

  while ( v != tree.end() ) {

    if ( v->first == "group" ) {
      ret.push_back ( mpi_group_p ( new mpi_group ( comm, v->second ) ) );
    }

    ++v;
  }

  return ret;
}


