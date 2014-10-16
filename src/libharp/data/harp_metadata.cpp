/*
  High Performance Astrophysical Reconstruction and Processing (HARP)
  (c) 2014-2015, The Regents of the University of California, 
  through Lawrence Berkeley National Laboratory.  See top
  level LICENSE file for details.
*/

#include <harp_data_internal.hpp>


using namespace std;
using namespace harp;


spec_p harp::load_spec ( boost::property_tree::ptree const & tree ) {
  // Get plugin registry
  plugin_registry & reg = plugin_registry::get();

  boost::property_tree::ptree props;
  string type = tree.get < string > ( "type" );
  props = tree.get_child ( "props" );

  return spec_p ( reg.create_spec ( type, props ) );
}


image_p harp::load_image ( boost::property_tree::ptree const & tree ) {
  // Get plugin registry
  plugin_registry & reg = plugin_registry::get();

  boost::property_tree::ptree props;
  string type = tree.get < string > ( "type" );
  props = tree.get_child ( "props" );

  return image_p ( reg.create_image ( type, props ) );
}


psf_p harp::load_psf ( boost::property_tree::ptree const & tree ) {
  // Get plugin registry
  plugin_registry & reg = plugin_registry::get();

  boost::property_tree::ptree props;
  string type = tree.get < string > ( "type" );
  props = tree.get_child ( "props" );

  return psf_p ( reg.create_psf ( type, props ) );
}


targets_p harp::load_targets ( boost::property_tree::ptree const & tree ) {
  // Get plugin registry
  plugin_registry & reg = plugin_registry::get();

  boost::property_tree::ptree props;
  string type = tree.get < string > ( "type" );
  props = tree.get_child ( "props" );

  return targets_p ( reg.create_targets ( type, props ) );
}


harp::group::group ( boost::property_tree::ptree const & tree ) {
  load ( tree );
}


void harp::group::load ( boost::property_tree::ptree const & tree ) {

  handle_.reset();

  boost::property_tree::ptree psf_tree = tree.get_child ( "psf" );
  handle_ = load_psf ( psf_tree );

  imgs_.clear();

  boost::property_tree::ptree::const_iterator v = tree.begin();

  while ( v != tree.end() ) {

    if ( v->first == "image" ) {
      imgs_.push_back ( load_image ( v->second ) );
    }

    ++v;
  }

  return;
}


std::list < group_p > harp::load_groups ( std::string const & path ) {
  // Read JSON into a property tree
  
  boost::property_tree::ptree tree;
  boost::property_tree::json_parser::read_json ( path, tree );

  std::list < group_p > ret;

  boost::property_tree::ptree::const_iterator v = tree.begin();

  while ( v != tree.end() ) {

    if ( v->first == "group" ) {
      ret.push_back ( group_p ( new group ( v->second ) ) );
    }

    ++v;
  }

  return ret;
}


