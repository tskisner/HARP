/*
  High Performance Astrophysical Reconstruction and Processing (HARP)
  (c) 2014-2015, The Regents of the University of California, 
  through Lawrence Berkeley National Laboratory.  See top
  level LICENSE file for details.
*/

#include <harp_data_internal.hpp>


using namespace std;
using namespace harp;


// instantiate serialization implementations for base classes

BOOST_CLASS_EXPORT(harp::image)
BOOST_CLASS_EXPORT(harp::spec)
BOOST_CLASS_EXPORT(harp::psf)
BOOST_CLASS_EXPORT(harp::object)
BOOST_CLASS_EXPORT(harp::targets)

// utilities

std::string harp::ptree_quote ( const std::string & s ) {
  return "\"" + s + "\"";
}

void harp::ptree_print ( const boost::property_tree::ptree & pt, int level ) {
  const std::string sep ( 2 * level, ' ' );
  BOOST_FOREACH ( const boost::property_tree::ptree::value_type & v, pt ) {
    std::cerr << sep << ptree_quote ( v.first ) << " : " << ptree_quote ( v.second.data() ) << "\n";
    ptree_print ( v.second, level + 1 );
  }
  return;
}

void harp::ptree_print ( const boost::property_tree::ptree & pt ) {
  ptree_print ( pt, 0 );
  return;
}


string const & harp::source_version ( ) {
  // Include the generated file that contains the git revision
  #include "git-version.cpp"
  return harp_revision_key;
}




