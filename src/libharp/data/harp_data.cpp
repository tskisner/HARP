// @COPYRIGHT@

#include <harp_data_internal.hpp>

#include <harp/plugin.hpp>
#ifdef HAVE_BOOST_MPI_HPP
#include <harp/plugin_mpi.hpp>
#endif

using namespace std;
using namespace harp;


// instantiate serialization implementations for base classes

BOOST_CLASS_EXPORT(harp::image)
BOOST_CLASS_EXPORT(harp::spec)
BOOST_CLASS_EXPORT(harp::psf)

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



