// @COPYRIGHT@

#include <harp_internal.hpp>


using namespace std;
using namespace harp;


harp::spec::spec ( string const & format, boost::property_tree::ptree const & props ) {
  format_ = format;
  props_ = props;
}


void harp::spec::cleanup ( ) {
  // nothing for now
  return;
}


string harp::spec::format ( ) {
  return format_;
}


spec * harp::spec::clone ( ) {
  return create ( format_, props_ );
}


spec * harp::spec::create ( std::string const & format, boost::property_tree::ptree const & props ) {
  
  #include "harp_spec_formats.cpp"
  
  std::ostringstream o;
  o << "Cannot create spec of unknown format (" << format << ")";
  MOAT_THROW( o.str().c_str() );

  return NULL;
  
}

