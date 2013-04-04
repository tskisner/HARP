// @COPYRIGHT@

#include <harp_internal.hpp>


using namespace std;
using namespace harp;


harp::spec::spec ( boost::property_tree::ptree const & props ) {
  props_ = props;
  format_ = props.get < string > ( "format" );
}


void harp::spec::cleanup ( ) {
  // nothing for now
  return;
}


string harp::spec::format ( ) {
  return format_;
}


spec * harp::spec::clone ( ) {
  return create ( props_ );
}


spec * harp::spec::create ( boost::property_tree::ptree const & props ) {

  string format = props.get < string > ( "format" );

  if ( format == "specter" ) {
    return static_cast < spec * > ( new spec_specter ( props ) );
  }

  if ( format == "sim" ) {
    return static_cast < spec * > ( new spec_sim ( props ) );
  }
  
  std::ostringstream o;
  o << "Cannot create spec of unknown format (" << format << ")";
  HARP_THROW( o.str().c_str() );

  return NULL;
  
}

