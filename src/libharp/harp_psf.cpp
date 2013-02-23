// @COPYRIGHT@

#include <harp_internal.hpp>


using namespace std;
using namespace harp;


harp::psf::psf ( boost::property_tree::ptree const & props ) {
  props_ = props;
  format_ = props.get < string > ( "format" );
}


void harp::psf::cleanup ( ) {
  // nothing for now
  return;
}


string harp::psf::format ( ) {
  return format_;
}


psf * harp::psf::clone ( ) {
  return create ( props_ );
}


psf * harp::psf::create ( boost::property_tree::ptree const & props ) {
  
  string format = props.get < string > ( "format" );

  #include "harp_psf_formats.cpp"
  
  std::ostringstream o;
  o << "Cannot create psf of unknown format (" << format << ")";
  HARP_THROW( o.str().c_str() );

  return NULL;
  
}


