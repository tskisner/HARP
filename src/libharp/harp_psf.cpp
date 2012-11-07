// @COPYRIGHT@

#include <harp_internal.hpp>


using namespace std;
using namespace harp;


harp::psf::psf ( string const & format, boost::property_tree::ptree const & props ) {
  format_ = format;
  props_ = props;
}


void harp::psf::cleanup ( ) {
  // nothing for now
  return;
}


string harp::psf::format ( ) {
  return format_;
}


psf * harp::psf::clone ( ) {
  return create ( format_, props_ );
}


psf * harp::psf::create ( std::string const & format, boost::property_tree::ptree const & props ) {
  
  #include "harp_psf_formats.cpp"
  
  std::ostringstream o;
  o << "Cannot create psf of unknown format (" << format << ")";
  HARP_THROW( o.str().c_str() );

  return NULL;
  
}


