// @COPYRIGHT@

#include <harp_internal.hpp>


using namespace std;
using namespace harp;


harp::psf::psf ( string const & format, std::map < std::string, std::string > const & params ) {
  format_ = format;
  params_ = params;
}


void harp::psf::cleanup ( ) {
  // nothing for now
  return;
}


string harp::psf::format ( ) {
  return format_;
}


psf * harp::psf::clone ( ) {
  return create ( format_, params_ );
}


psf * harp::psf::create ( std::string const & format, std::map < std::string, std::string > const & params ) {
  
  #include "harp_psf_formats.cpp"
  
  std::ostringstream o;
  o << "Cannot create psf of unknown format (" << format << ")";
  MOAT_THROW( o.str().c_str() );

  return NULL;
  
}


