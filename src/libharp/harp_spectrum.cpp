// @COPYRIGHT@

#include <harp_internal.hpp>


using namespace std;
using namespace harp;


harp::spectrum::spectrum ( string const & format, std::map < std::string, std::string > const & params ) {
  format_ = format;
  params_ = params;
}


void harp::spectrum::cleanup ( ) {
  // nothing for now
  return;
}


string harp::spectrum::format ( ) {
  return format_;
}


spectrum * harp::spectrum::clone ( ) {
  return create ( format_, params_ );
}


spectrum * harp::spectrum::create ( std::string const & format, std::map < std::string, std::string > const & params ) {
  
  #include "harp_spectrum_formats.cpp"
  
  std::ostringstream o;
  o << "Cannot create spectrum of unknown format (" << format << ")";
  MOAT_THROW( o.str().c_str() );

  return NULL;
  
}

