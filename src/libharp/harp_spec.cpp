// @COPYRIGHT@

#include <harp_internal.hpp>


using namespace std;
using namespace harp;


harp::spec::spec ( string const & format, std::map < std::string, std::string > const & params ) {
  format_ = format;
  params_ = params;
}


void harp::spec::cleanup ( ) {
  // nothing for now
  return;
}


string harp::spec::format ( ) {
  return format_;
}


spec * harp::spec::clone ( ) {
  return create ( format_, params_ );
}


spec * harp::spec::create ( std::string const & format, std::map < std::string, std::string > const & params ) {
  
  #include "harp_spec_formats.cpp"
  
  std::ostringstream o;
  o << "Cannot create spec of unknown format (" << format << ")";
  MOAT_THROW( o.str().c_str() );

  return NULL;
  
}

