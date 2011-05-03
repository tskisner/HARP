// @COPYRIGHT@

#include <harp_internal.hpp>


using namespace std;
using namespace harp;


harp::image::image ( string const & format, std::map < std::string, std::string > const & params ) {
  format_ = format;
  params_ = params;
}


void harp::image::cleanup ( ) {
  // nothing for now
  return;
}


string harp::image::format ( ) {
  return format_;
}


image * harp::image::clone ( ) {
  return create ( format_, params_ );
}


image * harp::image::create ( std::string const & format, std::map < std::string, std::string > const & params ) {
  
  #include "harp_image_formats.cpp"
  
  std::ostringstream o;
  o << "Cannot create image of unknown format (" << format << ")";
  MOAT_THROW( o.str().c_str() );

  return NULL;
  
}


