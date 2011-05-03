// @COPYRIGHT@

#include <harp.hpp>


using namespace std;
using namespace harp;


harp::image::image ( string const & format ) {
  format_ = format;
}


void harp::image::cleanup ( ) {
  // nothing for now
  return;
}


string harp::image::format ( ) {
  return format_;
}


image * harp::image::clone ( ) {
  return create ( format_, rows(), cols(), params_ );
}


image * harp::image::create ( std::string const & format, std::map < std::string, std::string > const & params ) {
  
  image * result = NULL;
  
  #include "harp_image_formats.cpp"
  
  std::ostringstream o;
  o << "Cannot create image of unknown format (" << format << ")";
  MOAT_THROW( o.str().c_str() );

  return result;
  
}


image * harp::image::create ( std::string const & format, std::string const & path, std::map < std::string, std::string > const & params ) {
  
  image * result = NULL;
  
  #include "harp_image-path_formats.cpp"
  
  std::ostringstream o;
  o << "Cannot create image of unknown format (" << format << ")";
  MOAT_THROW( o.str().c_str() );

  return result;
  
}

