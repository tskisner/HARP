// @COPYRIGHT@

#include <harp_internal.hpp>


using namespace std;
using namespace harp;


harp::image::image ( string const & format, boost::property_tree::ptree const & props ) {
  format_ = format;
  props_ = props;
}


void harp::image::cleanup ( ) {
  // nothing for now
  return;
}


string harp::image::format ( ) {
  return format_;
}


image * harp::image::clone ( ) {
  return create ( format_, props_ );
}


image * harp::image::create ( std::string const & format, boost::property_tree::ptree const & props ) {
  
  #include "harp_image_formats.cpp"
  
  std::ostringstream o;
  o << "Cannot create image of unknown format (" << format << ")";
  HARP_THROW( o.str().c_str() );

  return NULL;
  
}


