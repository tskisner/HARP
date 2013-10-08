// @COPYRIGHT@

#include <harp_internal.hpp>


using namespace std;
using namespace harp;


harp::image::image ( boost::property_tree::ptree const & props ) {
  props_ = props;
  format_ = props.get < string > ( "format" );
}


string harp::image::format ( ) {
  return format_;
}


image * harp::image::clone ( ) {
  return create ( props_ );
}


image * harp::image::create ( boost::property_tree::ptree const & props ) {

  string format = props.get < string > ( "format" );

  if ( format == "fits" ) {
    return static_cast < image * > ( new image_fits ( props ) );
  }

  if ( format == "sim" ) {
    return static_cast < image * > ( new image_sim ( props ) );
  }
  
  std::ostringstream o;
  o << "Cannot create image of unknown format (" << format << ")";
  HARP_THROW( o.str().c_str() );

  return NULL;
  
}


