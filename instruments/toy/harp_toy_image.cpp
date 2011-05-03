// @COPYRIGHT@

#include <harp_internal.hpp>


using namespace std;
using namespace harp;


static const char * toy_image_key_path = "path";
static const char * toy_image_key_name = "name";
static const char * toy_image_key_rows = "rows";
static const char * toy_image_key_cols = "cols";


harp::image_toy::image_toy ( std::map < std::string, std::string > const & params ) : image ( format_toy, params ) {
  
  map < std::string, std::string > :: const_iterator val;
  
  val = params.find( toy_image_key_path );
  
  if ( val == params.end() ) {
    
    // no path specified- must specify rows / cols
    
    val = params.find( toy_image_key_rows );
    if ( val == params.end() ) {
      MOAT_THROW( "toy_image: must specify rows if no path is given" );
    }
    
    rows_ = (size_t) atoll ( val->second.c_str() );
    
    val = params.find( toy_image_key_cols );
    if ( val == params.end() ) {
      MOAT_THROW( "toy_image: must specify cols if no path is given" );
    }
    
    cols_ = (size_t) atoll ( val->second.c_str() );
    
  } else {
    
    path_ = val->second;
    
    val = params.find( toy_image_key_name );
    if ( val == params.end() ) {
      MOAT_THROW( "toy_image: must specify image extension name in addition to path" );
    }
    
    name_ = val->second;
    
    // read rows / cols from the FITS header
    
    fitsfile *fp;

    fits::open_read ( fp, path_ );

    fits::img_seek ( fp, name_ );
    
    fits::img_dims ( fp, rows_, cols_ );
    
    fits::close ( fp );
    
  }  
  
}


harp::image_toy::~image_toy ( ) {
  
  
}


void harp::image_toy::read ( size_t startrow, size_t startcol, harp::dense_mat_view & data ) {
  
  
  
  return;
}
