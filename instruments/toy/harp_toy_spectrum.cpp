// @COPYRIGHT@

#include <harp_internal.hpp>


using namespace std;
using namespace harp;


static const char * toy_spectrum_key_path = "path";
static const char * toy_spectrum_key_hdu = "hdu";
static const char * toy_spectrum_key_size = "size";
static const char * toy_spectrum_key_pos = "pos";


harp::spectrum_toy::spectrum_toy ( std::map < std::string, std::string > const & params ) : spectrum ( format_toy, params ) {
  
  map < std::string, std::string > :: const_iterator val;
  
  val = params.find( toy_spectrum_key_pos );
  if ( val == params.end() ) {
    pos_ = 0;
  } else {
    pos_ = (size_t) atoll ( val->second.c_str() );
  }
  
  val = params.find( toy_spectrum_key_hdu );
  if ( val == params.end() ) {
    hdu_ = 1;
  } else {
    hdu_ = atoi ( val->second.c_str() );
  }
    
  val = params.find( toy_spectrum_key_path );
  
  if ( val == params.end() ) {
    
    // no path specified- must specify size
    
    val = params.find( toy_spectrum_key_size );
    if ( val == params.end() ) {
      MOAT_THROW( "toy_spectrum: must specify size if no path is given" );
    }
    
    size_ = (size_t) atoll ( val->second.c_str() );
    
  } else {
    
    path_ = val->second;
    
    // read size from the FITS header
    
    fitsfile *fp;

    fits::open_read ( fp, path_ );

    fits::img_seek ( fp, hdu_ );
    
    size_t rows, cols;
    
    fits::img_dims ( fp, rows, cols );
    
    size_ = cols;
    
    fits::close ( fp );
    
  }  
  
}


harp::spectrum_toy::~spectrum_toy ( ) {
  
  
}


void harp::spectrum_toy::read ( harp::data_vec & data ) {
  
  fitsfile *fp;

  fits::open_read ( fp, path_ );

  fits::img_seek ( fp, hdu_ );
  
  fits::img_read_row ( fp, pos_, data );
  
  fits::close ( fp );
  
  return;
}

