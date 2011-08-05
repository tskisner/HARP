// @COPYRIGHT@

#include <harp_internal.hpp>


using namespace std;
using namespace harp;

static const char * format_boss = "boss";

static const char * boss_spectrum_key_path = "path";
static const char * boss_spectrum_key_hdu = "hdu";
static const char * boss_spectrum_key_size = "size";
static const char * boss_spectrum_key_pos = "pos";


harp::spectrum_boss::spectrum_boss ( std::map < std::string, std::string > const & params ) : spectrum ( format_boss, params ) {
  
  map < std::string, std::string > :: const_iterator val;
  
  val = params.find( boss_spectrum_key_pos );
  if ( val == params.end() ) {
    pos_ = 0;
  } else {
    pos_ = (size_t) atoll ( val->second.c_str() );
  }
  
  val = params.find( boss_spectrum_key_hdu );
  if ( val == params.end() ) {
    hdu_ = 1;
  } else {
    hdu_ = atoi ( val->second.c_str() );
  }
    
  val = params.find( boss_spectrum_key_path );
  
  if ( val == params.end() ) {
    
    // no path specified- must specify size
    
    val = params.find( boss_spectrum_key_size );
    if ( val == params.end() ) {
      MOAT_THROW( "boss_spectrum: must specify size if no path is given" );
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


harp::spectrum_boss::~spectrum_boss ( ) {
  
  cleanup();
  
}


void harp::spectrum_boss::read ( harp::data_vec_view & data ) {
  
  fitsfile *fp;

  fits::open_read ( fp, path_ );

  fits::img_seek ( fp, hdu_ );
  
  fits::img_read_row ( fp, pos_, data );
  
  fits::close ( fp );
  
  return;
}


void harp::spectrum_boss::write ( string const & path, harp::data_vec_view & data ) {
  
  fitsfile *fp;

  fits::open_readwrite ( fp, path );

  fits::img_seek ( fp, hdu_ );
  
  fits::img_write_row ( fp, pos_, data );
  
  fits::close ( fp );
  
  return;
}



