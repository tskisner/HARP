// @COPYRIGHT@

#include <harp_internal.hpp>


using namespace std;
using namespace harp;

static const char * format_boss = "boss";

static const char * boss_image_key_path = "path";
static const char * boss_image_key_signal = "signal";
static const char * boss_image_key_noise = "noise";
static const char * boss_image_key_rows = "rows";
static const char * boss_image_key_cols = "cols";


harp::image_boss::image_boss ( std::map < std::string, std::string > const & params ) : image ( format_boss, params ) {
  
  map < std::string, std::string > :: const_iterator val;
  
  val = params.find( boss_image_key_signal );
  if ( val == params.end() ) {
    sighdu_ = 1;
  } else {
    sighdu_ = atoi ( val->second.c_str() );
  }
  
  val = params.find( boss_image_key_noise );
  if ( val == params.end() ) {
    nsehdu_ = 1;
  } else {
    nsehdu_ = atoi ( val->second.c_str() );
  }
  
  val = params.find( boss_image_key_path );
  
  if ( val == params.end() ) {
    
    // no path specified- must specify rows / cols
    
    val = params.find( boss_image_key_rows );
    if ( val == params.end() ) {
      MOAT_THROW( "boss_image: must specify rows if no path is given" );
    }
    
    rows_ = (size_t) atoll ( val->second.c_str() );
    
    val = params.find( boss_image_key_cols );
    if ( val == params.end() ) {
      MOAT_THROW( "boss_image: must specify cols if no path is given" );
    }
    
    cols_ = (size_t) atoll ( val->second.c_str() );
    
  } else {
    
    path_ = val->second;
    
    // read rows / cols from the FITS header
    
    fitsfile *fp;

    fits::open_read ( fp, path_ );

    fits::img_seek ( fp, sighdu_ );
    
    fits::img_dims ( fp, rows_, cols_ );
    
    fits::close ( fp );
    
  }  
  
}


harp::image_boss::~image_boss ( ) {
  
  cleanup();
  
}


void harp::image_boss::read ( size_t startrow, size_t startcol, harp::dense_rowmat_view & data ) {
  
  fitsfile *fp;

  fits::open_read ( fp, path_ );

  fits::img_seek ( fp, sighdu_ );
  
  fits::img_read ( fp, startrow, startcol, data );
  
  fits::close ( fp );
  
  return;
}


void harp::image_boss::write ( std::string const & path, size_t startrow, size_t startcol, harp::dense_rowmat_view & data ) {
  
  fitsfile *fp;
  
  fits::open_readwrite ( fp, path );
  
  int nh = fits::nhdus ( fp );

  if ( nh < sighdu_ ) {
    while ( nh < sighdu_ ) {
      fits::img_append ( fp, data.size1(), data.size2() );
      ++nh;
    }
  } else {
    fits::img_seek ( fp, sighdu_ );
  }
  
  fits::img_write ( fp, 0, 0, data );
  
  fits::close ( fp );
  
  return;
}


void harp::image_boss::read_noise ( size_t startrow, size_t startcol, harp::dense_rowmat_view & data ) {
  
  fitsfile *fp;

  fits::open_read ( fp, path_ );

  fits::img_seek ( fp, nsehdu_ );
  
  fits::img_read ( fp, startrow, startcol, data );
  
  fits::close ( fp );
  
  return;
}


void harp::image_boss::write_noise ( std::string const & path, size_t startrow, size_t startcol, harp::dense_rowmat_view & data ) {
  
  fitsfile *fp;
  
  fits::open_readwrite ( fp, path );
  
  int nh = fits::nhdus ( fp );

  if ( nh < nsehdu_ ) {
    while ( nh < nsehdu_ ) {
      fits::img_append ( fp, data.size1(), data.size2() );
      ++nh;
    }
  } else {
    fits::img_seek ( fp, nsehdu_ );
  }
  
  fits::img_write ( fp, 0, 0, data );
  
  fits::close ( fp );
  
  return;
}




