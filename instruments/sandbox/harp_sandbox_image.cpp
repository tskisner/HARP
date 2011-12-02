// @COPYRIGHT@

#include <harp_internal.hpp>


using namespace std;
using namespace harp;

static const char * format_sandbox = "sandbox";

static const char * sandbox_image_key_path = "path";
static const char * sandbox_image_key_signal = "signal";
static const char * sandbox_image_key_noise = "noise";
static const char * sandbox_image_key_rows = "rows";
static const char * sandbox_image_key_cols = "cols";


harp::image_sandbox::image_sandbox ( boost::ptree const & props ) : image ( format_sandbox, props ) {

  sighdu_ = props.get ( sandbox_image_key_signal, 1 );

  nsehdu_ = props.get ( sandbox_image_key_noise, 2 );

  path_ = props.get ( sandbox_image_key_path, "" );

  if ( path_ == "" ) {

    // no path specified- must specify rows / cols

    rows_ = props.get < size_t > ( sandbox_image_key_rows );

    cols_ = props.get < size_t > ( sandbox_image_key_cols );

  } else {
    
    // read rows / cols from the FITS header
    
    fitsfile *fp;

    fits::open_read ( fp, path_ );

    fits::img_seek ( fp, sighdu_ );
    
    fits::img_dims ( fp, rows_, cols_ );
    
    fits::close ( fp );
    
  }  
  
}


harp::image_sandbox::~image_sandbox ( ) {
  
  cleanup();
  
}


boost::ptree harp::image_sandbox::serialize ( ) {
  boost::ptree ret;

  ret.put ( "format", image::format() );

  if ( sighdu_ != 1 ) {
    ret.put ( sandbox_image_key_signal, sighdu_ );
  }

  if ( nsehdu_ != 2 ) {
    ret.put ( sandbox_image_key_noise, nsehdu_ );
  }

  if ( path_ == "" ) {
    ret.put ( sandbox_image_key_rows, rows_ );
    ret.put ( sandbox_image_key_cols, cols_ );
  } else {
    ret.put ( sandbox_image_key_path, path_ );
  }

  return ret;
}


void harp::image_sandbox::read ( size_t startrow, size_t startcol, mat_denserow & data ) {
  
  fitsfile *fp;

  fits::open_read ( fp, path_ );

  fits::img_seek ( fp, sighdu_ );
  
  fits::img_read ( fp, startrow, startcol, data );
  
  fits::close ( fp );
  
  return;
}


void harp::image_sandbox::read ( vec_dense & data ) {
  
  fitsfile *fp;

  fits::open_read ( fp, path_ );

  fits::img_seek ( fp, sighdu_ );
  
  fits::img_read ( fp, data );
  
  fits::close ( fp );
  
  return;
}


void harp::image_sandbox::write ( std::string const & path, size_t startrow, size_t startcol, mat_denserow & data ) {
  
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


void harp::image_sandbox::write ( std::string const & path, vec_dense & data ) {
  
  fitsfile *fp;
  
  fits::open_readwrite ( fp, path );
  
  int nh = fits::nhdus ( fp );

  if ( nh < sighdu_ ) {
    while ( nh < sighdu_ ) {
      fits::img_append ( fp, rows_, cols_ );
      ++nh;
    }
  } else {
    fits::img_seek ( fp, sighdu_ );
  }
  
  fits::img_write ( fp, data );
  
  fits::close ( fp );
  
  return;
}


void harp::image_sandbox::read_noise ( size_t startrow, size_t startcol, mat_denserow & data ) {
  
  fitsfile *fp;

  fits::open_read ( fp, path_ );

  fits::img_seek ( fp, nsehdu_ );
  
  fits::img_read ( fp, startrow, startcol, data );
  
  fits::close ( fp );
  
  return;
}


void harp::image_sandbox::read_noise ( vec_dense & data ) {
  
  fitsfile *fp;

  fits::open_read ( fp, path_ );

  fits::img_seek ( fp, nsehdu_ );
  
  fits::img_read ( fp, data );
  
  fits::close ( fp );
  
  return;
}


void harp::image_sandbox::write_noise ( std::string const & path, size_t startrow, size_t startcol, mat_denserow & data ) {
  
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


void harp::image_sandbox::write_noise ( std::string const & path, vec_dense & data ) {
  
  fitsfile *fp;
  
  fits::open_readwrite ( fp, path );
  
  int nh = fits::nhdus ( fp );

  if ( nh < nsehdu_ ) {
    while ( nh < nsehdu_ ) {
      fits::img_append ( fp, rows_, cols_ );
      ++nh;
    }
  } else {
    fits::img_seek ( fp, nsehdu_ );
  }
  
  fits::img_write ( fp, data );
  
  fits::close ( fp );
  
  return;
}



