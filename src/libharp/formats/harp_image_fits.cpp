// @COPYRIGHT@

#include <harp_internal.hpp>


using namespace std;
using namespace harp;

static const char * image_fits_key_path = "path";
static const char * image_fits_key_signal = "signal";
static const char * image_fits_key_sky = "sky";
static const char * image_fits_key_noise = "noise";
static const char * image_fits_key_rows = "rows";
static const char * image_fits_key_cols = "cols";


harp::image_fits::image_fits ( boost::property_tree::ptree const & props ) : image ( props ) {

  sighdu_ = props.get ( image_fits_key_signal, 1 );

  nsehdu_ = props.get ( image_fits_key_noise, 2 );

  skyhdu_ = props.get ( image_fits_key_sky, 3 );

  path_ = props.get ( image_fits_key_path, "" );

  if ( path_ == "" ) {

    // no path specified- must specify rows / cols

    rows_ = props.get < size_t > ( image_fits_key_rows );

    cols_ = props.get < size_t > ( image_fits_key_cols );

  } else {
    
    // read rows / cols from the FITS header

    fitsfile *fp;
    
    fits::open_read ( fp, path_ );

    fits::img_seek ( fp, sighdu_ );
    
    fits::img_dims ( fp, rows_, cols_ );

    size_t rowcheck;
    size_t colcheck;

    fits::img_seek ( fp, nsehdu_ );
    
    fits::img_dims ( fp, rowcheck, colcheck );

    if ( ( rowcheck != rows_ ) || ( colcheck != cols_ ) ) {
      HARP_THROW( "noise variance dimensions do not match data dimensions" );
    }

  }  
  
}


harp::image_fits::~image_fits ( ) {

}


void harp::image_fits::read ( vector_double & data, vector_double & invvar, std::vector < bool > & sky ) {

  fitsfile *fp;

  fits::open_read ( fp, path_ );

  fits::img_seek ( fp, sighdu_ );
    
  fits::img_read ( fp, data );

  fits::img_seek ( fp, nsehdu_ );
    
  fits::img_read ( fp, invvar );

  fits::bin_seek ( fp, skyhdu_ );

  specter_read_sky ( fp, sky );

  fits::close ( fp );

  return;
}


void harp::image_fits::write ( std::string const & path, vector_double & data, vector_double & invvar, std::vector < bool > & sky ) {

  fitsfile * fp;
    
  fits::create ( fp, path );
  
  fits::img_append ( fp, rows_, cols_ );
  
  fits::img_write ( fp, data );

  fits::img_append ( fp, rows_, cols_ );
  
  fits::img_write ( fp, invvar );

  specter_write_sky ( fp, sky );

  fits::close ( fp );

  return;
}

