// @COPYRIGHT@

#include <harp_internal.hpp>


using namespace std;
using namespace harp;


void harp::fits::check ( int status ) {
  if ( status ) {
    char msg[ FLEN_ERRMSG ];
    fits_get_errstatus ( status, msg );

    ostringstream o;
    o << "cfitsio library error: " << msg;
    MOAT_THROW( o.str().c_str() );
  }
  return;
}


void harp::fits::close ( fitsfile * fp ) {
  
  int ret;
  int status = 0;
  
  ret = fits_close_file ( fp, &status );
  fits::check ( status );
  
  return;
}


void harp::fits::open_read ( fitsfile * & fp, string const & path ) {
  
  int ret;
  int status = 0;

  ret = fits_open_file ( &fp, path.c_str(), READONLY, &status );
  fits::check ( status );  
  
  return;
}


void harp::fits::img_seek ( fitsfile * fp, string const & extname ) {
  
  int ret;
  int status = 0;
  
  char extcopy[FLEN_VALUE];
  strncpy ( extcopy, extname.c_str(), FLEN_VALUE );
  
  ret = fits_movnam_hdu ( fp, IMAGE_HDU, extcopy, 0, &status );
  fits::check ( status );
  
  return;
}


void harp::fits::img_dims ( fitsfile * fp, size_t & rows, size_t & cols ) {
  
  int ret;
  int status = 0;
  int naxis;
  
  ret = fits_get_img_dim ( fp, &naxis, &status );
  fits::check ( status );
  
  if ( naxis != 2 ) {
    ostringstream o;
    o << "FITS image has " << naxis << " dimensions instead of 2";
    MOAT_THROW( o.str().c_str() );
  }
  
  long naxes[2];
  
  ret = fits_get_img_size ( fp, 2, naxes, &status );
  fits::check ( status );
  
  cols = naxes[0];
  rows = naxes[1];
  
  return;
}


