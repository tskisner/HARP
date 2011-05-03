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


void harp::fits::img_seek ( fitsfile * fp, int hdu ) {
  
  int ret;
  int status = 0;
  int type;
  
  ret = fits_movabs_hdu ( fp, hdu, &type, &status );
  fits::check ( status );
  
  if ( type != IMAGE_HDU ) {
    ostringstream o;
    o << "FITS HDU " << hdu << " is not an image";
    MOAT_THROW( o.str().c_str() );
  }
  
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


void harp::fits::img_read ( fitsfile * fp, size_t frow, size_t fcol, dense_mat_view & data ) {
  
  int ret;
  int status = 0;

  long fpixel[2];
  fpixel[0] = (long)fcol + 1;
  
  long long width = data.size2();
  int anynul;

  for ( size_t i = 0; i < data.size1(); ++i ) {
    fpixel[1] = (long)( frow + i ) + 1;

    ret = fits_read_pix ( fp, TDOUBLE, fpixel, width, 0, &(data( i, 0 )), &anynul, &status );
    fits::check ( status );
  }
  
  return;
}


void img_read_row ( fitsfile * fp, size_t row, data_vec & data ) {
  
  int ret;
  int status = 0;

  long fpixel[2];
  fpixel[0] = 0;
  fpixel[1] = (long)row + 1;
  
  int anynul;

  ret = fits_read_pix ( fp, TDOUBLE, fpixel, (long long)data.size(), 0, &(data[0]), &anynul, &status );
  fits::check ( status );
  
  return;
}


