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


int harp::fits::img_seek ( fitsfile * fp, string const & extname ) {
  int hdu;
  
  int ret;
  int status = 0;
  
  char extcopy[FLEN_VALUE];
  strncpy ( extcopy, extname.c_str(), FLEN_VALUE );
  
  ret = fits_movnam_hdu ( fp, IMAGE_HDU, extcopy, 0, &status );
  fits::check ( status );
  
  ret = fits_get_hdu_num ( fp, &hdu );
  
  return hdu;
}


int harp::fits::img_seek ( fitsfile * fp, std::string const & keyname, std::string const & keyval ) {
  int hdu;
  
  int ret;
  int status = 0;
  
  char keycopy[FLEN_VALUE];
  strncpy ( keycopy, keyname.c_str(), FLEN_VALUE );
  
  char valcopy[FLEN_VALUE];
  strncpy ( valcopy, keyval.c_str(), FLEN_VALUE );
  
  char valcheck[FLEN_VALUE];
  char comment[FLEN_VALUE];
  
  int nhdu;
  
  ret = fits_get_num_hdus ( fp, &nhdu, &status );
  fits::check ( status );
  
  int type;
  
  for ( int i = 0; i < nhdu; ++i ) {
    hdu = 1 + i;
    
    ret = fits_movabs_hdu ( fp, hdu, &type, &status );
    fits::check ( status );
    
    if ( type == IMAGE_HDU ) {
      ret = fits_read_keyword ( fp, keycopy, valcheck, comment, &status );
      if ( status == 0 ) {
        // keyword exists
        if ( strncmp ( valcheck, valcopy, FLEN_VALUE ) == 0 ) {
          // a match!
          break;
        }
      }
      status = 0;
    }
  }
  
  return hdu;
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
  fpixel[1] = (long)fcol + 1;
  
  long long width = data.size2();
  int anynul;

  for ( size_t i = 0; i < data.size1(); ++i ) {
    fpixel[0] = (long)( frow + i ) + 1;

    ret = fits_read_pix ( fp, TDOUBLE, fpixel, width, 0, &(data( i, 0 )), &anynul, &status );
    fits::check ( status );
  }
  
  return;
}


void harp::fits::img_read_row ( fitsfile * fp, size_t row, data_vec & data ) {
  
  int ret;
  int status = 0;

  long fpixel[2];
  fpixel[1] = 1;
  fpixel[0] = (long)row + 1;
  
  int anynul;

  ret = fits_read_pix ( fp, TDOUBLE, fpixel, (long long)data.size(), 0, &(data[0]), &anynul, &status );
  fits::check ( status );
  
  return;
}


void harp::fits::img_read_row_int ( fitsfile * fp, size_t row, int_vec & data ) {
  
  int ret;
  int status = 0;

  long fpixel[2];
  fpixel[1] = 1;
  fpixel[0] = (long)row + 1;
  
  int anynul;

  ret = fits_read_pix ( fp, TINT, fpixel, (long long)data.size(), 0, &(data[0]), &anynul, &status );
  fits::check ( status );
  
  return;
}


