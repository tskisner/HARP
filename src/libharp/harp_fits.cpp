// @COPYRIGHT@

#include <harp_internal.hpp>

extern "C" {
  #include <sys/stat.h>
}

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


void harp::fits::open_readwrite ( fitsfile * & fp, string const & path ) {
  
  int ret;
  int status = 0;
  
  struct stat statbuf;
  int statret;
  
  statret = stat ( path.c_str(), &statbuf );

  if ( statret ) {
    // create file
    
    ret = fits_create_file ( &fp, path.c_str(), &status );
    fits::check ( status );
  } else {
    // just open it read-write
    
    ret = fits_open_file ( &fp, path.c_str(), READWRITE, &status );
    fits::check ( status );
  }
  
  return;
}


int harp::fits::nhdus ( fitsfile * fp ) {
  int nhdu;
  
  int ret;
  int status = 0;
  
  ret = fits_get_num_hdus ( fp, &nhdu, &status );
  fits::check ( status );
  
  return nhdu;
}


string harp::fits::key_string ( fitsfile * fp, std::string keyname ) {
  int ret;
  int status = 0;
  
  char keycopy[FLEN_VALUE];
  strncpy ( keycopy, keyname.c_str(), FLEN_VALUE );
  
  char value[FLEN_VALUE];
  char comment[FLEN_VALUE];
  
  ret = fits_read_key ( fp, TSTRING, keycopy, value, comment, &status );
  fits::check ( status );
  
  return string(value);
}


long harp::fits::key_long ( fitsfile * fp, std::string keyname ) {
  int ret;
  int status = 0;
  
  char keycopy[FLEN_VALUE];
  strncpy ( keycopy, keyname.c_str(), FLEN_VALUE );
  
  long value;
  char comment[FLEN_VALUE];
  
  ret = fits_read_key ( fp, TLONG, keycopy, &value, comment, &status );
  fits::check ( status );
  
  return value;
}


double harp::fits::key_double ( fitsfile * fp, std::string keyname ) {
  int ret;
  int status = 0;
  
  char keycopy[FLEN_VALUE];
  strncpy ( keycopy, keyname.c_str(), FLEN_VALUE );
  
  double value;
  char comment[FLEN_VALUE];
  
  ret = fits_read_key ( fp, TDOUBLE, keycopy, &value, comment, &status );
  fits::check ( status );
  
  return value;
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
      ret = fits_read_key ( fp, TSTRING, keycopy, valcheck, comment, &status );
      //cerr << "key compare " << keycopy << ": " << valcheck << " =? " << valcopy << endl;
      if ( status == 0 ) {
        // keyword exists
        if ( strncmp ( valcheck, valcopy, strlen ( valcopy ) ) == 0 ) {
          // a match!
          return hdu;
        }
      }
      status = 0;
    }
  }
  
  return -1;
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


void harp::fits::img_append ( fitsfile * fp, size_t rows, size_t cols ) {
  
  int ret;
  int status = 0;
  
  long naxes[2];
  naxes[0] = cols;
  naxes[1] = rows;
  
  ret = fits_create_img ( fp, -64, 2, naxes, &status );
  fits::check ( status );
  
  return;
}


void harp::fits::img_write ( fitsfile * fp, size_t frow, size_t fcol, mat_denserow & data ) {
  
  int ret;
  int status = 0;

  long fpixel[2];
  long lpixel[2];
  
  fpixel[0] = (long)fcol + 1;
  
  long width = data.size2();
  
  lpixel[0] = (fpixel[0] + width - 1) + 1;
  
  double * buffer = moat::double_alloc ( width );
  
  //cerr << "FITS writing image from columns " << fpixel[0]-1 << " to " << lpixel[0]-1 << endl;
  
  for ( size_t i = 0; i < data.size1(); ++i ) {
    fpixel[1] = (long)( frow + i ) + 1;
    lpixel[1] = fpixel[1];
    
    //cerr << "FITS writing row " << i << "(" << fpixel[1] << ")" << endl;

    for ( size_t j = 0; j < data.size2(); ++j ) {
      buffer[j] = data( i, j );
    }

    ret = fits_write_subset ( fp, TDOUBLE, fpixel, lpixel, buffer, &status );
    fits::check ( status );
  }
  
  free ( buffer );
  
  //cerr << "write complete" << endl;
  
  return;
}


void harp::fits::img_write ( fitsfile * fp, vec_dense & data ) {
  
  int ret;
  int status = 0;

  long fpixel[2];
  
  fpixel[0] = 1;
  fpixel[1] = 1;

  ret = fits_write_pix ( fp, TDOUBLE, fpixel, data.size(), &(data[0]), &status );
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


void harp::fits::img_read ( fitsfile * fp, size_t frow, size_t fcol, mat_denserow & data ) {
  
  int ret;
  int status = 0;

  long fpixel[2];
  long lpixel[2];
  long inc[2] = {1, 1};
  
  fpixel[0] = (long)fcol + 1;
  
  long width = data.size2();
  
  lpixel[0] = fpixel[0] + width - 1;
  
  int anynul;
  
  double * buffer = moat::double_alloc ( width );
  
  //cerr << "FITS reading image from columns " << fpixel[0]-1 << " to " << lpixel[0]-1 << endl;

  for ( size_t i = 0; i < data.size1(); ++i ) {
    fpixel[1] = (long)( frow + i ) + 1;
    lpixel[1] = fpixel[1];
    
    //cerr << "FITS reading row " << i << "(" << fpixel[1] << ")" << endl;
    
    ret = fits_read_subset ( fp, TDOUBLE, fpixel, lpixel, inc, 0, buffer, &anynul, &status );
    fits::check ( status );
    
    for ( size_t j = 0; j < data.size2(); ++j ) {
      data( i, j ) = buffer[j];
    }
  }
  
  free ( buffer );
  
  //cerr << "read complete" << endl;
  
  return;
}


void harp::fits::img_read ( fitsfile * fp, vec_dense & data ) {
  
  int ret;
  int status = 0;
  int anynul;

  long fpixel[2];

  fpixel[0] = 1;
  fpixel[1] = 1;

  ret = fits_read_pix ( fp, TDOUBLE, fpixel, data.size(), NULL, &(data[0]), &anynul, &status );
  fits::check ( status );
  
  return;
}


void harp::fits::img_read_row ( fitsfile * fp, size_t row, vec_dense & data ) {
  
  int ret;
  int status = 0;

  long fpixel[2];
  long lpixel[2];
  long inc[2] = {1, 1};
  
  fpixel[0] = 1;
  fpixel[1] = (long)row + 1;
  
  lpixel[0] = (long)data.size();
  lpixel[1] = fpixel[1];
  
  int anynul;
  
  double * buffer = moat::double_alloc ( data.size() );
  
  //cerr << "FITS reading image row " << fpixel[1]-1 << " from columns " << fpixel[0]-1 << " to " << lpixel[0]-1 << endl;
  
  ret = fits_read_subset ( fp, TDOUBLE, fpixel, lpixel, inc, 0, buffer, &anynul, &status );
  fits::check ( status );
  
  for ( size_t j = 0; j < data.size(); ++j ) {
    data( j ) = buffer[j];
  }
  
  free ( buffer );
  
  //cerr << "read row complete" << endl;
  
  return;
}


void harp::fits::img_read_row_int ( fitsfile * fp, size_t row, vec_int & data ) {
  
  int ret;
  int status = 0;

  long fpixel[2];
  long lpixel[2];
  long inc[2] = {1, 1};
  
  fpixel[0] = 1;
  fpixel[1] = (long)row + 1;
  
  lpixel[0] = (long)data.size();
  lpixel[1] = fpixel[1];
  
  int anynul;

  int * buffer = moat::int_alloc ( data.size() );
  
  //cerr << "FITS reading image int row " << fpixel[1]-1 << " from columns " << fpixel[0]-1 << " to " << lpixel[0]-1 << endl;

  ret = fits_read_subset ( fp, TINT, fpixel, lpixel, inc, 0, buffer, &anynul, &status );
  fits::check ( status );
  
  for ( size_t j = 0; j < data.size(); ++j ) {
    data( j ) = buffer[j];
  }
  
  free ( buffer );
  
  //cerr << "read int row complete" << endl;
  
  return;
}


void harp::fits::img_write_row ( fitsfile * fp, size_t row, vec_dense & data ) {
  
  int ret;
  int status = 0;

  long fpixel[2];
  long lpixel[2];
  long inc[2] = {1, 1};
  
  fpixel[0] = 1;
  fpixel[1] = (long)row + 1;
  
  lpixel[0] = (long)data.size();
  lpixel[1] = fpixel[1];
  
  double * buffer = moat::double_alloc ( data.size() );
  
  for ( size_t j = 0; j < data.size(); ++j ) {
    buffer[j] = data( j );
  }
  
  //cerr << "FITS reading image row " << fpixel[1]-1 << " from columns " << fpixel[0]-1 << " to " << lpixel[0]-1 << endl;
  
  ret = fits_write_subset ( fp, TDOUBLE, fpixel, lpixel, buffer, &status );
  fits::check ( status );
  
  free ( buffer );
  
  //cerr << "read row complete" << endl;
  
  return;
}


int harp::fits::bin_seek ( fitsfile * fp, string const & keyname, string const & keyval ) {
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
    
    if ( type == BINARY_TBL ) {
      ret = fits_read_key ( fp, TSTRING, keycopy, valcheck, comment, &status );
      //cerr << "key compare " << keycopy << ": " << valcheck << " =? " << valcopy << endl;
      if ( status == 0 ) {
        // keyword exists
        if ( strncmp ( valcheck, valcopy, strlen ( valcopy ) ) == 0 ) {
          // a match!
          return hdu;
        }
      }
      status = 0;
    }
  }
  
  return -1;
  
}


vector < int > harp::fits::bin_columns ( fitsfile * fp, vector < string > & names ) {
  
  int ret;
  int status = 0;
  
  vector < int > cols;
  
  vector < string > :: iterator nit;
  char namecopy[FLEN_VALUE];
  int id;
  
  for ( nit = names.begin(); nit != names.end(); ++nit ) {
    strncpy ( namecopy, nit->c_str(), FLEN_VALUE );
    ret = fits_get_colnum ( fp, CASEINSEN, namecopy, &id, &status );
    fits::check ( status );
    cols.push_back ( id-1 );
  }
  
  return cols;
} 


void harp::fits::bin_seek ( fitsfile * fp, int hdu ) {
  
  int ret;
  int status = 0;
  int type;
  
  ret = fits_movabs_hdu ( fp, hdu, &type, &status );
  fits::check ( status );
  
  if ( type != BINARY_TBL ) {
    ostringstream o;
    o << "FITS HDU " << hdu << " is not a binary table";
    MOAT_THROW( o.str().c_str() );
  }
  
  return;
}


void harp::fits::bin_read ( fitsfile * fp, size_t firstrow, size_t lastrow, vector < int > & columns, vector < vec_dense > & data ) {
  
  int ret;
  int status = 0;
  long offset = (long)firstrow;
  long nread = (long)lastrow - offset + 1;
  
  // find optimal rowsize
  
  long optimal;  
  ret = fits_get_rowsize ( fp, &optimal, &status );
  fits::check ( status );
  
  optimal--;  // decrement by one, just to be safe
  if (optimal < 1) {
    optimal = 1;
  }
  
  // get table dimensions
  
  long nrows;
  int tfields;
  char fitsval[FLEN_VALUE];
  long pcount;
  
  ret = fits_read_btblhdr ( fp, 100, &nrows, &tfields, NULL, NULL, NULL, fitsval, &pcount, &status );
  fits::check ( status );
  
  if ( offset + nread > nrows ) {
    MOAT_THROW( "binary read range is beyond end of table" );
  } 
  
  // check that column numbers are in range
  
  vector < int > :: iterator it;
  
  int cur = 0;
  for ( it = columns.begin(); it != columns.end(); ++it ) {
    if ( ( (*it) >= tfields ) || ( (*it) < 0 ) ) {
      ostringstream o;
      o << "cannot read (zero-based) column " << (*it) << " from binary table with " << tfields << " columns";
      MOAT_THROW( o.str().c_str() );
    }
    data[ cur ].resize ( nread );
    ++cur;
  }
    
  // read data in a buffered way
  
  int anynul;
  long n = optimal;
  long dataoffset = 0;

  while ( n == optimal ) {
    if ( dataoffset + optimal > nread ) {
      n = nrows - dataoffset;
    }
    
    cur = 0;
    for ( it = columns.begin(); it != columns.end(); ++it ) {
      ret = fits_read_col_dbl ( fp, (*it) + 1, offset + 1, 1, n, 0, &((data[cur])[dataoffset]), &anynul, &status );
      fits::check ( status );
      //cerr << "fits::bin_read col " << (*it)+1 << " rows " << offset+1 << "-" << offset+1+n << " = " << (data[cur])[dataoffset] << ", " << (data[cur])[dataoffset+1] << ", " << (data[cur])[dataoffset+2] << endl;
      ++cur;
    }
    
    offset += optimal;
    dataoffset += optimal;
  }
  
  return;
}


void harp::fits::bin_write ( fitsfile * fp, size_t firstrow, size_t lastrow, std::vector < int > & columns, std::vector < vec_dense > & data ) {
  
  MOAT_THROW( "binary FITS table writing not yet implemented" );
  
  return;
}


void harp::fits::test ( string const & datadir ) {
  
  cerr << "Testing FITS operations..." << endl;
  
  string imgfile = datadir + "/" + "fits_test_img.fits.out";
  
  size_t rows = 100;
  size_t cols = 500;
  
  mat_denserow data ( rows, cols );
  
  for ( size_t i = 0; i < rows; ++i ) {
    for ( size_t j = 0; j < cols; ++j ) {
      data( i, j ) = (double)(i * cols + j);
    }
  }

  
  fitsfile * fp;
  
  fits::open_readwrite ( fp, imgfile );
  
  fits::img_append ( fp, rows, cols );
  
  fits::img_write ( fp, 0, 0, data );

  fits::close ( fp );

  
  fits::open_read ( fp, imgfile );

  fits::img_seek ( fp, 1 );
  
  size_t checkrows;
  size_t checkcols;
  
  fits::img_dims ( fp, checkrows, checkcols );

  if ( ( checkrows != rows ) || ( checkcols != cols ) ) {
    cerr << "  (FAILED): img dimensions wrong (" << checkrows << "x" << checkcols << ") != (" << rows << "x" << cols << ")" << endl;
    exit(1);
  }

  mat_denserow checkdata ( rows, cols );

  fits::img_read ( fp, 0, 0, checkdata );
  
  for ( size_t i = 0; i < rows; ++i ) {
    for ( size_t j = 0; j < cols; ++j ) {
      if ( checkdata( i, j ) != data( i, j ) ) {
        cerr << "  (FAILED): img element (" << i << ", " << j << ") has wrong value (" << checkdata( i, j ) << " != " << data( i, j ) << ")" << endl;
        exit(1);
      }
    }
  }
  
  vec_dense checkrow ( cols );
  vec_int checkrowint ( cols );
  
  for ( size_t i = 0; i < rows; ++i ) {
    
    fits::img_read_row ( fp, i, checkrow );
    
    for ( size_t j = 0; j < cols; ++j ) {
      if ( checkrow( j ) != data( i, j ) ) {
        cerr << "  (FAILED): img row " << i << ", element " << j << " has wrong value (" << checkrow( j ) << " != " << data( i, j ) << ")" << endl;
        exit(1);
      }
    }
  }
  
  for ( size_t i = 0; i < rows; ++i ) {
    
    fits::img_read_row_int ( fp, i, checkrowint );
    
    for ( size_t j = 0; j < cols; ++j ) {
      if ( checkrowint( j ) != (int)data( i, j ) ) {
        cerr << "  (FAILED): img INT row " << i << ", element " << j << " has wrong value (" << checkrowint( j ) << " != " << (int)data( i, j ) << ")" << endl;
        exit(1);
      }
    }
  }

  fits::close ( fp );
  
  
  cerr << "  (PASSED)" << endl;
  
  return;
}


