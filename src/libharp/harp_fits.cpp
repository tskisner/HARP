// @COPYRIGHT@

#include <harp_internal.hpp>

#include <cstdio>

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
    HARP_THROW( o.str().c_str() );
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

  if ( statret != 0 ) {
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


void harp::fits::create ( fitsfile * & fp, string const & path ) {
  
  int ret;
  int status = 0;
  
  struct stat statbuf;
  int statret;
  
  statret = stat ( path.c_str(), &statbuf );

  if ( statret == 0 ) {
    // delete file
    ret = remove ( path.c_str() );
  }

  ret = fits_create_file ( &fp, path.c_str(), &status );
  fits::check ( status );
  
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


void harp::fits::read_key ( fitsfile * fp, std::string const & keyname, std::string & keyval ) {
  int ret;
  int status = 0;
  
  char keycopy[FLEN_VALUE];
  strncpy ( keycopy, keyname.c_str(), FLEN_VALUE );
  
  char value[FLEN_VALUE];
  char comment[FLEN_VALUE];
  
  ret = fits_read_key ( fp, TSTRING, keycopy, value, comment, &status );
  fits::check ( status );

  keyval = value;
  
  return;
}


void harp::fits::read_key ( fitsfile * fp, std::string const & keyname, long & keyval ) {
  int ret;
  int status = 0;
  
  char keycopy[FLEN_VALUE];
  strncpy ( keycopy, keyname.c_str(), FLEN_VALUE );
  
  long value;
  char comment[FLEN_VALUE];
  
  ret = fits_read_key ( fp, TLONG, keycopy, &value, comment, &status );
  fits::check ( status );

  keyval = value;
  
  return;
}


void harp::fits::read_key ( fitsfile * fp, std::string const & keyname, double & keyval ) {
  int ret;
  int status = 0;
  
  char keycopy[FLEN_VALUE];
  strncpy ( keycopy, keyname.c_str(), FLEN_VALUE );
  
  double value;
  char comment[FLEN_VALUE];
  
  ret = fits_read_key ( fp, TDOUBLE, keycopy, &value, comment, &status );
  fits::check ( status );
  
  keyval = value;

  return;
}


void harp::fits::write_key ( fitsfile * fp, std::string const & keyname, std::string const & keyval, std::string const & keycom ) {
  int ret;
  int status = 0;
  
  char keycopy[FLEN_VALUE];
  strncpy ( keycopy, keyname.c_str(), FLEN_VALUE );
  
  char comment[FLEN_VALUE];
  strncpy ( comment, keycom.c_str(), FLEN_VALUE );

  char value[FLEN_VALUE];
  strncpy ( value, keyval.c_str(), FLEN_VALUE );  

  ret = fits_update_key ( fp, TSTRING, keycopy, value, comment, &status );
  fits::check ( status );

  return;
}


void harp::fits::write_key ( fitsfile * fp, std::string const & keyname, long const & keyval, std::string const & keycom ) {
  int ret;
  int status = 0;
  
  char keycopy[FLEN_VALUE];
  strncpy ( keycopy, keyname.c_str(), FLEN_VALUE );

  char comment[FLEN_VALUE];
  strncpy ( comment, keycom.c_str(), FLEN_VALUE );
  
  long value = keyval;
  
  ret = fits_update_key ( fp, TLONG, keycopy, &value, comment, &status );
  fits::check ( status );
  
  return;
}


void harp::fits::write_key ( fitsfile * fp, std::string const & keyname, double const & keyval, std::string const & keycom ) {
  int ret;
  int status = 0;
  
  char keycopy[FLEN_VALUE];
  strncpy ( keycopy, keyname.c_str(), FLEN_VALUE );

  char comment[FLEN_VALUE];
  strncpy ( comment, keycom.c_str(), FLEN_VALUE );
  
  double value = keyval;
  
  ret = fits_update_key ( fp, TDOUBLE, keycopy, &value, comment, &status );
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
    //cerr << "image seek, looking for " << valcopy << endl;
    
    ret = fits_movabs_hdu ( fp, hdu, &type, &status );
    fits::check ( status );
    
    if ( type == IMAGE_HDU ) {
      ret = fits_read_key ( fp, TSTRING, keycopy, valcheck, comment, &status );
      //cerr << "key compare " << keycopy << ": " << valcheck << " =? " << valcopy << endl;
      if ( status == 0 ) {
        // keyword exists
        if ( strncasecmp ( valcheck, valcopy, strlen ( valcopy ) ) == 0 ) {
          // a match!
          //cerr << "  match! hdu = " << hdu << endl;
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
    HARP_THROW( o.str().c_str() );
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


void harp::fits::img_write ( fitsfile * fp, matrix_local & data ) {
  
  int ret;
  int status = 0;

  long fpixel[2];
  
  fpixel[0] = 1;
  fpixel[1] = 1;

  long npix = data.Height() * data.Width();

  ret = fits_write_pix ( fp, TDOUBLE, fpixel, npix, data.Buffer(), &status );
  fits::check ( status );
  
  return;
}


void harp::fits::img_dims ( fitsfile * fp, size_t & rows, size_t & cols ) {
  
  int ret;
  int status = 0;
  int naxis;
  
  ret = fits_get_img_dim ( fp, &naxis, &status );
  fits::check ( status );
  
  if ( (naxis > 2) || (naxis <= 0) ) {
    ostringstream o;
    o << "FITS image has " << naxis << " dimensions instead of 1 or 2";
    HARP_THROW( o.str().c_str() );
  }
  
  long naxes[2];
  
  ret = fits_get_img_size ( fp, naxis, naxes, &status );
  fits::check ( status );

  cols = naxes[0];
  if ( naxis == 1 ) {
    rows = 1;
  } else {
    rows = naxes[1];
  }
  
  return;
}


void harp::fits::img_read ( fitsfile * fp, matrix_local & data ) {
  
  int ret;
  int status = 0;
  int anynul;

  int naxis;

  ret = fits_get_img_dim ( fp, &naxis, &status );
  fits::check ( status );
  
  if ( (naxis > 2) || (naxis <= 0) ) {
    ostringstream o;
    o << "FITS image has " << naxis << " dimensions instead of 1 or 2";
    HARP_THROW( o.str().c_str() );
  }
  
  long naxes[2];
  
  ret = fits_get_img_size ( fp, naxis, naxes, &status );
  fits::check ( status );

  long nelem = 1;
  for ( int i = 0; i < naxis; ++i ) {
    nelem *= naxes[i];
  }

  if ( data.MemorySize() < nelem ) {
    ostringstream o;
    o << "FITS image has " << nelem << " elements, but matrix has space for only " << data.MemorySize();
    HARP_THROW( o.str().c_str() );
  }

  long fpixel[2];

  fpixel[0] = 1;
  fpixel[1] = 1;

  ret = fits_read_pix ( fp, TDOUBLE, fpixel, nelem, NULL, data.Buffer(), &anynul, &status );
  fits::check ( status );
  
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
        if ( strncasecmp ( valcheck, valcopy, strlen ( valcopy ) ) == 0 ) {
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
    HARP_THROW( o.str().c_str() );
  }
  
  return;
}


void harp::fits::bin_read ( fitsfile * fp, size_t firstrow, size_t lastrow, vector < int > & columns, vector < matrix_local > & data ) {
  
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
    HARP_THROW( "binary read range is beyond end of table" );
  } 
  
  // check that column numbers are in range
  
  vector < int > :: iterator it;
  
  int cur = 0;
  for ( it = columns.begin(); it != columns.end(); ++it ) {
    if ( ( (*it) >= tfields ) || ( (*it) < 0 ) ) {
      ostringstream o;
      o << "cannot read (zero-based) column " << (*it) << " from binary table with " << tfields << " columns";
      HARP_THROW( o.str().c_str() );
    }
    data[ cur ].ResizeTo ( nread, 1 );
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
      ret = fits_read_col_dbl ( fp, (*it) + 1, offset + 1, 1, n, 0, &(( data[cur].Buffer() )[dataoffset]), &anynul, &status );
      fits::check ( status );
      //cerr << "fits::bin_read col " << (*it)+1 << " rows " << offset+1 << "-" << offset+1+n << " = " << (data[cur])[dataoffset] << ", " << (data[cur])[dataoffset+1] << ", " << (data[cur])[dataoffset+2] << endl;
      ++cur;
    }
    
    offset += optimal;
    dataoffset += optimal;
  }
  
  return;
}


void harp::fits::bin_read_strings ( fitsfile * fp, size_t firstrow, size_t lastrow, int col, vector < string > & data ) {
  
  int ret;
  int status = 0;
  long offset = (long)firstrow;
  long nread = (long)lastrow - offset + 1;
  
  // get table dimensions
  
  long nrows;
  int tfields;
  char fitsval[FLEN_VALUE];
  long pcount;
  
  ret = fits_read_btblhdr ( fp, 100, &nrows, &tfields, NULL, NULL, NULL, fitsval, &pcount, &status );
  fits::check ( status );
  
  if ( offset + nread > nrows ) {
    HARP_THROW( "binary read range is beyond end of table" );
  } 
  
  // check that column number is in range
  
  if ( ( col >= tfields ) || ( col < 0 ) ) {
    ostringstream o;
    o << "cannot read (zero-based) column " << col << " from binary table with " << tfields << " columns";
    HARP_THROW( o.str().c_str() );
  }

  // temporary char buffer

  char ** charray;
  charray = (char**) malloc ( nread * sizeof( char* ) );
  if ( ! charray ) {
    HARP_THROW( "cannot allocate temp char array" );
  }
  for ( long i = 0; i < nread; ++i ) {
    charray[i] = (char*) malloc ( 50 * sizeof (char) );
    if ( ! charray[i] ) {
      HARP_THROW( "cannot allocate temp char array member" );
    }
  }

  // read data in one shot

  int anynul;

  char empty[5];
  strcpy ( empty, "" );

  ret = fits_read_col_str ( fp, col + 1, offset + 1, 1, nread, empty, charray, &anynul, &status );
  fits::check ( status );

  // copy into output vector

  data.resize ( nread );
  for ( long i = 0; i < nread; ++i ) {
    data[i] = charray[i];
    free ( charray[i] );
  }
  free ( charray );

  return;
}


void harp::fits::bin_write ( fitsfile * fp, size_t firstrow, size_t lastrow, std::vector < int > & columns, std::vector < matrix_local > & data ) {
  
  HARP_THROW( "binary FITS table writing not yet implemented" );
  
  return;
}


void harp::fits::test ( string const & datadir ) {

  int np;
  int myp;

  MPI_Comm_size ( MPI_COMM_WORLD, &np );
  MPI_Comm_rank ( MPI_COMM_WORLD, &myp );
  
  if ( myp == 0 ) {
    cerr << "Testing FITS operations..." << endl;
  
    string imgfile = datadir + "/" + "fits_test_img.fits.out";
    
    size_t rows = 100;
    size_t cols = 500;
    size_t nelems = rows * cols;
    
    matrix_local data ( nelems, 1 );
    
    for ( size_t i = 0; i < cols; ++i ) {
      for ( size_t j = 0; j < rows; ++j ) {
        data.Set ( i*rows+j, 0, (double)(i * rows + j) );
      }
    }

    fitsfile * fp;
    
    fits::open_readwrite ( fp, imgfile );
    
    fits::img_append ( fp, rows, cols );
    
    fits::img_write ( fp, data );

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

    matrix_local checkdata ( nelems, 1 );

    fits::img_read ( fp, checkdata );
    
    for ( size_t i = 0; i < cols; ++i ) {
      for ( size_t j = 0; j < rows; ++j ) {
        if ( checkdata.Get( i*rows+j, 0 ) != data.Get( i*rows+j, 0 ) ) {
          cerr << "  (FAILED): img element (" << i << ", " << j << ") has wrong value (" << checkdata.Get( i*rows+j, 0 ) << " != " << data.Get( i*rows+j, 0 ) << ")" << endl;
          exit(1);
        }
      }
    }

    fits::close ( fp );

    cerr << "  (PASSED)" << endl;
  }
  
  return;
}


