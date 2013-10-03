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


int harp::fits::seek ( fitsfile * fp, string const & extname ) {
  int hdu;
  
  int ret;
  int status = 0;
  
  char extcopy[FLEN_VALUE];
  strncpy ( extcopy, extname.c_str(), FLEN_VALUE );
  
  ret = fits_movnam_hdu ( fp, ANY_HDU, extcopy, 0, &status );
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


void harp::fits::bin_create ( fitsfile * fp, std::string extname, size_t nrows, std::vector < std::string > colnames, std::vector < std::string > coltypes, std::vector < std::string > colunits ) {

  size_t ncols = colnames.size();

  if ( ( ncols != coltypes.size() ) || ( ncols != colunits.size() ) ) {
    HARP_THROW( "vectors of column names, types, and units do not have equal lengths" );
  }

  char ** ttype;
  char ** tform;
  char ** tunit;

  ttype = (char**) malloc ( ncols * sizeof(char*) );
  tform = (char**) malloc ( ncols * sizeof(char*) );
  tunit = (char**) malloc ( ncols * sizeof(char*) );

  if ( ! ( ttype && tform && tunit ) ) {
    ostringstream o;
    o << "cannot allocate column info for binary table";
    HARP_THROW( o.str().c_str() );
  }

  for ( size_t c = 0; c < ncols; ++c ) {
    ttype[c] = (char*) malloc ( FLEN_VALUE * sizeof(char) );
    tform[c] = (char*) malloc ( FLEN_VALUE * sizeof(char) );
    tunit[c] = (char*) malloc ( FLEN_VALUE * sizeof(char) );
    if ( ! ( ttype[c] && tform[c] && tunit[c] ) ) {
      ostringstream o;
      o << "cannot allocate column info for binary table";
      HARP_THROW( o.str().c_str() );
    }
  }
  for ( size_t i = 0; i < ncols; ++i ) {
    strncpy ( ttype[i], colnames[i].c_str(), FLEN_VALUE );
    strncpy ( tform[i], coltypes[i].c_str(), FLEN_VALUE );
    strncpy ( tunit[i], colunits[i].c_str(), FLEN_VALUE );
  }

  int status = 0;

  char cextname[ FLEN_VALUE ];
  strncpy ( cextname, extname.c_str(), FLEN_VALUE );

  int ret = fits_create_tbl ( fp, BINARY_TBL, (long)nrows, (int)ncols, ttype, tform, tunit, cextname, &status );
  fits::check ( status );

  for ( size_t c = 0; c < ncols; ++c ) {
    free ( ttype[c] );
    free ( tform[c] );
    free ( tunit[c] );
  }
  free ( ttype );
  free ( tform );
  free ( tunit );

  return;
}


void harp::fits::bin_info ( fitsfile * fp, size_t & nrows, vector < string > & colnames ) {

  int ret;
  int status = 0;

  // get table dimensions
  
  long lnrows;
  int tfields;
  char fitsval[FLEN_VALUE];
  long pcount;
  
  ret = fits_read_btblhdr ( fp, 100, &lnrows, &tfields, NULL, NULL, NULL, fitsval, &pcount, &status );
  fits::check ( status );
  
  nrows = lnrows;

  colnames.resize ( tfields );

  char templt[FLEN_VALUE];
  strcpy ( templt, "*" );

  int colnum;
  char colnam[FLEN_VALUE];

  for ( int i = 0; i < tfields; ++i ) {
    ret = fits_get_colname ( fp, CASEINSEN, templt, colnam, &colnum, &status );
    if ( colnum != i + 1 ) {
      HARP_THROW( "error in binary table column name lookup" );
    }
    colnames[i] = colnam;
  }
  status = 0;

  return;
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


void harp::fits::bin_read_column_strings ( fitsfile * fp, size_t firstrow, size_t lastrow, int col, vector < string > & data ) {
  
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
    charray[i] = (char*) malloc ( FLEN_VALUE * sizeof (char) );
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


void harp::fits::bin_write_column_strings ( fitsfile * fp, size_t firstrow, size_t lastrow, int col, vector < string > & data ) {
  
  int ret;
  int status = 0;
  long offset = (long)firstrow;
  long nwrite = (long)lastrow - offset + 1;
  
  // get table dimensions
  
  long nrows;
  int tfields;
  char fitsval[FLEN_VALUE];
  long pcount;
  
  ret = fits_read_btblhdr ( fp, 100, &nrows, &tfields, NULL, NULL, NULL, fitsval, &pcount, &status );
  fits::check ( status );
  
  if ( offset + nwrite > nrows ) {
    HARP_THROW( "binary write range is beyond end of table" );
  } 
  
  // check that column number is in range
  
  if ( ( col >= tfields ) || ( col < 0 ) ) {
    ostringstream o;
    o << "cannot write (zero-based) column " << col << " from binary table with " << tfields << " columns";
    HARP_THROW( o.str().c_str() );
  }

  // temporary char buffer

  char ** charray;
  charray = (char**) malloc ( nwrite * sizeof( char* ) );
  if ( ! charray ) {
    HARP_THROW( "cannot allocate temp char array" );
  }
  for ( long i = 0; i < nwrite; ++i ) {
    charray[i] = (char*) malloc ( FLEN_VALUE * sizeof (char) );
    if ( ! charray[i] ) {
      HARP_THROW( "cannot allocate temp char array member" );
    }
  }

  // copy into output vector

  for ( long i = 0; i < nwrite; ++i ) {
    strncpy ( charray[i], data[i].c_str(), FLEN_VALUE );
  }

  // write data in one shot

  ret = fits_write_col_str ( fp, col + 1, offset + 1, 1, nwrite, charray, &status );
  fits::check ( status );

  // free memory

  for ( long i = 0; i < nwrite; ++i ) {
    free ( charray[i] );
  }
  free ( charray );

  return;
}




void harp::fits::test ( string const & datadir ) {
  
  cerr << "Testing FITS operations..." << endl;

  string imgfile = datadir + "/" + "fits_test_img.fits.out";
  
  size_t rows = 100;
  size_t cols = 500;
  size_t nelems = rows * cols;
  
  matrix_double data ( cols, rows );
  
  for ( size_t i = 0; i < cols; ++i ) {
    for ( size_t j = 0; j < rows; ++j ) {
      data( i, j ) = (double)(i * rows + j);
    }
  }

  fitsfile * fp;
  
  fits::open_readwrite ( fp, imgfile );
  
  fits::img_append < double > ( fp, rows, cols );
  
  fits::img_write ( fp, data );

  size_t tabrows = 10;

  vector < string > colnames ( 7 );
  vector < string > coltypes ( 7 );
  vector < string > colunits ( 7 );

  colnames[0] = "DOUBLE";
  coltypes[0] = "1" + fits::ftype < double > :: coltype();
  colunits[0] = "None";

  colnames[1] = "FLOAT";
  coltypes[1] = "1" + fits::ftype < float > :: coltype();
  colunits[1] = "None";

  colnames[2] = "LONGLONG";
  coltypes[2] = "1" + fits::ftype < long long int > :: coltype();
  colunits[2] = "None";

  colnames[3] = "INT";
  coltypes[3] = "1" + fits::ftype < int > :: coltype();
  colunits[3] = "None";

  colnames[4] = "SHORT";
  coltypes[4] = "1" + fits::ftype < short int > :: coltype();
  colunits[4] = "None";

  colnames[5] = "CHAR";
  coltypes[5] = "1" + fits::ftype < unsigned char > :: coltype();
  colunits[5] = "None";

  colnames[6] = "STRING";
  coltypes[6] = "6A";
  colunits[6] = "None";


  fits::bin_create ( fp, string("TEST"), tabrows, colnames, coltypes, colunits );

  boost::numeric::ublas::vector < double > data_double ( tabrows );
  boost::numeric::ublas::vector < float > data_float ( tabrows );
  boost::numeric::ublas::vector < long long > data_longlong ( tabrows );
  boost::numeric::ublas::vector < int > data_int ( tabrows );
  boost::numeric::ublas::vector < short int > data_short ( tabrows );
  boost::numeric::ublas::vector < unsigned char > data_char ( tabrows );
  vector < string > names ( rows );

  for ( size_t i = 0; i < tabrows; ++i ) {
    names[i] = "Blah";
    data_double[i] = (double)i;
    data_float[i] = (float)i;
    data_longlong[i] = (long long)i;
    data_int[i] = (int)i;
    data_short[i] = (short)i;
    data_char[i] = (unsigned char)i;
  }

  fits::bin_write_column ( fp, 0, tabrows - 1, 0, data_double );
  fits::bin_write_column ( fp, 0, tabrows - 1, 1, data_float );
  fits::bin_write_column ( fp, 0, tabrows - 1, 2, data_longlong );
  fits::bin_write_column ( fp, 0, tabrows - 1, 3, data_int );
  fits::bin_write_column ( fp, 0, tabrows - 1, 4, data_short );
  fits::bin_write_column ( fp, 0, tabrows - 1, 5, data_char ); 

  fits::bin_write_column_strings ( fp, 0, tabrows - 1, 6, names );

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

  matrix_double checkdata ( cols, rows );

  fits::img_read ( fp, checkdata );
  
  for ( size_t i = 0; i < cols; ++i ) {
    for ( size_t j = 0; j < rows; ++j ) {
      if ( checkdata( i, j ) != data( i, j ) ) {
        cerr << "  (FAILED): img element (" << i << ", " << j << ") has wrong value (" << checkdata( i, j ) << " != " << data( i, j ) << ")" << endl;
        exit(1);
      }
    }
  }

  fits::close ( fp );



  cerr << "  (PASSED)" << endl;
  
  return;
}


