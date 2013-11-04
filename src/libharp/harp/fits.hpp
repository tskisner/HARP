// @COPYRIGHT@

#ifndef HARP_FITS_HPP
#define HARP_FITS_HPP


extern "C" {
  #include <fitsio.h>
}


namespace harp { namespace fits {

  // general manipulation of files and keywords

  void check ( int status );
  
  void open_read ( fitsfile * & fp, std::string const & path );
  
  void open_readwrite ( fitsfile * & fp, std::string const & path );

  void create ( fitsfile * & fp, std::string const & path );
  
  void close ( fitsfile * fp );
  
  int nhdus ( fitsfile * fp );

  void read_key ( fitsfile * fp, std::string const & keyname, std::string & keyval );
  void read_key ( fitsfile * fp, std::string const & keyname, long & keyval );
  void read_key ( fitsfile * fp, std::string const & keyname, double & keyval );

  void write_key ( fitsfile * fp, std::string const & keyname, std::string const & keyval, std::string const & keycom );
  void write_key ( fitsfile * fp, std::string const & keyname, long const & keyval, std::string const & keycom );
  void write_key ( fitsfile * fp, std::string const & keyname, double const & keyval, std::string const & keycom );
  
  // FITS data type map

  template < typename T >
  struct ftype {
    static int datatype () { return 0; }
    static int bitpix () { return 0; }
    static std::string coltype () { return std::string(""); }
  };

  template < >
  struct ftype < double > {
    static int datatype () { return TDOUBLE; }
    static int bitpix () { return -64; }
    static std::string coltype () { return std::string("D"); }
  };

  template < >
  struct ftype < float > {
    static int datatype () { return TFLOAT; }
    static int bitpix () { return -32; }
    static std::string coltype () { return std::string("E"); }
  };

  template < >
  struct ftype < long long int > {
    static int datatype () { return TLONGLONG; }
    static int bitpix () { return 64; }
    static std::string coltype () { return std::string("K"); }
  };

  template < >
  struct ftype < int > {
    static int datatype () { return TINT32BIT; }
    static int bitpix () { return 32; }
    static std::string coltype () { return std::string("J"); }
  };

  template < >
  struct ftype < short int > {
    static int datatype () { return TSHORT; }
    static int bitpix () { return 16; }
    static std::string coltype () { return std::string("I"); }
  };

  template < >
  struct ftype < unsigned char > {
    static int datatype () { return TBYTE; }
    static int bitpix () { return 8; }
    static std::string coltype () { return std::string("B"); }
  };

  // seek operations

  int seek ( fitsfile * fp, std::string const & extname );
  
  int img_seek ( fitsfile * fp, std::string const & keyname, std::string const & keyval );
  
  void img_seek ( fitsfile * fp, int hdu );

  int bin_seek ( fitsfile * fp, std::string const & keyname, std::string const & keyval );
  
  void bin_seek ( fitsfile * fp, int hdu );


  // image operations

  bool img_colmajor ( fitsfile * fp );


  void img_dims ( fitsfile * fp, size_t & rows, size_t & cols );


  template < typename T >
  void img_append ( fitsfile * fp, size_t rows, size_t cols, bool col_major = false ) {

    int ret;
    int status = 0;
    
    long naxes[2];
    if ( col_major ) {
      naxes[0] = rows;
      naxes[1] = cols;
    } else {
      naxes[0] = cols;
      naxes[1] = rows;
    }
    
    ret = fits_create_img ( fp, ftype< T >::bitpix(), 2, naxes, &status );
    fits::check ( status );

    if ( col_major ) {
      write_key ( fp, "HRPCOLMJ", "TRUE", "Whether the data is column-major" );
    }

    return;
  }


  template < class V >
  void img_write ( fitsfile * fp, boost::numeric::ublas::vector_expression < V > const & data ) {

    typedef V vector_type;
    typedef typename vector_type::value_type value_type;

    int ret;
    int status = 0;

    long fpixel[2];
    
    fpixel[0] = 1;
    fpixel[1] = 1;

    long npix = (long) data().size();

    size_t rows;
    size_t cols;
    img_dims ( fp, rows, cols );

    if ( rows * cols != data().size() ) {
      HARP_THROW( "data vector size not consistent with image dimensions" );
    }

    int fitstype = ftype < value_type > :: datatype();

    bool colmajor = img_colmajor ( fp );

    // copy data to a buffer to work around stupid CFITSIO non-const API

    value_type * buffer = (value_type*) malloc ( npix * sizeof( value_type ) );

    if ( ! buffer ) {
      HARP_THROW( "cannot allocate img buffer" );
    }

    if ( colmajor ) {
      // we can just copy
      for ( long i = 0; i < npix; ++i ) {
        buffer[i] = data()[i];
      }
    } else {
      // we need to remap the memory

      for ( size_t i = 0; i < rows; ++i ) {
        for ( size_t j = 0; j < cols; ++j ) {
          buffer[ i * cols + j ] = data()[ j * rows + i ];
        }
      }

    }

    ret = fits_write_pix ( fp, fitstype, fpixel, npix, buffer, &status );
    fits::check ( status );

    free ( buffer );

    return;
  }


  template < class M >
  void img_write ( fitsfile * fp, boost::numeric::ublas::matrix_expression < M > const & data ) {

    typedef M matrix_type;
    typedef typename matrix_type::value_type value_type;

    boost::numeric::ublas::vector < value_type > buffer ( data().size1() * data().size2() );

    // copy data

    size_t elem = 0;

    for ( size_t col = 0; col < data().size1(); ++col ) {
      for ( size_t row = 0; row < data().size2(); ++row ) {
        buffer[ elem ] = data()( col, row );
        ++elem;
      }
    }

    img_write ( fp, buffer );

    return;
  }
  

  template < class V >
  void img_read ( fitsfile * fp, boost::numeric::ublas::vector_expression < V > & data ) {

    typedef V vector_type;
    typedef typename vector_type::value_type value_type;

    int ret;
    int status = 0;
    int anynul;

    size_t rows;
    size_t cols;
    img_dims ( fp, rows, cols );

    long nelem = (long)( rows * cols );

    data().resize ( nelem );

    long fpixel[2];

    fpixel[0] = 1;
    fpixel[1] = 1;

    int fitstype = ftype < value_type > :: datatype();

    bool colmajor = img_colmajor ( fp );

    if ( colmajor ) {
      // we can just read it directly

      ret = fits_read_pix ( fp, fitstype, fpixel, nelem, NULL, &( data()[0] ), &anynul, &status );
      fits::check ( status );

    } else {
      // we need to remap the memory

      value_type * buffer = (value_type*) malloc ( nelem * sizeof( value_type ) );

      if ( ! buffer ) {
        HARP_THROW( "cannot allocate img buffer" );
      }

      ret = fits_read_pix ( fp, fitstype, fpixel, nelem, NULL, buffer, &anynul, &status );
      fits::check ( status );

      for ( size_t i = 0; i < rows; ++i ) {
        for ( size_t j = 0; j < cols; ++j ) {
          data()[ j * rows + i ] = buffer[ i * cols + j ];
        }
      }

      free ( buffer );

    }

    return;
  }


  template < class M >
  void img_read ( fitsfile * fp, boost::numeric::ublas::matrix_expression < M > & data ) {

    typedef M matrix_type;
    typedef typename matrix_type::value_type value_type;

    size_t rows;
    size_t cols;
    img_dims ( fp, rows, cols );

    data().resize ( cols, rows );

    boost::numeric::ublas::vector < value_type > buffer ( cols * rows );

    img_read ( fp, buffer );

    // copy data

    size_t elem = 0;
    for ( size_t col = 0; col < data().size1(); ++col ) {
      for ( size_t row = 0; row < data().size2(); ++row ) {
        data()( col, row ) = buffer[ elem ];
        ++elem;
      }
    }

    return;
  }


  // binary table operations


  void bin_create ( fitsfile * fp, std::string extname, size_t nrows, std::vector < std::string > colnames, std::vector < std::string > coltypes, std::vector < std::string > colunits );


  void bin_info ( fitsfile * fp, size_t & nrows, std::vector < std::string > & colnames );

  
  std::vector < int > bin_columns ( fitsfile * fp, std::vector < std::string > & names );


  template < class V >
  void bin_read_column ( fitsfile * fp, size_t firstrow, size_t lastrow, int column, boost::numeric::ublas::vector_expression < V > & data ) {

    typedef V vector_type;
    typedef typename vector_type::value_type value_type;

    int fitstype = ftype < value_type > :: datatype();

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
    
    if ( ( column >= tfields ) || ( column < 0 ) ) {
      std::ostringstream o;
      o << "cannot read (zero-based) column " << column << " from binary table with " << tfields << " columns";
      HARP_THROW( o.str().c_str() );
    }

    data().resize ( nread );

    int anynul;

    ret = fits_read_col ( fp, fitstype, column + 1, offset + 1, 1, nread, 0, &(data()[0]), &anynul, &status );
    fits::check ( status );

    return;
  }


  template < class V >
  void bin_write_column ( fitsfile * fp, size_t firstrow, size_t lastrow, int column, boost::numeric::ublas::vector_expression < V > const & data ) {

    typedef V vector_type;
    typedef typename vector_type::value_type value_type;

    int fitstype = ftype < value_type > :: datatype();

    int ret;
    int status = 0;
    long offset = (long)firstrow;
    long nwrite = (long)lastrow - offset + 1;

    if ( data().size() < nwrite ) {
      HARP_THROW( "data vector is shorter than number of rows to write" );
    }

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
    
    if ( ( column >= tfields ) || ( column < 0 ) ) {
      std::ostringstream o;
      o << "cannot write (zero-based) column " << column << " from binary table with " << tfields << " columns";
      HARP_THROW( o.str().c_str() );
    }

    // copy data to a buffer to work around stupid CFITSIO non-const API

    value_type * buffer = (value_type*) malloc ( nwrite * sizeof( value_type ) );

    if ( ! buffer ) {
      HARP_THROW( "cannot allocate bintable buffer" );
    }
    for ( long i = 0; i < nwrite; ++i ) {
      buffer[i] = data()[i];
    }

    ret = fits_write_col ( fp, fitstype, column + 1, offset + 1, 1, nwrite, buffer, &status );
    fits::check ( status );

    free ( buffer );

    return;
  }


  void bin_read_column_strings ( fitsfile * fp, size_t firstrow, size_t lastrow, int col, std::vector < std::string > & data );


  void bin_write_column_strings ( fitsfile * fp, size_t firstrow, size_t lastrow, int col, std::vector < std::string > & data );
  
  
  template < class V >
  void bin_read ( fitsfile * fp, size_t firstrow, size_t lastrow, std::vector < int > & columns, std::vector < boost::numeric::ublas::vector_expression < V > > & data ) {

    typedef V vector_type;
    typedef typename vector_type::value_type value_type;

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
    
    std::vector < int > :: iterator it;
    
    int cur = 0;
    for ( it = columns.begin(); it != columns.end(); ++it ) {
      if ( ( (*it) >= tfields ) || ( (*it) < 0 ) ) {
        std::ostringstream o;
        o << "cannot read (zero-based) column " << (*it) << " from binary table with " << tfields << " columns";
        HARP_THROW( o.str().c_str() );
      }
      data[ cur ]().resize ( nread );
      ++cur;
    }
      
    // read data in a buffered way

    int fitstype = ftype < value_type > :: datatype();
    
    int anynul;
    long n = optimal;
    long dataoffset = 0;

    while ( n == optimal ) {
      if ( dataoffset + optimal > nread ) {
        n = nrows - dataoffset;
      }
      
      cur = 0;
      for ( it = columns.begin(); it != columns.end(); ++it ) {
        ret = fits_read_col ( fp, fitstype, (*it) + 1, offset + 1, 1, n, 0, &((data[cur])[dataoffset]), &anynul, &status );
        fits::check ( status );
        //cerr << "fits::bin_read col " << (*it)+1 << " rows " << offset+1 << "-" << offset+1+n << " = " << (data[cur])[dataoffset] << ", " << (data[cur])[dataoffset+1] << ", " << (data[cur])[dataoffset+2] << endl;
        ++cur;
      }
      
      offset += optimal;
      dataoffset += optimal;
    }

    return;
  }



  // sanity tests...
  
  void test ( std::string const & datadir );

} }

#endif

