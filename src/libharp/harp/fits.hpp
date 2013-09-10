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
  };

  template < >
  struct ftype < double > {
    static int datatype () { return TDOUBLE; }
    static int bitpix () { return -64; }
  };

  template < >
  struct ftype < float > {
    static int datatype () { return TFLOAT; }
    static int bitpix () { return -32; }
  };

  template < >
  struct ftype < long long int > {
    static int datatype () { return TLONGLONG; }
    static int bitpix () { return 64; }
  };

  template < >
  struct ftype < int > {
    static int datatype () { return TINT32BIT; }
    static int bitpix () { return 32; }
  };

  template < >
  struct ftype < short int > {
    static int datatype () { return TSHORT; }
    static int bitpix () { return 16; }
  };

  template < >
  struct ftype < unsigned char > {
    static int datatype () { return TBYTE; }
    static int bitpix () { return 8; }
  };

  // seek operations

  int seek ( fitsfile * fp, std::string const & extname );
  
  int img_seek ( fitsfile * fp, std::string const & keyname, std::string const & keyval );
  
  void img_seek ( fitsfile * fp, int hdu );

  int bin_seek ( fitsfile * fp, std::string const & keyname, std::string const & keyval );
  
  void bin_seek ( fitsfile * fp, int hdu );


  // image operations


  void img_dims ( fitsfile * fp, size_t & rows, size_t & cols );


  template < typename T >
  void img_append ( fitsfile * fp, size_t rows, size_t cols ) {

    int ret;
    int status = 0;
    
    long naxes[2];
    naxes[0] = cols;
    naxes[1] = rows;
    
    ret = fits_create_img ( fp, ftype< T >::bitpix(), 2, naxes, &status );
    fits::check ( status );

    return;
  }


  template < class V >
  void img_write ( fitsfile * fp, boost::numeric::ublas::vector_expression < V > & data ) {

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

    ret = fits_write_pix ( fp, fitstype, fpixel, npix, &( data()[0] ), &status );
    fits::check ( status );

    return;
  }


  template < class M >
  void img_write ( fitsfile * fp, boost::numeric::ublas::matrix_expression < M > const & data ) {

    typedef M matrix_type;
    typedef typename matrix_type::value_type value_type;

    check_column_major ( data() );

    boost::numeric::ublas::vector < value_type > buffer ( data().size1() * data().size2() );

    // transpose data when copying

    for ( size_t col = 0; col < data().size1(); ++col ) {
      for ( size_t row = 0; row < data().size2(); ++row ) {
        buffer[ row * data().size1() + col ] = data()( col, row );
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

    ret = fits_read_pix ( fp, fitstype, fpixel, nelem, NULL, &( data()[0] ), &anynul, &status );
    fits::check ( status );

    return;
  }


  template < class M >
  void img_read ( fitsfile * fp, boost::numeric::ublas::matrix_expression < M > & data ) {

    typedef M matrix_type;
    typedef typename matrix_type::value_type value_type;

    check_column_major ( data() );

    size_t rows;
    size_t cols;
    img_dims ( fp, rows, cols );

    data().resize ( cols, rows );

    boost::numeric::ublas::vector < value_type > buffer ( cols * rows );

    img_read ( fp, buffer );

    // transpose data when copying

    for ( size_t col = 0; col < data().size1(); ++col ) {
      for ( size_t row = 0; row < data().size2(); ++row ) {
        data()( col, row ) = buffer[ row * data().size1() + col ];
      }
    }

    return;
  }


  // binary table operations
  
  std::vector < int > bin_columns ( fitsfile * fp, std::vector < std::string > & names );

  void bin_read_strings ( fitsfile * fp, size_t firstrow, size_t lastrow, int col, std::vector < std::string > & data );
  
  
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

  
  template < class V >
  void bin_write ( fitsfile * fp, size_t firstrow, size_t lastrow, std::vector < int > & columns, std::vector < boost::numeric::ublas::vector_expression < V > > & data ) {

    typedef V vector_type;
    typedef typename vector_type::value_type value_type;

    HARP_THROW( "binary FITS table writing not yet implemented" );

    return;
  }


  // sanity tests...
  
  void test ( std::string const & datadir );

} }

#endif

