// @COPYRIGHT@

#ifndef HARP_FITS_HPP
#define HARP_FITS_HPP


extern "C" {
  #include <fitsio.h>
}


namespace harp { namespace fits {

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
  
  void img_append ( fitsfile * fp, size_t rows, size_t cols );
  
  void img_write ( fitsfile * fp, matrix_local & data );
  
  int img_seek ( fitsfile * fp, std::string const & extname );
  
  int img_seek ( fitsfile * fp, std::string const & keyname, std::string const & keyval );
  
  void img_seek ( fitsfile * fp, int hdu );

  void img_dims ( fitsfile * fp, size_t & rows, size_t & cols );
  
  void img_read ( fitsfile * fp, matrix_local & data );
  
  std::vector < int > bin_columns ( fitsfile * fp, std::vector < std::string > & names );
  
  int bin_seek ( fitsfile * fp, std::string const & keyname, std::string const & keyval );
  
  void bin_seek ( fitsfile * fp, int hdu );
  
  void bin_read ( fitsfile * fp, size_t firstrow, size_t lastrow, std::vector < int > & columns, std::vector < matrix_local > & data );

  void bin_read_strings ( fitsfile * fp, size_t firstrow, size_t lastrow, int col, std::vector < std::string > & data );
  
  void bin_write ( fitsfile * fp, size_t firstrow, size_t lastrow, std::vector < int > & columns, std::vector < matrix_local > & data );
  
  void test ( std::string const & datadir );

} }

#endif

