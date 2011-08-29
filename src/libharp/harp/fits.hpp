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
  
  void close ( fitsfile * fp );
  
  int nhdus ( fitsfile * fp );
  
  std::string key_string ( fitsfile * fp, std::string keyname );
  
  long key_long ( fitsfile * fp, std::string keyname );
  
  double key_double ( fitsfile * fp, std::string keyname );
  
  void img_append ( fitsfile * fp, size_t rows, size_t cols );
  
  void img_write ( fitsfile * fp, size_t frow, size_t fcol, mat_denserow & data );
  
  void img_write ( fitsfile * fp, vec_dense & data );
  
  int img_seek ( fitsfile * fp, std::string const & extname );
  
  int img_seek ( fitsfile * fp, std::string const & keyname, std::string const & keyval );
  
  void img_seek ( fitsfile * fp, int hdu );

  void img_dims ( fitsfile * fp, size_t & rows, size_t & cols );
  
  void img_read ( fitsfile * fp, size_t frow, size_t fcol, mat_denserow & data );
  
  void img_read ( fitsfile * fp, vec_dense & data );
  
  void img_read_row ( fitsfile * fp, size_t row, vec_dense & data );
  
  void img_write_row ( fitsfile * fp, size_t row, vec_dense & data );
  
  void img_read_row_int ( fitsfile * fp, size_t row, vec_int & data );
  
  std::vector < int > bin_columns ( fitsfile * fp, std::vector < std::string > & names );
  
  int bin_seek ( fitsfile * fp, std::string const & keyname, std::string const & keyval );
  
  void bin_seek ( fitsfile * fp, int hdu );
  
  void bin_read ( fitsfile * fp, size_t firstrow, size_t lastrow, std::vector < int > & columns, std::vector < vec_dense > & data );
  
  void bin_write ( fitsfile * fp, size_t firstrow, size_t lastrow, std::vector < int > & columns, std::vector < vec_dense > & data );
  
  void test ( std::string const & datadir );

} }

#endif

