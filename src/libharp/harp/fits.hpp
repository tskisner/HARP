// @COPYRIGHT@

#ifndef HARP_FITS_HPP
#define HARP_FITS_HPP


extern "C" {
  #include <fitsio.h>
}


namespace harp { namespace fits {

  void check ( int status );
  
  void open_read ( fitsfile * & fp, std::string const & path );
  
  void close ( fitsfile * fp );
  
  void img_seek ( fitsfile * fp, std::string const & extname );
  
  void img_seek ( fitsfile * fp, int hdu );

  void img_dims ( fitsfile * fp, size_t & rows, size_t & cols );
  
  void img_read ( fitsfile * fp, size_t frow, size_t fcol, dense_mat_view & data );
  
  void img_read_row ( fitsfile * fp, size_t row, data_vec & data );

} }

#endif

