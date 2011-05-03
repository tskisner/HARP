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

  void img_dims ( fitsfile * fp, size_t & rows, size_t & cols );

} }

#endif

