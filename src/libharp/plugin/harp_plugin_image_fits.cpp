// @COPYRIGHT@

#include <harp_data_internal.hpp>

#include <harp/plugin.hpp>
//#ifdef HAVE_BOOST_MPI_HPP
//#include <harp/plugin_mpi.hpp>
//#endif

#include <harp/static_plugins.hpp>

using namespace std;
using namespace harp;

static const char * image_fits_key_path = "path";
static const char * image_fits_key_signal = "signal";
static const char * image_fits_key_noise = "noise";
static const char * image_fits_key_rows = "rows";
static const char * image_fits_key_cols = "cols";


harp::image_fits::image_fits ( boost::property_tree::ptree const & props ) : image ( "fits", props ) {

  sighdu_ = props.get ( image_fits_key_signal, 1 );

  nsehdu_ = props.get ( image_fits_key_noise, 2 );

  path_ = props.get ( image_fits_key_path, "" );

  if ( path_ == "" ) {

    // no path specified- must specify rows / cols

    rows_ = props.get < size_t > ( image_fits_key_rows );

    cols_ = props.get < size_t > ( image_fits_key_cols );

  } else {
    
    // read rows / cols from the FITS header

    fitsfile *fp;
    
    fits::open_read ( fp, path_ );

    fits::img_seek ( fp, sighdu_ );
    
    fits::img_dims ( fp, cols_, rows_ );

    size_t rowcheck;
    size_t colcheck;

    fits::img_seek ( fp, nsehdu_ );
    
    fits::img_dims ( fp, colcheck, rowcheck );

    if ( ( rowcheck != rows_ ) || ( colcheck != cols_ ) ) {
      HARP_THROW( "noise variance dimensions do not match data dimensions" );
    }

  }  
  
}


harp::image_fits::~image_fits ( ) {

}


boost::property_tree::ptree harp::image_fits::metadata ( ) const {

  return boost::property_tree::ptree();
}


void harp::image_fits::values ( vector_double & data ) const {

  fitsfile *fp;

  fits::open_read ( fp, path_ );

  fits::img_seek ( fp, sighdu_ );

  vector_double buffer;

  fits::img_read ( fp, buffer );

  data.resize ( buffer.size() );
  fits::img_transpose ( cols_, rows_, buffer, data );

  fits::close ( fp );

  return;
}


void harp::image_fits::inv_variance ( vector_double & invvar ) const {

  fitsfile *fp;

  fits::open_read ( fp, path_ );

  fits::img_seek ( fp, nsehdu_ );

  vector_double buffer;
    
  fits::img_read ( fp, buffer );

  invvar.resize ( buffer.size() );
  fits::img_transpose ( cols_, rows_, buffer, invvar );

  fits::close ( fp );

  return;
}


void harp::image_fits::write ( std::string const & path, vector_double & data, vector_double & invvar ) {

  if ( data.size() != rows_ * cols_ ) {
    HARP_THROW( "data size does not match image dimensions" );
  }

  if ( invvar.size() != rows_ * cols_ ) {
    HARP_THROW( "inverse variance size does not match image dimensions" );
  }

  fitsfile * fp;
    
  fits::create ( fp, path );
  
  fits::img_append < double > ( fp, cols_, rows_ );
  
  vector_double buffer ( data.size() );

  fits::img_transpose ( rows_, cols_, data, buffer );

  fits::img_write ( fp, buffer );

  fits::img_append < double > ( fp, cols_, rows_ );
  
  fits::img_transpose ( rows_, cols_, invvar, buffer );

  fits::img_write ( fp, buffer );

  fits::close ( fp );

  return;
}


void harp::image_fits::write ( std::string const & path, matrix_double & data, matrix_double & invvar ) {

  size_t nelem = rows_ * cols_;

  if ( ( rows_ != data.size2() ) || ( cols_ != data.size1() ) ) {
    HARP_THROW( "data size does not match image dimensions" );
  }

  if ( ( rows_ != invvar.size2() ) || ( cols_ != invvar.size1() ) ) {
    HARP_THROW( "inverse variance size does not match image dimensions" );
  }

  vector_double tempdata ( nelem );
  vector_double tempvar ( nelem );

  for ( size_t i = 0; i < cols_; ++i ) {
    for ( size_t j = 0; j < rows_; ++j ) {
      tempdata[ i * rows_ + j ] = data( i, j );
      tempvar[ i * rows_ + j ] = invvar( i, j );
    }
  }

  write ( path, tempdata, tempvar );

  return;
}


BOOST_CLASS_EXPORT(harp::image_fits)


image * harp::image_fits_create ( boost::property_tree::ptree const & props ) {
  return new image_fits ( props );
}



