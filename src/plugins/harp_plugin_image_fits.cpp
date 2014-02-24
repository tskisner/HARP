// @COPYRIGHT@

#ifdef USE_MPI
#  include <harp_mpi.hpp>
#else
#  include <harp.hpp>
#endif

#include <harp/plugin.hpp>
#ifdef USE_MPI
#  include <harp/mpi_plugin.hpp>
#endif


#include <internal_image_plugins.hpp>

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
    
    fits::img_dims ( fp, rows_, cols_ );

    size_t rowcheck;
    size_t colcheck;

    fits::img_seek ( fp, nsehdu_ );
    
    fits::img_dims ( fp, rowcheck, colcheck );

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

  fits::img_read ( fp, data, true );

  fits::close ( fp );

  return;
}


void harp::image_fits::inv_variance ( vector_double & invvar ) const {

  fitsfile *fp;

  fits::open_read ( fp, path_ );

  fits::img_seek ( fp, nsehdu_ );

  fits::img_read ( fp, invvar, true );

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
  
  fits::img_append < double > ( fp, rows_, cols_ );
  
  fits::img_write ( fp, data, true );

  fits::img_append < double > ( fp, rows_, cols_ );
  
  fits::img_write ( fp, invvar, true );

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



