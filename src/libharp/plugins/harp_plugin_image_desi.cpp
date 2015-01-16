/*
  High Performance Astrophysical Reconstruction and Processing (HARP)
  (c) 2014-2015, The Regents of the University of California, 
  through Lawrence Berkeley National Laboratory.  See top
  level LICENSE file for details.
*/


#include <harp/static_plugins.hpp>

using namespace std;
using namespace harp;

static const char * image_desi_key_path = "path";
static const char * image_desi_key_camera = "CAMERA";
static const char * image_desi_key_vspecter = "VSPECTER";
static const char * image_desi_key_exptime = "EXPTIME";
static const char * image_desi_key_rdnoise = "RDNOISE";
static const char * image_desi_key_flavor = "FLAVOR";



harp::image_desi::image_desi ( boost::property_tree::ptree const & props ) : image ( "desi", props ) {

  sighdu_ = 1;
  nsehdu_ = 2;
  mskhdu_ = 3;

  path_ = props.get ( image_desi_key_path, "" );

  fitsfile *fp;
  
  fits::open_read ( fp, path_ );

  // process first HDU

  fits::img_seek ( fp, sighdu_ );
  fits::img_dims ( fp, rows_, cols_ );

  fits::key_read ( fp, string(image_desi_key_camera), camera_ );

  fits::key_read ( fp, string(image_desi_key_vspecter), vspecter_ );

  fits::key_read ( fp, string(image_desi_key_exptime), exptime_ );

  fits::key_read ( fp, string(image_desi_key_rdnoise), rdnoise_ );

  fits::key_read ( fp, string(image_desi_key_flavor), flavor_ );

  // check second HDU

  size_t rowcheck;
  size_t colcheck;
  string extcheck;

  fits::img_seek ( fp, nsehdu_ );

  fits::key_read ( fp, string("EXTNAME"), extcheck );

  if ( extcheck != "IVAR" ) {
    HARP_THROW( "second HDU is not named \"IVAR\"" );
  }
  
  fits::img_dims ( fp, rowcheck, colcheck );

  if ( ( rowcheck != rows_ ) || ( colcheck != cols_ ) ) {
    HARP_THROW( "inverse noise variance dimensions do not match data dimensions" );
  }

  // check third HDU

  fits::img_seek ( fp, mskhdu_ );

  fits::key_read ( fp, string("EXTNAME"), extcheck );

  if ( extcheck != "MASK" ) {
    HARP_THROW( "third HDU is not named \"MASK\"" );
  }
  
  fits::img_dims ( fp, rowcheck, colcheck );

  if ( ( rowcheck != rows_ ) || ( colcheck != cols_ ) ) {
    HARP_THROW( "mask dimensions do not match data dimensions" );
  }

}


harp::image_desi::~image_desi ( ) {

}


void harp::image_desi::values ( vector_double & data ) const {

  fitsfile *fp;

  fits::open_read ( fp, path_ );

  fits::img_seek ( fp, sighdu_ );

  fits::img_read ( fp, data, true );

  fits::close ( fp );

  return;
}


void harp::image_desi::inv_variance ( vector_double & invvar ) const {

  fitsfile *fp;

  fits::open_read ( fp, path_ );

  fits::img_seek ( fp, nsehdu_ );

  fits::img_read ( fp, invvar, true );

  fits::close ( fp );

  return;
}


void harp::image_desi::mask ( vector_mask & msk ) const {

  fitsfile *fp;

  fits::open_read ( fp, path_ );

  fits::img_seek ( fp, mskhdu_ );

  fits::img_read ( fp, msk, true );

  fits::close ( fp );

  return;
}


void harp::image_desi::write ( std::string const & path, boost::property_tree::ptree const & meta, size_t rows, vector_double & data, vector_double & invvar, vector_mask & msk ) {

  size_t cols = (size_t)( data.size() / rows );

  if ( data.size() != rows * cols ) {
    HARP_THROW( "data size does not match number of rows" );
  }

  if ( invvar.size() != data.size() ) {
    HARP_THROW( "inverse variance size does not match data" );
  }

  if ( msk.size() != data.size() ) {
    HARP_THROW( "mask size does not match data" );
  }

  fitsfile * fp;
    
  fits::create ( fp, path );

  // write signal HDU
  
  fits::img_append < double > ( fp, rows, cols );

  if ( meta.count ( image_desi_key_camera ) ) {
    fits::key_write ( fp, image_desi_key_camera, meta.get < string > ( image_desi_key_camera ), "Spectograph Camera" );
  } else {
    HARP_THROW( "you must specify the camera name" );
  }

  if ( meta.count ( image_desi_key_vspecter ) ) {
    fits::key_write ( fp, image_desi_key_vspecter, meta.get < string > ( image_desi_key_vspecter ), "Specter Version" );
  } else {
    HARP_THROW( "you must specify the specter version" );
  }

  if ( meta.count ( image_desi_key_exptime ) ) {
    fits::key_write ( fp, image_desi_key_exptime, meta.get < float > ( image_desi_key_exptime ), "Exposure time [sec]" );
  } else {
    HARP_THROW( "you must specify the exposure time" );
  }

  if ( meta.count ( image_desi_key_rdnoise ) ) {
    fits::key_write ( fp, image_desi_key_rdnoise, meta.get < float > ( image_desi_key_rdnoise ), "Read noise [electrons]" );
  } else {
    HARP_THROW( "you must specify the read noise" );
  }

  if ( meta.count ( image_desi_key_flavor ) ) {
    fits::key_write ( fp, image_desi_key_flavor, meta.get < string > ( image_desi_key_flavor ), "Exposure type (arc, flat, science)" );
  } else {
    HARP_THROW( "you must specify the exposure flavor" );
  }

  fits::img_write ( fp, data, true );

  // write inverse variance

  fits::img_append < double > ( fp, rows, cols );

  fits::key_write ( fp, "EXTNAME", "IVAR", "" );
  fits::key_write ( fp, image_desi_key_rdnoise, meta.get < float > ( image_desi_key_rdnoise ), "Read noise [electrons]" );

  fits::img_write ( fp, invvar, true );

  // write the mask

  fits::img_append < int > ( fp, rows, cols );

  fits::key_write ( fp, "EXTNAME", "MASK", "" );

  fits::img_write ( fp, msk, true );

  // close file

  fits::close ( fp );

  return;
}


void harp::image_desi::write ( std::string const & path, boost::property_tree::ptree const & meta, matrix_double & data, matrix_double & invvar, matrix_mask & msk ) {

  size_t rows = data.size2();
  size_t cols = data.size1();

  size_t nelem = rows * cols;

  if ( ( rows != invvar.size2() ) || ( cols != invvar.size1() ) ) {
    HARP_THROW( "inverse variance size does not match data dimensions" );
  }

  vector_double tempdata ( nelem );
  vector_double tempvar ( nelem );
  vector_mask tempmask ( nelem );

  for ( size_t i = 0; i < cols; ++i ) {
    for ( size_t j = 0; j < rows; ++j ) {
      tempdata[ i * rows + j ] = data( i, j );
      tempvar[ i * rows + j ] = invvar( i, j );
      tempmask[ i * rows + j ] = msk( i, j );
    }
  }

  write ( path, meta, rows, tempdata, tempvar, tempmask );

  return;
}


BOOST_CLASS_EXPORT(harp::image_desi)

image * harp::image_desi_create ( boost::property_tree::ptree const & props ) {
  return new image_desi ( props );
}

