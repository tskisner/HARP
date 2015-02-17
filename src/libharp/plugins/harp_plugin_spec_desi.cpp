/*
  High Performance Astrophysical Reconstruction and Processing (HARP)
  (c) 2014-2015, The Regents of the University of California, 
  through Lawrence Berkeley National Laboratory.  See top
  level LICENSE file for details.
*/


#include <harp/static_plugins.hpp>

using namespace std;
using namespace harp;


static const char * spec_desi_key_path = "path";
static const char * spec_desi_key_crval = "CRVAL1";
static const char * spec_desi_key_cdelt = "CDELT1";
static const char * spec_desi_key_airorvac = "AIRORVAC";
static const char * spec_desi_key_loglam = "LOGLAM";
static const char * spec_desi_key_simfile = "SIMFILE";
static const char * spec_desi_key_camera = "CAMERA";
static const char * spec_desi_key_harp = "HARP";
static const char * spec_desi_key_exptime = "EXPTIME";
static const char * spec_desi_key_rdnoise = "RDNOISE";
static const char * spec_desi_key_flavor = "FLAVOR";
static const char * spec_desi_key_inpsf = "IN_PSF";
static const char * spec_desi_key_inimg = "IN_IMG";


harp::spec_desi::spec_desi ( ) : spec () {
  nspec_ = 0;
  nlambda_ = 0;
  nglobal_ = 0;
  path_ = "";
  spechdu_ = 1;
  invvarhdu_ = 2;
  lambdahdu_ = 3;
  crval = 0.0;
  cdelt = 0.0;
  airorvac = "vac";
  loglam = 0;
  simfile = "";
  camera = "";
  exptime = 0.0;
  rdnoise = 0.0;
  flavor = "";
  in_psf = "";
  in_img = "";
}


harp::spec_desi::spec_desi ( boost::property_tree::ptree const & props ) : spec ( "desi", props ) {

  path_ = props.get < string > ( spec_desi_key_path, "" );

  spechdu_ = 1;

  invvarhdu_ = 2;

  lambdahdu_ = 3;

  fitsfile * fp;

  fits::open_read ( fp, path_ );

  fits::img_seek ( fp, spechdu_ );

  fits::img_dims ( fp, nspec_, nlambda_ );

  // read keywords

  meta_ = fits::key_read_all ( fp );

  fits::key_parse ( meta_, spec_desi_key_crval, crval );

  fits::key_parse ( meta_, spec_desi_key_cdelt, cdelt );

  fits::key_parse ( meta_, spec_desi_key_airorvac, airorvac );

  fits::key_parse ( meta_, spec_desi_key_loglam, loglam );

  fits::key_parse ( meta_, spec_desi_key_simfile, simfile );

  fits::key_parse ( meta_, spec_desi_key_camera, camera );

  fits::key_parse ( meta_, spec_desi_key_exptime, exptime );

  fits::key_parse ( meta_, spec_desi_key_rdnoise, rdnoise );

  fits::key_parse ( meta_, spec_desi_key_flavor, flavor );
  
  fits::key_parse ( meta_, spec_desi_key_inpsf, in_psf );
  
  fits::key_parse ( meta_, spec_desi_key_inimg, in_img );

  fits::close ( fp );

  nglobal_ = nspec_ * nlambda_;
  
}


harp::spec_desi::~spec_desi ( ) {
  
}


size_t harp::spec_desi::n_spec ( ) const {
  return nspec_;
}


size_t harp::spec_desi::n_lambda ( ) const {
  return nlambda_;
}


void harp::spec_desi::values ( vector_double & data ) const {

  data.resize ( nglobal_ );
  data.clear();

  fitsfile * fp;

  fits::open_read ( fp, path_ );

  // read the spectral data

  fits::img_seek ( fp, spechdu_ );

  fits::img_read ( fp, data, false );

  fits::close ( fp );

  return;
}


void harp::spec_desi::inv_variance ( vector_double & data ) const {

  data.resize ( nglobal_ );
  data.clear();

  fitsfile * fp;

  fits::open_read ( fp, path_ );

  // read the spectral data

  fits::img_seek ( fp, invvarhdu_ );

  fits::img_read ( fp, data, false );

  fits::close ( fp );

  return;
}


void harp::spec_desi::lambda ( vector_double & lambda_vals ) const {

  lambda_vals.resize ( nlambda_ );

  fitsfile * fp;

  fits::open_read ( fp, path_ );

  // read the wavelength vector

  fits::img_seek ( fp, lambdahdu_ );
  fits::img_read ( fp, lambda_vals, false );

  fits::close ( fp );

  return;
}


boost::property_tree::ptree harp::spec_desi::meta () const {
  return meta_;
}


void harp::spec_desi::write ( std::string const & path, boost::property_tree::ptree const & meta, vector_double const & data, vector_double const & invvar, vector_double const & lambda ) {

  size_t nlambda = lambda.size();
  size_t nspec = (size_t)( data.size() / nlambda );

  fitsfile * fp;
    
  fits::create ( fp, path );

  fits::img_append < double > ( fp, nspec, nlambda );

  // check required keywords and write

  if ( meta.count ( spec_desi_key_crval ) ) {
    fits::key_require ( meta, spec_desi_key_crval, "F" );
  } else {
    HARP_THROW( "you must specify the starting wavelength" );
  }

  if ( meta.count ( spec_desi_key_cdelt ) ) {
    fits::key_require ( meta, spec_desi_key_cdelt, "F" );
  } else {
    HARP_THROW( "you must specify the wavelength step" );
  }

  if ( meta.count ( spec_desi_key_airorvac ) ) {
    fits::key_require ( meta, spec_desi_key_airorvac, "C" );
  } else {
    HARP_THROW( "you must specify air or vac" );
  }

  if ( meta.count ( spec_desi_key_loglam ) ) {
    fits::key_require ( meta, spec_desi_key_loglam, "I" );
  } else {
    HARP_THROW( "you must specify log lambda" );
  }

  if ( meta.count ( spec_desi_key_camera ) ) {
    fits::key_require ( meta, spec_desi_key_camera, "C" );
  } else {
    HARP_THROW( "you must specify the camera name" );
  }

  if ( meta.count ( spec_desi_key_exptime ) ) {
    fits::key_require ( meta, spec_desi_key_exptime, "F" );
  } else {
    HARP_THROW( "you must specify the exposure time" );
  }

  if ( meta.count ( spec_desi_key_rdnoise ) ) {
    fits::key_require ( meta, spec_desi_key_rdnoise, "F" );
  } else {
    HARP_THROW( "you must specify the read noise" );
  }

  if ( meta.count ( spec_desi_key_flavor ) ) {
    fits::key_require ( meta, spec_desi_key_flavor, "C" );
  } else {
    HARP_THROW( "you must specify the exposure flavor" );
  }

  if ( meta.count ( spec_desi_key_inpsf ) ) {
    fits::key_require ( meta, spec_desi_key_inpsf, "C" );
  } else {
    HARP_THROW( "you must specify the input PSF" );
  }

  if ( meta.count ( spec_desi_key_inimg ) ) {
    fits::key_require ( meta, spec_desi_key_inimg, "C" );
  } else {
    HARP_THROW( "you must specify the input image" );
  }

  fits::key_write ( fp, spec_desi_key_harp, source_version(), "HARP Version" );

  fits::key_write_all ( fp, meta );


  fits::img_write ( fp, data, false );


  fits::img_append < double > ( fp, nspec, nlambda );

  fits::key_write ( fp, "EXTNAME", string("IVAR"), "" );
  
  fits::img_write ( fp, invvar, false );


  fits::img_append < double > ( fp, 1, nlambda );

  fits::key_write ( fp, "EXTNAME", string("WAVELENGTH"), "" );

  fits::img_write ( fp, lambda, false );

  fits::close ( fp );

  return;
}


BOOST_CLASS_EXPORT(harp::spec_desi)

spec * harp::spec_desi_create ( boost::property_tree::ptree const & props ) {
  return new spec_desi ( props );
}
