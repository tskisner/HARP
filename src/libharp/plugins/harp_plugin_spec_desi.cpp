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

  fits::key_read ( fp, string(spec_desi_key_crval), crval );

  fits::key_read ( fp, string(spec_desi_key_cdelt), cdelt );

  fits::key_read ( fp, string(spec_desi_key_airorvac), airorvac );

  fits::key_read ( fp, string(spec_desi_key_loglam), loglam );

  fits::key_read ( fp, string(spec_desi_key_simfile), simfile );

  fits::key_read ( fp, string(spec_desi_key_camera), camera );

  fits::key_read ( fp, string(spec_desi_key_exptime), exptime );

  fits::key_read ( fp, string(spec_desi_key_rdnoise), rdnoise );

  fits::key_read ( fp, string(spec_desi_key_flavor), flavor );

  fits::key_read ( fp, string(spec_desi_key_inpsf), in_psf );

  fits::key_read ( fp, string(spec_desi_key_inimg), in_img );

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


void harp::spec_desi::write ( std::string const & path, boost::property_tree::ptree const & meta, vector_double const & data, vector_double const & invvar, vector_double const & lambda ) {

  size_t nlambda = lambda.size();
  size_t nspec = (size_t)( data.size() / nlambda );

  fitsfile * fp;
    
  fits::create ( fp, path );

  fits::img_append < double > ( fp, nspec, nlambda );

  if ( meta.count ( spec_desi_key_crval ) ) {
    fits::key_write ( fp, spec_desi_key_crval, meta.get < float > ( spec_desi_key_crval ), "Starting wavelength [Angstroms]" );
  } else {
    HARP_THROW( "you must specify the starting wavelength" );
  }

  if ( meta.count ( spec_desi_key_cdelt ) ) {
    fits::key_write ( fp, spec_desi_key_cdelt, meta.get < float > ( spec_desi_key_cdelt ), "Wavelength step [Angstroms]" );
  } else {
    HARP_THROW( "you must specify the wavelength step" );
  }

  if ( meta.count ( spec_desi_key_airorvac ) ) {
    fits::key_write ( fp, spec_desi_key_airorvac, meta.get < string > ( spec_desi_key_airorvac ), "Air or Vacuum wavelengths" );
  } else {
    HARP_THROW( "you must specify air or vac" );
  }

  if ( meta.count ( spec_desi_key_loglam ) ) {
    fits::key_write ( fp, spec_desi_key_loglam, meta.get < int > ( spec_desi_key_loglam ), "linear wavelength steps, not log10" );
  } else {
    HARP_THROW( "you must specify log lambda" );
  }

  if ( meta.count ( spec_desi_key_simfile ) ) {
    fits::key_write ( fp, spec_desi_key_simfile, meta.get < string > ( spec_desi_key_simfile ), "Input simulation file" );
  }

  if ( meta.count ( spec_desi_key_camera ) ) {
    fits::key_write ( fp, spec_desi_key_camera, meta.get < string > ( spec_desi_key_camera ), "Spectograph Camera" );
  } else {
    HARP_THROW( "you must specify the camera name" );
  }

  fits::key_write ( fp, spec_desi_key_harp, source_version(), "HARP Version" );

  if ( meta.count ( spec_desi_key_exptime ) ) {
    fits::key_write ( fp, spec_desi_key_exptime, meta.get < float > ( spec_desi_key_exptime ), "Exposure time [sec]" );
  } else {
    HARP_THROW( "you must specify the exposure time" );
  }

  if ( meta.count ( spec_desi_key_rdnoise ) ) {
    fits::key_write ( fp, spec_desi_key_rdnoise, meta.get < float > ( spec_desi_key_rdnoise ), "Read noise [electrons]" );
  } else {
    HARP_THROW( "you must specify the read noise" );
  }

  if ( meta.count ( spec_desi_key_flavor ) ) {
    fits::key_write ( fp, spec_desi_key_flavor, meta.get < string > ( spec_desi_key_flavor ), "Exposure type (arc, flat, science)" );
  } else {
    HARP_THROW( "you must specify the exposure flavor" );
  }

  if ( meta.count ( spec_desi_key_inpsf ) ) {
    fits::key_write ( fp, spec_desi_key_inpsf, meta.get < string > ( spec_desi_key_inpsf ), "Input spectral PSF" );
  } else {
    HARP_THROW( "you must specify the input PSF" );
  }

  if ( meta.count ( spec_desi_key_inimg ) ) {
    fits::key_write ( fp, spec_desi_key_inimg, meta.get < string > ( spec_desi_key_inimg ), "Input image" );
  } else {
    HARP_THROW( "you must specify the input image" );
  }


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
