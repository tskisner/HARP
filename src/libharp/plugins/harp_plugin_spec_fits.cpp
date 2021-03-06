/*
  High Performance Astrophysical Reconstruction and Processing (HARP)
  (c) 2014-2015, The Regents of the University of California, 
  through Lawrence Berkeley National Laboratory.  See top
  level LICENSE file for details.
*/


#include <harp/static_plugins.hpp>

using namespace std;
using namespace harp;


static const char * spec_fits_key_path = "path";
static const char * spec_fits_key_nspec = "nspec";
static const char * spec_fits_key_nlambda = "nlambda";
static const char * spec_fits_key_hduspec = "hduspec";
static const char * spec_fits_key_hduinvvar = "hduinvvar";
static const char * spec_fits_key_hdulambda = "hdulambda";


harp::spec_fits::spec_fits ( ) : spec () {
  nspec_ = 0;
  nlambda_ = 0;
  nglobal_ = 0;
  path_ = "";
  spechdu_ = -1;
  invvarhdu_ = -1;
  lambdahdu_ = -1;
}


harp::spec_fits::spec_fits ( boost::property_tree::ptree const & props ) : spec ( "fits", props ) {

  path_ = props.get < string > ( spec_fits_key_path, "" );

  spechdu_ = props.get < int > ( spec_fits_key_hduspec, 1 );

  invvarhdu_ = props.get < int > ( spec_fits_key_hduinvvar, 2 );

  lambdahdu_ = props.get < int > ( spec_fits_key_hdulambda, 3 );

  if ( path_ == "" ) {

    nspec_ = props.get < size_t > ( spec_fits_key_nspec );

    nlambda_ = props.get < size_t > ( spec_fits_key_nlambda );

  } else {

    fitsfile * fp;

    fits::open_read ( fp, path_ );

    fits::img_seek ( fp, spechdu_ );

    fits::img_dims ( fp, nspec_, nlambda_ );

    fits::close ( fp );

  }

  nglobal_ = nspec_ * nlambda_;
  
}


harp::spec_fits::~spec_fits ( ) {
  
}


size_t harp::spec_fits::n_spec ( ) const {
  return nspec_;
}


size_t harp::spec_fits::n_lambda ( ) const {
  return nlambda_;
}


void harp::spec_fits::values ( vector_double & data ) const {

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


void harp::spec_fits::inv_variance ( vector_double & data ) const {

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


void harp::spec_fits::lambda ( vector_double & lambda_vals ) const {

  lambda_vals.resize ( nlambda_ );

  fitsfile * fp;

  fits::open_read ( fp, path_ );

  // read the wavelength vector

  fits::img_seek ( fp, lambdahdu_ );
  fits::img_read ( fp, lambda_vals, false );

  fits::close ( fp );

  return;
}


void harp::spec_fits::write ( std::string const & path, vector_double const & data, vector_double const & invvar, vector_double const & lambda ) {

  size_t nlambda = lambda.size();
  size_t nspec = (size_t)( data.size() / nlambda );

  fitsfile * fp;
    
  fits::create ( fp, path );

  fits::img_append < double > ( fp, nspec, nlambda );
  fits::key_write ( fp, "EXTNAME", string("FLUX"), "" );
  fits::img_write ( fp, data, false );

  fits::img_append < double > ( fp, nspec, nlambda );
  fits::key_write ( fp, "EXTNAME", string("INV_VAR"), "" );
  fits::img_write ( fp, invvar, false );

  fits::img_append < double > ( fp, 1, nlambda );
  fits::key_write ( fp, "EXTNAME", string("WAVELENGTH"), "" );
  fits::img_write ( fp, lambda, false );

  fits::close ( fp );

  return;
}


BOOST_CLASS_EXPORT(harp::spec_fits)

spec * harp::spec_fits_create ( boost::property_tree::ptree const & props ) {
  return new spec_fits ( props );
}
