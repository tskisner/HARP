/*
  High Performance Astrophysical Reconstruction and Processing (HARP)
  (c) 2014-2015, The Regents of the University of California, 
  through Lawrence Berkeley National Laboratory.  See top
  level LICENSE file for details.
*/


#include <harp/static_plugins.hpp>

using namespace std;
using namespace harp;


static const char * spec_desisim_key_path = "path";
static const char * spec_desisim_key_crval = "CRVAL1";
static const char * spec_desisim_key_cdelt = "CDELT1";
static const char * spec_desisim_key_airorvac = "AIRORVAC";
static const char * spec_desisim_key_loglam = "LOGLAM";


harp::spec_desisim::spec_desisim ( ) : spec () {
  nspec_ = 0;
  nlambda_ = 0;
  nglobal_ = 0;
  path_ = "";
  objhdu_ = 1;
  skyhdu_ = 2;
  crval = 0.0;
  cdelt = 0.0;
  airorvac = "vac";
  loglam = 0;
}


harp::spec_desisim::spec_desisim ( boost::property_tree::ptree const & props ) : spec ( "desi", props ) {

  path_ = props.get < string > ( spec_desisim_key_path, "" );

  objhdu_ = 1;

  skyhdu_ = 2;

  fitsfile * fp;

  fits::open_read ( fp, path_ );

  fits::img_seek ( fp, objhdu_ );

  fits::img_dims ( fp, nspec_, nlambda_ );

  // read keywords

  meta_ = fits::key_read_all ( fp );

  fits::key_parse ( meta_, spec_desisim_key_crval, crval );

  fits::key_parse ( meta_, spec_desisim_key_cdelt, cdelt );

  fits::key_parse ( meta_, spec_desisim_key_airorvac, airorvac );

  fits::key_parse ( meta_, spec_desisim_key_loglam, loglam );

  fits::close ( fp );

  nglobal_ = nspec_ * nlambda_;
  
}


harp::spec_desisim::~spec_desisim ( ) {
  
}


size_t harp::spec_desisim::n_spec ( ) const {
  return nspec_;
}


size_t harp::spec_desisim::n_lambda ( ) const {
  return nlambda_;
}


void harp::spec_desisim::values ( vector_double & data ) const {

  data.resize ( nglobal_ );
  data.clear();

  fitsfile * fp;

  fits::open_read ( fp, path_ );

  // read the object flux

  fits::img_seek ( fp, objhdu_ );

  fits::img_read ( fp, data, false );

  // read the sky flux and sum

  vector_double skyflux ( data.size() );

  fits::img_seek ( fp, skyhdu_ );

  fits::img_read ( fp, skyflux, false );

  fits::close ( fp );

  for ( size_t i = 0; i < data.size(); ++i ) {
    data[i] += skyflux[i];
  }

  return;
}


void harp::spec_desisim::sky ( vector_double & data ) const {

  data.resize ( nglobal_ );
  data.clear();

  fitsfile * fp;

  fits::open_read ( fp, path_ );

  // read the sky flux

  fits::img_seek ( fp, skyhdu_ );

  fits::img_read ( fp, data, false );

  fits::close ( fp );

  return;
}


void harp::spec_desisim::inv_variance ( vector_double & data ) const {

  data.resize ( nglobal_ );
  data.clear();

  return;
}


void harp::spec_desisim::lambda ( vector_double & lambda_vals ) const {

  lambda_vals.resize ( nlambda_ );

  for ( size_t i = 0; i < nlambda_; ++i ) {
    lambda_vals[i] = crval + cdelt * (double)i;
  }

  return;
}


boost::property_tree::ptree harp::spec_desisim::meta () const {
  return meta_;
}



BOOST_CLASS_EXPORT(harp::spec_desisim)

spec * harp::spec_desisim_create ( boost::property_tree::ptree const & props ) {
  return new spec_desisim ( props );
}
