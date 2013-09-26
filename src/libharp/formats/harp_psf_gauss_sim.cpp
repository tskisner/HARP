// @COPYRIGHT@

#include <harp_internal.hpp>

#include <boost/optional.hpp>

#include <cmath>

using namespace std;
using namespace harp;


// parameters

static const char * psf_gauss_sim_key_corr = "corr";
static const char * psf_gauss_sim_key_rows = "imgrows";
static const char * psf_gauss_sim_key_cols = "imgcols";

static const char * psf_gauss_sim_key_lambda_spec = "lambda_spec";
static const char * psf_gauss_sim_key_lambda_n = "lambda_n";
static const char * psf_gauss_sim_key_lambda_start = "lambda_start";
static const char * psf_gauss_sim_key_lambda_stop = "lambda_stop";

static const char * psf_gauss_sim_key_bundle_size = "bundle_size";
static const char * psf_gauss_sim_key_nbundle = "nbundle";
static const char * psf_gauss_sim_key_fwhm = "fwhm";
static const char * psf_gauss_sim_key_margin = "margin";
static const char * psf_gauss_sim_key_gap = "gap";

// HDU names for writing

static const char * psf_gauss_sim_hdu_x = "X";
static const char * psf_gauss_sim_hdu_y = "Y";
static const char * psf_gauss_sim_hdu_lambda = "LogLam";
static const char * psf_gauss_sim_hdu_amp = "Amplitude";
static const char * psf_gauss_sim_hdu_maj = "MajorAxis";
static const char * psf_gauss_sim_hdu_min = "MinorAxis";
static const char * psf_gauss_sim_hdu_ang = "Angle";


harp::psf_gauss_sim::psf_gauss_sim ( boost::property_tree::ptree const & props ) : psf ( props ) {

  // check to see if we are getting the wavelength solution from a spec

  if ( props.count ( psf_gauss_sim_key_lambda_spec ) ) {

    lambda_spec_props_ = props.get_child ( psf_gauss_sim_key_lambda_spec );

    spec_p child_spec ( spec::create ( lambda_spec_props_ ) );

    nlambda_ = child_spec->n_lambda();

    vector_double child_spec_data;
    std::vector < bool > child_spec_sky;

    child_spec->read ( child_spec_data, lambda_, child_spec_sky );

    first_lambda_ = lambda_[0];
    last_lambda_ = lambda_[ nlambda_ - 1 ];

  } else {
    // in this case, the props must specify the number of lambda points and spacing

    last_lambda_ = props.get < double > ( psf_gauss_sim_key_lambda_stop, 9808.0 );

    first_lambda_ = props.get < double > ( psf_gauss_sim_key_lambda_start, 7460.0 );

    nlambda_ = props.get < size_t > ( psf_gauss_sim_key_lambda_n, 4697 );

  }

  // fiber spacing

  n_bundle_ = props.get < size_t > ( psf_gauss_sim_key_nbundle, 20 );

  bundle_size_ = props.get < size_t > ( psf_gauss_sim_key_bundle_size, 25 );

  nspec_ = n_bundle_ * bundle_size_;

  pix_margin_ = props.get < double > ( psf_gauss_sim_key_margin, 10.0 );

  pix_gap_ = props.get < double > ( psf_gauss_sim_key_gap, 7.0 );

  pix_bundle_ = 2.0 * pix_margin_ + (double)(bundle_size_ - 1) * pix_gap_ + (double)bundle_size_;

  // response fwhm

  psf_fwhm_ = props.get < double > ( psf_gauss_sim_key_fwhm, 2.2 );

  // response correlation length in pixels

  pixcorr_ = props.get < size_t > ( psf_gauss_sim_key_corr, 10 );

  cols_ = (size_t)( pix_bundle_ * (double)n_bundle_ + 1.0 );
  pix_offset_ = 0.5 * ( (double)cols_ - ( pix_bundle_ * (double)n_bundle_ ) );

  rows_ = nlambda_ + 2 * pixcorr_;

  // HDU list for writing
  
  hdus_[ psf_gauss_sim_hdu_x ] = 1;
  hdus_[ psf_gauss_sim_hdu_y ] = 2;
  hdus_[ psf_gauss_sim_hdu_lambda ] = 3;
  hdus_[ psf_gauss_sim_hdu_amp ] = 4;
  hdus_[ psf_gauss_sim_hdu_maj ] = 5;
  hdus_[ psf_gauss_sim_hdu_min ] = 6;
  hdus_[ psf_gauss_sim_hdu_ang ] = 7;
 
  // generate response parameters

  nglobal_ = nspec_ * nlambda_;

  npix_ = rows_ * cols_;

  resp_.resize ( nglobal_ );

  double xpix;
  double ypix;

  for ( size_t spec = 0; spec < nspec_; ++spec ) {
    
    for ( size_t lambda = 0; lambda < nlambda_; ++lambda ) {

      size_t bin = spec * nlambda_ + lambda;

      spec2pix ( spec, lambda, ypix, xpix );

      resp_[ bin ].x = xpix;
      resp_[ bin ].y = ypix;
      resp_[ bin ].lambda = lambda_[ lambda ];
      resp_[ bin ].amp = 1.0;
      resp_[ bin ].maj = psf_fwhm_ / 2.0;
      resp_[ bin ].min = psf_fwhm_ / 2.0;
      resp_[ bin ].ang = 0.0;

    }

  }

}


harp::psf_gauss_sim::~psf_gauss_sim ( ) {
  
}


void harp::psf_gauss_sim::spec2pix ( size_t spec, size_t specbin, double & x, double & y ) {

  size_t bundle = (size_t)( spec / bundle_size_ );
  size_t bundle_spec = spec - bundle * bundle_size_;

  x = pix_offset_ + ( (double)bundle * pix_bundle_ ) + pix_margin_ + ( (double)bundle_spec * ( pix_gap_ + 1.0 ) );

  y = (double)( specbin + pixcorr_ );

  return;
}


void harp::psf_gauss_sim::response ( size_t spec, size_t lambda, size_t & x_offset, size_t & y_offset, matrix_double & patch ) {

  size_t bin = spec * nlambda_ + lambda;

  resp_[ bin ].sample ( x_offset, y_offset, patch );

  return;
}




