// @COPYRIGHT@


#include <harp/static_plugins.hpp>


using namespace std;
using namespace harp;


// parameters

static const char * psf_gauss_sim_key_corr = "corr";
static const char * psf_gauss_sim_key_rows = "imgrows";
static const char * psf_gauss_sim_key_cols = "imgcols";

static const char * psf_gauss_sim_key_lambda_spec = "lambda_spec";
static const char * psf_gauss_sim_key_lambda_spec_type = "lambda_spec_type";
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


harp::psf_gauss_sim::psf_gauss_sim ( boost::property_tree::ptree const & props ) : psf ( "gauss_sim", props ) {

  // check to see if we are getting the wavelength solution from a spec

  if ( props.count ( psf_gauss_sim_key_lambda_spec ) ) {

    lambda_spec_props_ = props.get_child ( psf_gauss_sim_key_lambda_spec );

    lambda_spec_type_ = props.get < string > ( psf_gauss_sim_key_lambda_spec_type );

    plugin_registry & reg = plugin_registry::get();

    spec_p child_spec ( reg.create_spec ( lambda_spec_type_, lambda_spec_props_ ) );

    nlambda_ = child_spec->n_lambda();

    child_spec->lambda ( lambda_ );

    first_lambda_ = lambda_[0];
    last_lambda_ = lambda_[ nlambda_ - 1 ];

  } else {
    // in this case, the props must specify the number of lambda points and spacing

    last_lambda_ = props.get < double > ( psf_gauss_sim_key_lambda_stop, 9808.0 );

    first_lambda_ = props.get < double > ( psf_gauss_sim_key_lambda_start, 7460.0 );

    nlambda_ = props.get < size_t > ( psf_gauss_sim_key_lambda_n, 4697 );

    lambda_.resize ( nlambda_ );
    double incr = ( last_lambda_ - first_lambda_ ) / (double)( nlambda_ - 1 );
    for ( size_t i = 0; i < nlambda_; ++i ) {
      lambda_[i] = incr * (double)i;
    }

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

      spec2pix ( spec, lambda, xpix, ypix );

      resp_[ bin ].x = xpix;
      resp_[ bin ].y = ypix;
      resp_[ bin ].amp = 1.0;
      resp_[ bin ].lambda = lambda_[ lambda ];
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


size_t harp::psf_gauss_sim::response_nnz_estimate ( ) const {
  size_t est = 2 * pixcorr_ + 1;
  est *= est;
  return est;
}



void harp::psf_gauss_sim::extent ( size_t spec, size_t lambda, size_t & x_offset, size_t & y_offset, size_t & n_x, size_t & n_y ) const {

  size_t bin = spec * nlambda_ + lambda;

  double xmin = resp_[ bin ].x - (double)pixcorr_;
  double xmax = resp_[ bin ].x + (double)pixcorr_;

  double ymin = resp_[ bin ].y - (double)pixcorr_;
  double ymax = resp_[ bin ].y + (double)pixcorr_;

  if ( xmin < 0.0 ) {
    x_offset = 0;
  } else {
    x_offset = (size_t)xmin;
  }

  if ( ymin < 0.0 ) {
    y_offset = 0;
  } else {
    y_offset = (size_t)ymin;
  }

  if ( xmax > (double)(cols_ - 1) ) {
    n_x = ( cols_ - 1 ) - x_offset + 1;
  } else {
    n_x = (size_t)xmax - x_offset + 1;
  }

  if ( ymax > (double)(rows_ - 1) ) {
    n_y = ( rows_ - 1 ) - y_offset + 1;
  } else {
    n_y = (size_t)ymax - y_offset + 1;
  }

  return;
}


void harp::psf_gauss_sim::response ( size_t spec, size_t lambda, size_t & x_offset, size_t & y_offset, matrix_double & patch ) const {

  size_t x_size;
  size_t y_size;

  extent ( spec, lambda, x_offset, y_offset, x_size, y_size );

  size_t bin = spec * nlambda_ + lambda;

  patch.resize ( y_size, x_size );
  patch.clear();

  resp_[ bin ].sample ( x_offset, y_offset, patch );

  return;
}


void harp::psf_gauss_sim::write ( std::string const & path ) {

  fitsfile *fp;

  fits::create ( fp, path );

  vector_double buffer ( nglobal_ );

  for ( size_t i = 0; i < nglobal_; ++i ) {
    buffer[i] = resp_[i].x;
  }
  fits::img_append < double > ( fp, nspec_, nlambda_ );
  fits::key_write ( fp, "EXTNAME", psf_gauss_sim_hdu_x, "" );
  fits::img_write ( fp, buffer, false );

  for ( size_t i = 0; i < nglobal_; ++i ) {
    buffer[i] = resp_[i].y;
  }
  fits::img_append < double > ( fp, nspec_, nlambda_ );
  fits::key_write ( fp, "EXTNAME", psf_gauss_sim_hdu_y, "" );
  fits::img_write ( fp, buffer, false );

  for ( size_t i = 0; i < nglobal_; ++i ) {
    buffer[i] = resp_[i].lambda;
  }
  fits::img_append < double > ( fp, nspec_, nlambda_ );
  fits::key_write ( fp, "EXTNAME", psf_gauss_sim_hdu_lambda, "" );
  fits::img_write ( fp, buffer, false );

  for ( size_t i = 0; i < nglobal_; ++i ) {
    buffer[i] = resp_[i].amp;
  }
  fits::img_append < double > ( fp, nspec_, nlambda_ );
  fits::key_write ( fp, "EXTNAME", psf_gauss_sim_hdu_amp, "" );
  fits::img_write ( fp, buffer, false );

  for ( size_t i = 0; i < nglobal_; ++i ) {
    buffer[i] = resp_[i].maj;
  }
  fits::img_append < double > ( fp, nspec_, nlambda_ );
  fits::key_write ( fp, "EXTNAME", psf_gauss_sim_hdu_maj, "" );
  fits::img_write ( fp, buffer, false );

  for ( size_t i = 0; i < nglobal_; ++i ) {
    buffer[i] = resp_[i].min;
  }
  fits::img_append < double > ( fp, nspec_, nlambda_ );
  fits::key_write ( fp, "EXTNAME", psf_gauss_sim_hdu_min, "" );
  fits::img_write ( fp, buffer, false );

  for ( size_t i = 0; i < nglobal_; ++i ) {
    buffer[i] = resp_[i].ang;
  }
  fits::img_append < double > ( fp, nspec_, nlambda_ );
  fits::key_write ( fp, "EXTNAME", psf_gauss_sim_hdu_ang, "" );
  fits::img_write ( fp, buffer, false );

  fits::close ( fp );

  return;
}


BOOST_CLASS_EXPORT(harp::psf_gauss_sim)

psf * harp::psf_gauss_sim_create ( boost::property_tree::ptree const & props ) {
  return new psf_gauss_sim ( props );
}

