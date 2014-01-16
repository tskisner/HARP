// @COPYRIGHT@

#include <harp_internal.hpp>

#include <boost/optional.hpp>

#include <cmath>

using namespace std;
using namespace harp;


void harp::psf_gauss_resp::sample ( size_t x_offset, size_t y_offset, matrix_double & patch ) const {

  double PI = std::atan2 ( 0.0, -1.0 );
  
  double cang = cos ( ang );
  double sang = sin ( ang );
  
  double invmaj = 1.0 / maj;
  double invmin = 1.0 / min;
  
  double xt, yt, exparg;

  double xrel, yrel;
  
  double norm = 0.0;

  double val;

  for ( size_t patch_col = 0; patch_col < patch.size2(); ++patch_col ) {

    for ( size_t patch_row = 0; patch_row < patch.size1(); ++patch_row ) {

      xrel = (double)( x_offset + patch_col ) - x;
      yrel = (double)( y_offset + patch_row ) - y;

      xt = xrel * cang + yrel * sang;
      yt = - xrel * sang + yrel * cang;
      exparg = - 0.5 * ( xt * xt * invmaj * invmaj + yt * yt * invmin * invmin );
      val = amp * exp ( exparg );
      norm += val;

      patch ( patch_row, patch_col ) = val;

    }

  }

  norm = 1.0 / norm;

  for ( size_t patch_col = 0; patch_col < patch.size2(); ++patch_col ) {
    for ( size_t patch_row = 0; patch_row < patch.size1(); ++patch_row ) {

      val = patch ( patch_row, patch_col );
      val *= norm;
      patch ( patch_row, patch_col ) = val;

    }
  }

  return;
}


// parameters

static const char * psf_gauss_key_path = "path";
static const char * psf_gauss_key_corr = "corr";
static const char * psf_gauss_key_rows = "imgrows";
static const char * psf_gauss_key_cols = "imgcols";

static const char * psf_gauss_key_nspec = "n_spec";
static const char * psf_gauss_key_nlambda = "n_lambda";

static const char * psf_gauss_key_name = "EXTNAME";

// HDU names for writing

static const char * psf_gauss_hdu_x = "X";
static const char * psf_gauss_hdu_y = "Y";
static const char * psf_gauss_hdu_lambda = "LogLam";
static const char * psf_gauss_hdu_amp = "Amplitude";
static const char * psf_gauss_hdu_maj = "MajorAxis";
static const char * psf_gauss_hdu_min = "MinorAxis";
static const char * psf_gauss_hdu_ang = "Angle";


harp::psf_gauss::psf_gauss ( boost::property_tree::ptree const & props ) : psf ( props ) {

  path_ = props.get < string > ( psf_gauss_key_path, "" );

  pixcorr_ = props.get < int > ( psf_gauss_key_corr );

  rows_ = props.get < size_t > ( psf_gauss_key_rows );

  cols_ = props.get < size_t > ( psf_gauss_key_cols );

  if ( path_ == "" ) {
    // if a path is not specified, it probably means that the calling code is going to manually set
    // the elliptical gaussian response using psf_gauss::parameters().  In that case, we set up the
    // HDU order to the default.  In this case, the calling code must also specify the number
    // of spectra and wavelength points.

    nspec_ = props.get < size_t > ( psf_gauss_key_nspec );

    nlambda_ = props.get < size_t > ( psf_gauss_key_nlambda );

    hdus_[ psf_gauss_hdu_x ] = 1;
    hdus_[ psf_gauss_hdu_y ] = 2;
    hdus_[ psf_gauss_hdu_lambda ] = 3;
    hdus_[ psf_gauss_hdu_amp ] = 4;
    hdus_[ psf_gauss_hdu_maj ] = 5;
    hdus_[ psf_gauss_hdu_min ] = 6;
    hdus_[ psf_gauss_hdu_ang ] = 7;

    // global dimensions

    nglobal_ = nspec_ * nlambda_;

    npix_ = rows_ * cols_;

    resp_.resize ( nglobal_ );

    lambda_.resize ( nlambda_ );

  } else {

    // look up HDU indices

    fitsfile *fp;

    fits::open_read ( fp, path_ );

    int hdu = fits::img_seek ( fp, psf_gauss_key_name, psf_gauss_hdu_x );
    hdus_[ psf_gauss_hdu_x ] = hdu;
    fits::img_dims ( fp, nspec_, nlambda_ );

    hdus_[ psf_gauss_hdu_y ] = hdu_info ( fp, psf_gauss_hdu_y );
    hdus_[ psf_gauss_hdu_lambda ] = hdu_info ( fp, psf_gauss_hdu_lambda );
    hdus_[ psf_gauss_hdu_amp ] = hdu_info ( fp, psf_gauss_hdu_amp );
    hdus_[ psf_gauss_hdu_maj ] = hdu_info ( fp, psf_gauss_hdu_maj );
    hdus_[ psf_gauss_hdu_min ] = hdu_info ( fp, psf_gauss_hdu_min );
    hdus_[ psf_gauss_hdu_ang ] = hdu_info ( fp, psf_gauss_hdu_ang );

    // global dimensions

    nglobal_ = nspec_ * nlambda_;

    npix_ = rows_ * cols_;

    // read the data

    resp_.resize ( nglobal_ );

    // verify that the wavelength solution for all spectra is equal
    // this is a requirement for HARP in order to allow for operations
    // that work simultaneously across all spectra.

    lambda_.resize ( nlambda_ );

    vector_double buffer;

    fits::img_seek ( fp, hdus_[ psf_gauss_hdu_lambda ] );      
    fits::img_read ( fp, buffer );

    for ( size_t j = 0; j < nlambda_; ++j ) {
      lambda_[j] = buffer[j];
      resp_[j].lambda = lambda_[j];
    }

    for ( size_t i = 1; i < nspec_; ++i ) {
      for ( size_t j = 0; j < nlambda_; ++j ) {
        if ( fabs ( buffer[ i * nlambda_ + j ] - lambda_[j] ) / lambda_[j] > 1.0e-6 ) {
          HARP_THROW( "wavelength solution is not constant across all spectra" );
        }
        resp_[ i * nlambda_ + j ].lambda = lambda_[j];
      }
    }

    fits::img_seek ( fp, hdus_[ psf_gauss_hdu_x ] );      
    fits::img_read ( fp, buffer );

    for ( size_t i = 0; i < nglobal_; ++i ) {
      resp_[i].x = buffer[i];
    }

    fits::img_seek ( fp, hdus_[ psf_gauss_hdu_y ] );      
    fits::img_read ( fp, buffer );

    for ( size_t i = 0; i < nglobal_; ++i ) {
      resp_[i].y = buffer[i];
    }

    fits::img_seek ( fp, hdus_[ psf_gauss_hdu_amp ] );
    fits::img_read ( fp, buffer );

    for ( size_t i = 0; i < nglobal_; ++i ) {
      resp_[i].amp = buffer[i];
    }

    fits::img_seek ( fp, hdus_[ psf_gauss_hdu_maj ] );      
    fits::img_read ( fp, buffer );

    for ( size_t i = 0; i < nglobal_; ++i ) {
      resp_[i].maj = buffer[i];
    }

    fits::img_seek ( fp, hdus_[ psf_gauss_hdu_min ] );      
    fits::img_read ( fp, buffer );

    for ( size_t i = 0; i < nglobal_; ++i ) {
      resp_[i].min = buffer[i];
    }

    fits::img_seek ( fp, hdus_[ psf_gauss_hdu_ang ] );      
    fits::img_read ( fp, buffer );

    for ( size_t i = 0; i < nglobal_; ++i ) {
      resp_[i].ang = buffer[i];
    }

    fits::close ( fp );

  }

}


harp::psf_gauss::~psf_gauss ( ) {
  
}


int harp::psf_gauss::hdu_info ( fitsfile *fp, const char * psf_gauss_hdu ) {
  int hdu = fits::img_seek ( fp, psf_gauss_key_name, psf_gauss_hdu );
  if ( hdu < 1 ) {
    std::ostringstream o;
    o << "could not find HDU named \"" << psf_gauss_hdu << "\"";
    HARP_THROW( o.str().c_str() );
  }
  size_t rows, cols;
  fits::img_dims ( fp, rows, cols );
  if ( ( rows != nspec_ ) || ( cols != nlambda_ ) ) {
    HARP_THROW( "psf_gauss: PSF file must have identical dimensions for all HDUs" );
  }
  return hdu;
}


boost::property_tree::ptree harp::psf_gauss::metadata ( ) const {

  return boost::property_tree::ptree();
}


size_t harp::psf_gauss::response_nnz_estimate ( ) const {
  size_t est = 2 * pixcorr_ + 1;
  est *= est;
  return est;
}


void harp::psf_gauss::response ( size_t spec, size_t lambda, size_t & x_offset, size_t & y_offset, matrix_double & patch ) const {

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

  size_t x_size;
  size_t y_size;

  if ( xmax > (double)(cols_ - 1) ) {
    x_size = ( cols_ - 1 ) - x_offset;
  } else {
    x_size = (size_t)xmax - x_offset;
  }

  if ( ymax > (double)(rows_ - 1) ) {
    y_size = ( rows_ - 1 ) - y_offset;
  } else {
    y_size = (size_t)ymax - y_offset;
  }

  patch.resize ( y_size, x_size );
  patch.clear();

  resp_[ bin ].sample ( x_offset, y_offset, patch );

  return;
}


void harp::psf_gauss::write ( std::string const & path ) {

  fitsfile *fp;

  fits::create ( fp, path );

  vector_double buffer ( nglobal_ );

  for ( size_t i = 0; i < nglobal_; ++i ) {
    buffer[i] = resp_[i].x;
  }
  fits::img_append < double > ( fp, nspec_, nlambda_ );
  fits::key_write ( fp, psf_gauss_key_name, psf_gauss_hdu_x, "" );
  fits::img_write ( fp, buffer );

  for ( size_t i = 0; i < nglobal_; ++i ) {
    buffer[i] = resp_[i].y;
  }
  fits::img_append < double > ( fp, nspec_, nlambda_ );
  fits::key_write ( fp, psf_gauss_key_name, psf_gauss_hdu_y, "" );
  fits::img_write ( fp, buffer );

  for ( size_t i = 0; i < nglobal_; ++i ) {
    buffer[i] = resp_[i].lambda;
  }
  fits::img_append < double > ( fp, nspec_, nlambda_ );
  fits::key_write ( fp, psf_gauss_key_name, psf_gauss_hdu_lambda, "" );
  fits::img_write ( fp, buffer );

  for ( size_t i = 0; i < nglobal_; ++i ) {
    buffer[i] = resp_[i].amp;
  }
  fits::img_append < double > ( fp, nspec_, nlambda_ );
  fits::key_write ( fp, psf_gauss_key_name, psf_gauss_hdu_amp, "" );
  fits::img_write ( fp, buffer );

  for ( size_t i = 0; i < nglobal_; ++i ) {
    buffer[i] = resp_[i].maj;
  }
  fits::img_append < double > ( fp, nspec_, nlambda_ );
  fits::key_write ( fp, psf_gauss_key_name, psf_gauss_hdu_maj, "" );
  fits::img_write ( fp, buffer );

  for ( size_t i = 0; i < nglobal_; ++i ) {
    buffer[i] = resp_[i].min;
  }
  fits::img_append < double > ( fp, nspec_, nlambda_ );
  fits::key_write ( fp, psf_gauss_key_name, psf_gauss_hdu_min, "" );
  fits::img_write ( fp, buffer );

  for ( size_t i = 0; i < nglobal_; ++i ) {
    buffer[i] = resp_[i].ang;
  }
  fits::img_append < double > ( fp, nspec_, nlambda_ );
  fits::key_write ( fp, psf_gauss_key_name, psf_gauss_hdu_ang, "" );
  fits::img_write ( fp, buffer );

  fits::close ( fp );

  return;
}





