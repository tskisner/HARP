// @COPYRIGHT@

#include <harp_internal.hpp>

#include <boost/optional.hpp>

#include <cmath>

using namespace std;
using namespace harp;


void harp::psf_gauss_resp::sample ( matrix_local & vals, matrix_local & xrel, matrix_local & yrel ) {

  double PI = std::atan2 ( 0.0, -1.0 );
  
  amp /= maj_ * min_ * 2.0 * PI;
  
  size_t nvals = xrel.Height();
  
  double cang = cos ( ang );
  double sang = sin ( ang );
  
  double invmaj = 1.0 / maj;
  double invmin = 1.0 / min;
  
  double xt, yt, exparg;
  
  size_t i;

  double norm = 0.0;
  
  for ( i = 0; i < nvals; ++i ) {
    xt = xrel.Get( i, 0 ) * cang + yrel.Get( i, 0 ) * sang;
    yt = - xrel.Get( i, 0 ) * sang + yrel.Get( i, 0 ) * cang;
    exparg = - 0.5 * ( xt * xt * invmaj * invmaj + yt * yt * invmin * invmin );
    vals.Set( i, 0, amp * exp ( exparg ) );
    norm += vals.Get ( i, 0 );
  }

  norm = 1.0 / norm;

  double temp;

  for ( i = 0; i < nvals; ++i ) {
    temp = vals.Get ( i, 0 );
    temp *= norm;
    vals.Set ( i, 0, temp );
  }  
  
  return;
}


static const char * psf_gauss_key_path = "path";
static const char * psf_gauss_key_corr = "corr";
static const char * psf_gauss_key_rows = "imgrows";
static const char * psf_gauss_key_cols = "imgcols";

static const char * psf_gauss_key_name = "EXTNAME";

static const char * psf_gauss_hdu_x = "X";
static const char * psf_gauss_hdu_y = "Y";
static const char * psf_gauss_hdu_lambda = "LogLam";
static const char * psf_gauss_hdu_amp = "Amplitude";
static const char * psf_gauss_hdu_maj = "MajorAxis";
static const char * psf_gauss_hdu_min = "MinorAxis";
static const char * psf_gauss_hdu_ang = "Angle";

static const char * psf_gauss_key_fake = "fake";
static const char * psf_gauss_key_fake_spec = "spec";
static const char * psf_gauss_key_fake_bundle_size = "bundle_size";
static const char * psf_gauss_key_fake_nbundle = "nbundle";
static const char * psf_gauss_key_fake_fwhm = "fwhm";
static const char * psf_gauss_key_fake_margin = "margin";
static const char * psf_gauss_key_fake_gap = "gap";


harp::psf_gauss::psf_gauss ( boost::property_tree::ptree const & props ) : psf ( props ) {

  boost::optional < string > fakeval = props.get_optional < string > ( psf_gauss_key_fake );
  string fake = boost::get_optional_value_or ( fakeval, "FALSE" );
  dofake_ = ( fake == "TRUE" );

  if ( dofake_ ) {

    // in this case, we initialize the PSF to a symmetric gaussian

    path_ = "";

    fake_spec_props_ = props.get_child ( psf_gauss_key_fake_spec );

    spec_p child_spec ( spec::create ( fake_spec_props_ ) );
    size_t spec_nspec = child_spec->n_spec();

    nlambda_ = child_spec->n_lambda();

    fake_n_bundle_ = props.get < size_t > ( psf_gauss_key_fake_nbundle, 20 );
    fake_bundle_size_ = props.get < size_t > ( psf_gauss_key_fake_bundle_size, 25 );

    nspec_ = fake_n_bundle_ * fake_bundle_size_;

    if ( nspec_ != spec_nspec ) {
      std::ostringstream o;
      o << "number spectra in simulated PSF spec (" << spec_nspec << ") does not match the number from bundle size/count (" << nspec_ << ")";
      HARP_THROW( o.str().c_str() );
    }

    fake_pix_margin_ = props.get < size_t > ( psf_gauss_key_fake_margin, 10 );

    fake_pix_gap_ = props.get < size_t > ( psf_gauss_key_fake_gap, 7 );

    fake_pix_bundle_ = 2 * fake_pix_margin_ + (fake_bundle_size_ - 1) * fake_pix_gap_ + fake_bundle_size_;

    // response fwhm

    fake_psf_fwhm_ = props.get < double > ( psf_gauss_key_fake_fwhm, 2.2 );

    // response correlation length in pixels

    pixcorr_ = props.get < size_t > ( psf_gauss_key_corr, 10 );

    cols_ = fake_pix_bundle_ * fake_n_bundle_;
    rows_ = nlambda_ + 2 * pixcorr_;

    // wavelength solution

    lambda_.resize ( nlambda_ );

    matrix_dist specdata ( nlambda_ * nspec_, 1 );
    vector < bool > specsky;
    child_spec->read ( specdata, lambda_, specsky );

    fake_first_lambda_ = lambda_[0];
    fake_last_lambda_ = lambda_[ nlambda_ - 1 ];

    hdus_[ psf_gauss_hdu_x ] = 1;
    hdus_[ psf_gauss_hdu_y ] = 2;
    hdus_[ psf_gauss_hdu_lambda ] = 3;
    hdus_[ psf_gauss_hdu_amp ] = 4;
    hdus_[ psf_gauss_hdu_maj ] = 5;
    hdus_[ psf_gauss_hdu_min ] = 6;
    hdus_[ psf_gauss_hdu_ang ] = 7;
 
  } else {

    path_ = props.get < string > ( psf_gauss_key_path );

    pixcorr_ = props.get < int > ( psf_gauss_key_corr );

    rows_ = props.get < size_t > ( psf_gauss_key_rows );

    cols_ = props.get < size_t > ( psf_gauss_key_cols );

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

    lambda_.resize ( nlambda_ );
    fits::img_seek ( fp, hdus_[ psf_gauss_hdu_lambda ] );      
    fits::img_read ( fp, lambda_ );

    fits::close ( fp );

  }

  nglobal_ = nspec_ * nlambda_;

  npix_ = rows_ * cols_;

  // either read the data or generate it

  resp_.resize ( nglobal_ );

  if ( dofake_ ) {

    size_t xpix;
    size_t ypix;

    for ( size_t spec = 0; spec < nspec_; ++spec ) {
      for ( size_t lambda = 0; lambda < nlambda_; ++lambda ) {
        bin = spec * nlambda_ + lambda;

        fake_spec2pix ( spec, specbin, ypix, xpix );

        resp_[i].x = (double)xpix;
        resp_[i].y = (double)ypix;
        resp_[i].lambda = lambda_[ lambda ];
        resp_[i].amp = 1.0;
        resp_[i].maj = fake_psf_fwhm_ / 2.0;
        resp_[i].min = fake_psf_fwhm_ / 2.0;
        resp_[i].ang = 0.0;

      }
    }

  } else {

    vector_double buffer;

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

  }

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


harp::psf_gauss::~psf_gauss ( ) {
  
}


void harp::psf_gauss::fake_spec2pix ( size_t spec, size_t bin, size_t & row, size_t & col ) {
  size_t bundle = (size_t)( spec / fake_bundle_size_ );
  size_t bundle_spec = spec - bundle * fake_bundle_size_;
  col = (bundle * fake_pix_bundle_) + fake_pix_margin_ + ( bundle_spec * ( fake_pix_gap_ + 1 ) );
  row = bin + pixcorr_;
  return;
}


void harp::psf_gauss::extent ( size_t firstspec, size_t lastspec, size_t firstbin, size_t lastbin, size_t & firstcol, size_t & firstrow, size_t & lastcol, size_t & lastrow ) {
  
  double minX = 1000000000.0;
  double minY = 1000000000.0;
  double maxX = 0.0;
  double maxY = 0.0;

  if ( dofake_ ) {

    size_t firstcol, firstrow;
    size_t lastcol, lastrow;

    fake_spec2pix ( firstspec, firstbin, firstrow, firstcol );
    fake_spec2pix ( lastspec, lastbin, lastrow, lastcol );

    minY = (double)firstrow;
    maxY = (double)lastrow;
    minX = (double)firstcol;
    maxX = (double)lastcol;

  } else {
    
    for ( size_t curspec = firstspec; curspec <= lastspec; ++curspec ) {
      for ( size_t curbin = firstbin; curbin <= lastbin; ++curbin ) {
        if ( resp_[ curspec ]->x.Get( curbin, 0 ) < minX ) {
          minX = resp_[ curspec ]->x.Get( curbin, 0 );
        }
        if ( resp_[ curspec ]->x.Get( curbin, 0 ) > maxX ) {
          maxX = resp_[ curspec ]->x.Get( curbin, 0 );
        }
        if ( resp_[ curspec ]->y.Get( curbin, 0 ) < minY ) {
          minY = resp_[ curspec ]->y.Get( curbin, 0 );
        }
        if ( resp_[ curspec ]->y.Get( curbin, 0 ) > maxY ) {
          maxY = resp_[ curspec ]->y.Get( curbin, 0 );
        }
      }
    }

  }


  //cerr << "center for spec " << firstspec << ", bin " << firstbin << " to " << lastspec << ", " << lastbin << " is " << firstrow << "," << firstcol << " - " << lastrow << "," << lastcol << endl;
  
  if ( (int)minX - (int)pixcorr_ < 0 ) {
    firstcol = 0;
  } else if ( (int)minX - (int)pixcorr_ >= (int)cols_ ) {
    firstcol = cols_ - 1;
  } else {
    firstcol = minX - pixcorr_;
  }
  
  if ( (int)minY - (int)pixcorr_ < 0 ) {
    firstrow = 0;
  } else if ( (int)minY - (int)pixcorr_ >= (int)rows_ ) {
    firstrow = rows_ - 1;
  } else {
    firstrow = minY - pixcorr_;
  }
  
  if ( (int)maxX + (int)pixcorr_ >= (int)cols_ ) {
    lastcol = cols_ - 1;
  } else if ( (int)maxX + (int)pixcorr_ < 0 ) {
    lastcol = 0;
  } else {
    lastcol = maxX + pixcorr_;
  }

  if ( (int)maxY + (int)pixcorr_ >= (int)rows_ ) {
    lastrow = rows_ - 1;
  } else if ( (int)maxY + (int)pixcorr_ < 0 ) {
    lastrow = 0;
  } else {
    lastrow = maxY + pixcorr_;
  }

  //cerr << "extent for spec " << firstspec << ", bin " << firstbin << " to " << lastspec << ", " << lastbin << " is " << firstrow << "," << firstcol << " - " << lastrow << "," << lastcol << endl;

  return;
}





void harp::psf_gauss::projection ( size_t first_spec, size_t last_spec, size_t first_lambda, size_t last_lambda, matrix_sparse & AT ) {

  if ( ( first_spec >= nspec_ ) || ( last_spec >= nspec_ ) ) {
    HARP_THROW( "specs out of range" );
  }

  if ( ( first_lambda >= nlambda_ ) || ( last_lambda >= nlambda_ ) ) {
    HARP_THROW( "lambda points out of range" );
  }

  size_t block_nspec = last_spec - first_spec + 1;
  size_t block_nlambda = last_lambda - first_lambda + 1;
  size_t block_nbins = block_nspec * block_nlambda;
  
  AT.ResizeTo ( block_nbins, npix_ );

  // Get local matrix row range

  size_t first_loc_row = AT.FirstLocalRow();
  size_t loc_height = AT.LocalHeight();

  // cache data

  size_t first_spec_local = (size_t)( first_loc_row / block_nlambda );
  size_t last_spec_local = (size_t)( (first_loc_row + loc_height - 1) / block_nlambda );

  cache ( 0, nspec_ - 1 );

  // populate matrix

  AT.StartAssembly();

  size_t nnz = (2 * pixcorr_ + 1) * (2 * pixcorr_ + 1);

  AT.Reserve ( nnz * loc_height );

  size_t bin_firstrow;
  size_t bin_lastrow;
  size_t bin_firstcol;
  size_t bin_lastcol;

  double amp;
  double maj;
  double min;
  double ang;
  double xcenter;
  double ycenter;

  for ( size_t loc_bin = 0; loc_bin < loc_height; ++loc_bin ) {

    size_t block_bin = loc_bin + first_loc_row;
    size_t block_spec = (size_t)( block_bin / block_nlambda );
    size_t block_specbin = block_bin - ( block_spec * block_nlambda );

    size_t spec = first_spec + block_spec;
    size_t specbin = first_lambda + block_specbin;
    size_t bin = ( spec * nlambda_ ) + specbin;

    extent ( spec, spec, specbin, specbin, bin_firstcol, bin_firstrow, bin_lastcol, bin_lastrow );

    size_t bin_nrow = bin_lastrow - bin_firstrow + 1;
    size_t bin_ncol = bin_lastcol - bin_firstcol + 1;
    size_t bin_npix = bin_nrow * bin_ncol;

    if ( dofake_ ) {

      size_t xpix;
      size_t ypix;

      fake_spec2pix ( spec, specbin, ypix, xpix );

      xcenter = (double)xpix;
      ycenter = (double)ypix;

      amp = 1.0;
      maj = fake_psf_fwhm_ / 2.0;
      min = fake_psf_fwhm_ / 2.0;
      ang = 0.0;

    } else {

      amp = resp_[ spec ]->amp.Get( specbin, 0 );
      maj = resp_[ spec ]->maj.Get( specbin, 0 );
      min = resp_[ spec ]->min.Get( specbin, 0 );
      ang = resp_[ spec ]->ang.Get( specbin, 0 );
      xcenter = resp_[ spec ]->x.Get( specbin, 0 );
      ycenter = resp_[ spec ]->y.Get( specbin, 0 );

    }

    //cerr << "bin " << loc_bin << ": [" << xcenter << "," << ycenter << "] ( " << bin_firstcol << " - " << bin_lastcol << ", " << bin_firstrow << " - " << bin_lastrow << " )" << endl;

    // compute pixel distances from center

    matrix_local fxdist ( bin_npix, 1 );
    matrix_local fydist ( bin_npix, 1 );
  
    size_t p = 0;
    for ( size_t imgrow = bin_firstrow; imgrow <= bin_lastrow; ++imgrow ) {
    
      for ( size_t imgcol = bin_firstcol; imgcol <= bin_lastcol; ++imgcol ) {
    
        fxdist.Set( p, 0, (double)imgcol - xcenter );
        fydist.Set( p, 0, (double)imgrow - ycenter );

        ++p;
      }
    
    }
    
    matrix_local vals ( bin_npix, 1 );

    gauss_sample ( vals, fxdist, fydist, amp, maj, min, ang );

    p = 0;
    for ( size_t imgrow = bin_firstrow; imgrow <= bin_lastrow; ++imgrow ) {
    
      for ( size_t imgcol = bin_firstcol; imgcol <= bin_lastcol; ++imgcol ) {

        size_t pix = imgrow * cols_ + imgcol;

        AT.Update ( block_bin, pix, vals.Get( p, 0 ) );

        ++p;
      }

    }

  }

  AT.StopAssembly();
  
  return;
}






