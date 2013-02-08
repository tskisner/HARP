// @COPYRIGHT@

#include <harp_internal.hpp>

#include <boost/optional.hpp>

#include <cmath>

using namespace std;
using namespace harp;

static const char * format_sandbox = "sandbox";

static const char * sandbox_psf_key_path = "path";
static const char * sandbox_psf_key_corr = "corr";
static const char * sandbox_psf_key_rows = "imgrows";
static const char * sandbox_psf_key_cols = "imgcols";

static const char * sandbox_psf_key_name = "PSFPARAM";

static const char * sandbox_psf_hdu_x = "X";
static const char * sandbox_psf_hdu_y = "Y";
static const char * sandbox_psf_hdu_lambda = "Wavelength";
static const char * sandbox_psf_hdu_amp = "Amplitude";
static const char * sandbox_psf_hdu_maj = "MajorAxis";
static const char * sandbox_psf_hdu_min = "MinorAxis";
static const char * sandbox_psf_hdu_ang = "Angle";

static const char * sandbox_psf_key_fake = "FAKE";


harp::psf_sandbox::psf_sandbox ( boost::property_tree::ptree const & props ) : psf ( format_sandbox, props ) {

  boost::optional < string > fakeval = props.get_optional < string > ( sandbox_psf_key_fake );
  string fake = boost::get_optional_value_or ( fakeval, "FALSE" );
  dofake_ = ( fake == "TRUE" );

  if ( dofake_ ) {

    path_ = "";
    
    fake_n_bundle_ = 1;
    fake_bundle_size_ = 25;

    nspec_ = fake_n_bundle_ * fake_bundle_size_;
    specsize_ = 100;

    fake_pix_margin_ = 5;
    fake_pix_gap_ = 7;

    fake_pix_bundle_ = 2 * fake_pix_margin_ + (fake_bundle_size_ - 1) * fake_pix_gap_ + fake_bundle_size_;

    cols_ = fake_pix_bundle_ * fake_n_bundle_;
    rows_ = specsize_;

    // magnitude of continuum
    fake_background_ = 10.0;

    // amplitude of atmosphere "spikes"
    fake_peak_amp_ = 50.0 * fake_background_;

    // spacing of spikes
    fake_peak_space_ = 20;
    //fake_peak_space_ = (size_t)( specsize_ / 10 );

    // object size
    fake_peak_obj_ = 8.0 * fake_background_;

    // response fwhm
    fake_psf_fwhm_ = 1.5;

    // response correlation length in pixels
    pixcorr_ = (int)( 5.0 * (double)fake_psf_fwhm_ );

  } else {

    path_ = props.get < string > ( sandbox_psf_key_path );

    pixcorr_ = props.get < int > ( sandbox_psf_key_corr );

    rows_ = props.get < size_t > ( sandbox_psf_key_rows );

    cols_ = props.get < size_t > ( sandbox_psf_key_cols );

    int np;
    int myp;

    MPI_Comm_size ( MPI_COMM_WORLD, &np );
    MPI_Comm_rank ( MPI_COMM_WORLD, &myp );

    if ( myp == 0 ) {
  
      fitsfile *fp;

      fits::open_read ( fp, path_ );

      int hdu = fits::img_seek ( fp, sandbox_psf_key_name, sandbox_psf_hdu_x );
      hdus_[ sandbox_psf_hdu_x ] = hdu;
      fits::img_dims ( fp, nspec_, specsize_ );

      hdus_[ sandbox_psf_hdu_y ] = hdu_info ( fp, sandbox_psf_hdu_y );
      hdus_[ sandbox_psf_hdu_lambda ] = hdu_info ( fp, sandbox_psf_hdu_lambda );
      hdus_[ sandbox_psf_hdu_amp ] = hdu_info ( fp, sandbox_psf_hdu_amp );
      hdus_[ sandbox_psf_hdu_maj ] = hdu_info ( fp, sandbox_psf_hdu_maj );
      hdus_[ sandbox_psf_hdu_min ] = hdu_info ( fp, sandbox_psf_hdu_min );
      hdus_[ sandbox_psf_hdu_ang ] = hdu_info ( fp, sandbox_psf_hdu_ang );
      
      fits::close ( fp );

    }

    // broadcast 

    int ret = MPI_Bcast ( (void*)(&nspec_), 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD );
    mpi_check ( MPI_COMM_WORLD, ret );

    ret = MPI_Bcast ( (void*)(&specsize_), 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD );
    mpi_check ( MPI_COMM_WORLD, ret );

  }

  cached_ = false;
  cached_first_ = 0;
  cached_last_ = 0;

  nglobal_ = nspec_ * specsize_;

  npix_ = rows_ * cols_;
  
}


int harp::psf_sandbox::hdu_info ( fitsfile *fp, const char * sandbox_psf_hdu ) {
  int hdu = fits::img_seek ( fp, sandbox_psf_key_name, sandbox_psf_hdu );
  size_t rows, cols;
  fits::img_dims ( fp, rows, cols );
  if ( ( rows != nspec_ ) || ( cols != specsize_ ) ) {
    HARP_THROW( "sandbox_psf: PSF file must have identical dimensions for all HDUs" );
  }
  return hdu;
}


harp::psf_sandbox::~psf_sandbox ( ) {
  
  cleanup();
  
}


boost::property_tree::ptree harp::psf_sandbox::serialize ( ) {
  boost::property_tree::ptree ret;

  ret.put ( "format", psf::format() );

  if ( dofake_ ) {

    ret.put ( sandbox_psf_key_fake, "TRUE" );

  } else {

    ret.put ( sandbox_psf_key_path, path_ );

    ret.put ( sandbox_psf_key_corr, pixcorr_ );

  }

  return ret;
}


void harp::psf_sandbox::cache ( size_t firstlocal, size_t nlocal ) {

  size_t firstspec = firstlocal % specsize_;
  size_t lastspec = ( firstlocal + nlocal - 1 ) % specsize_;
  size_t nspec = lastspec - firstspec + 1;

  if ( ! cached_ ) {

    cached_ = true;
    cached_first_ = firstspec;
    cached_last_ = lastspec;

    if ( dofake_ ) {
      return;
    }

    int np;
    int myp;
    int ret;

    MPI_Comm_size ( MPI_COMM_WORLD, &np );
    MPI_Comm_rank ( MPI_COMM_WORLD, &myp );
    
    fitsfile *fp;

    if ( myp == 0 ) {
      fits::open_read ( fp, path_ );
    }

    matrix_local iobuffer ( nspec_, specsize_ );

    if ( myp == 0 ) {
      fits::img_seek ( fp, hdus_[ sandbox_psf_hdu_x ] );      
      fits::img_read ( fp, iobuffer );
    }
    ret = MPI_Bcast ( (void*)(iobuffer.Buffer()), nglobal_, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    mpi_check ( MPI_COMM_WORLD, ret );
    for ( size_t i = 0; i < nspec; ++i ) {
      size_t spec = i + firstspec;
      resp_[ spec ].x.ResizeTo ( specsize_, 1 );
      for ( size_t j = 0; j < specsize_; ++j ) {
        resp_[ spec ].x.Set ( j, 0, iobuffer.Get ( spec, j ) );
      }
    }

    if ( myp == 0 ) {
      fits::img_seek ( fp, hdus_[ sandbox_psf_hdu_y ] );      
      fits::img_read ( fp, iobuffer );
    }
    ret = MPI_Bcast ( (void*)(iobuffer.Buffer()), nglobal_, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    mpi_check ( MPI_COMM_WORLD, ret );
    for ( size_t i = 0; i < nspec; ++i ) {
      size_t spec = i + firstspec;
      resp_[ spec ].y.ResizeTo ( specsize_, 1 );
      for ( size_t j = 0; j < specsize_; ++j ) {
        resp_[ spec ].y.Set ( j, 0, iobuffer.Get ( spec, j ) );
      }
    }

    if ( myp == 0 ) {
      fits::img_seek ( fp, hdus_[ sandbox_psf_hdu_lambda ] );      
      fits::img_read ( fp, iobuffer );
    }
    ret = MPI_Bcast ( (void*)(iobuffer.Buffer()), nglobal_, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    mpi_check ( MPI_COMM_WORLD, ret );
    for ( size_t i = 0; i < nspec; ++i ) {
      size_t spec = i + firstspec;
      resp_[ spec ].lambda.ResizeTo ( specsize_, 1 );
      for ( size_t j = 0; j < specsize_; ++j ) {
        resp_[ spec ].lambda.Set ( j, 0, iobuffer.Get ( spec, j ) );
      }
    }

    if ( myp == 0 ) {
      fits::img_seek ( fp, hdus_[ sandbox_psf_hdu_amp ] );      
      fits::img_read ( fp, iobuffer );
    }
    ret = MPI_Bcast ( (void*)(iobuffer.Buffer()), nglobal_, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    mpi_check ( MPI_COMM_WORLD, ret );
    for ( size_t i = 0; i < nspec; ++i ) {
      size_t spec = i + firstspec;
      resp_[ spec ].amp.ResizeTo ( specsize_, 1 );
      for ( size_t j = 0; j < specsize_; ++j ) {
        resp_[ spec ].amp.Set ( j, 0, iobuffer.Get ( spec, j ) );
      }
    }

    if ( myp == 0 ) {
      fits::img_seek ( fp, hdus_[ sandbox_psf_hdu_maj ] );      
      fits::img_read ( fp, iobuffer );
    }
    ret = MPI_Bcast ( (void*)(iobuffer.Buffer()), nglobal_, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    mpi_check ( MPI_COMM_WORLD, ret );
    for ( size_t i = 0; i < nspec; ++i ) {
      size_t spec = i + firstspec;
      resp_[ spec ].maj.ResizeTo ( specsize_, 1 );
      for ( size_t j = 0; j < specsize_; ++j ) {
        resp_[ spec ].maj.Set ( j, 0, iobuffer.Get ( spec, j ) );
      }
    }

    if ( myp == 0 ) {
      fits::img_seek ( fp, hdus_[ sandbox_psf_hdu_min ] );      
      fits::img_read ( fp, iobuffer );
    }
    ret = MPI_Bcast ( (void*)(iobuffer.Buffer()), nglobal_, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    mpi_check ( MPI_COMM_WORLD, ret );
    for ( size_t i = 0; i < nspec; ++i ) {
      size_t spec = i + firstspec;
      resp_[ spec ].min.ResizeTo ( specsize_, 1 );
      for ( size_t j = 0; j < specsize_; ++j ) {
        resp_[ spec ].min.Set ( j, 0, iobuffer.Get ( spec, j ) );
      }
    }

    if ( myp == 0 ) {
      fits::img_seek ( fp, hdus_[ sandbox_psf_hdu_ang ] );      
      fits::img_read ( fp, iobuffer );
    }
    ret = MPI_Bcast ( (void*)(iobuffer.Buffer()), nglobal_, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    mpi_check ( MPI_COMM_WORLD, ret );
    for ( size_t i = 0; i < nspec; ++i ) {
      size_t spec = i + firstspec;
      resp_[ spec ].ang.ResizeTo ( specsize_, 1 );
      for ( size_t j = 0; j < specsize_; ++j ) {
        resp_[ spec ].ang.Set ( j, 0, iobuffer.Get ( spec, j ) );
      }
    }
    
    if ( myp == 0 ) {
      fits::close ( fp );
    }

  } else {
    if ( ( cached_first_ != firstspec ) || ( cached_last_ != lastspec ) ) {
      std::ostringstream o;
      o << "attempting to change which spectra are cached locally (" << cached_first_ << "-" << cached_last_ << ") --> (" << firstspec << "-" << lastspec << ")";
      HARP_THROW( o.str().c_str() );
    }
  }
  
  return;
}


void harp::psf_sandbox::fake_spec2pix ( size_t spec, size_t bin, size_t & row, size_t & col ) {
  size_t bundle = (size_t)( spec / fake_bundle_size_ );
  size_t bundle_spec = spec - bundle * fake_bundle_size_;
  col = (bundle * fake_pix_bundle_) + fake_pix_margin_ + ( bundle_spec * ( fake_pix_gap_ + 1 ) );
  row = bin;
  return;
}


void harp::psf_sandbox::extent ( size_t firstspec, size_t lastspec, size_t firstbin, size_t lastbin, size_t & firstcol, size_t & firstrow, size_t & lastcol, size_t & lastrow ) {
  
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
        if ( resp_[ curspec ].x.Get( curbin, 0 ) < minX ) {
          minX = resp_[ curspec ].x.Get( curbin, 0 );
        }
        if ( resp_[ curspec ].x.Get( curbin, 0 ) > maxX ) {
          maxX = resp_[ curspec ].x.Get( curbin, 0 );
        }
        if ( resp_[ curspec ].y.Get( curbin, 0 ) < minY ) {
          minY = resp_[ curspec ].y.Get( curbin, 0 );
        }
        if ( resp_[ curspec ].y.Get( curbin, 0 ) > maxY ) {
          maxY = resp_[ curspec ].y.Get( curbin, 0 );
        }
      }
    }

  }


  //cerr << "center for spec " << firstspec << ", bin " << firstbin << " to " << lastspec << ", " << lastbin << " is " << firstrow << "," << firstcol << " - " << lastrow << "," << lastcol << endl;
  
  if ( (int)minX - (int)pixcorr_ < 0 ) {
    firstcol = 0;
  } else {
    firstcol = minX - pixcorr_;
  }
  
  if ( (int)minY - (int)pixcorr_ < 0 ) {
    firstrow = 0;
  } else {
    firstrow = minY - pixcorr_;
  }
  
  if ( maxX + pixcorr_ >= cols_ ) {
    lastcol = cols_ - 1;
  } else {
    lastcol = maxX + pixcorr_;
  }

  if ( maxY + pixcorr_ >= rows_ ) {
    lastrow = rows_ - 1;
  } else {
    lastrow = maxY + pixcorr_;
  }

  //cerr << "extent for spec " << firstspec << ", bin " << firstbin << " to " << lastspec << ", " << lastbin << " is " << firstrow << "," << firstcol << " - " << lastrow << "," << lastcol << endl;

  return;
}


void harp::psf_sandbox::gauss_sample ( matrix_local & vals, matrix_local & xrel, matrix_local & yrel, double amp, double maj, double min, double ang ) {

  double PI = std::atan2 ( 0.0, -1.0 );
  
  amp /= maj * min * 2.0 * PI;
  
  size_t nvals = xrel.Height();
  
  double cang = cos ( ang );
  double sang = sin ( ang );
  
  double invmaj = 1.0 / maj;
  double invmin = 1.0 / min;
  
  double xt, yt, exparg;
  
  size_t i;
  
  for ( i = 0; i < nvals; ++i ) {
    xt = xrel.Get( i, 0 ) * cang + yrel.Get( i, 0 ) * sang;
    yt = - xrel.Get( i, 0 ) * sang + yrel.Get( i, 0 ) * cang;
    exparg = - 0.5 * ( xt * xt * invmaj * invmaj + yt * yt * invmin * invmin );
    vals.Set( i, 0, amp * exp ( exparg ) );
  }
  
  return;
}


void harp::psf_sandbox::projection ( size_t firstbin, size_t lastbin, matrix_sparse & AT ) {

  if ( ( firstbin >= specsize_ ) || ( lastbin >= specsize_ ) ) {
    HARP_THROW( "spectral bins out of range" );
  }

  size_t nflux = nspec_ * ( lastbin - firstbin + 1 );
  
  AT.ResizeTo ( nflux, npix_ );

  // Get local matrix row range

  size_t first_loc_row = AT.FirstLocalRow();
  size_t loc_height = AT.LocalHeight();

  // cache data

  cache ( first_loc_row, loc_height );

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

    size_t glob_bin = loc_bin + first_loc_row;

    size_t spec = (size_t) ( glob_bin / specsize_ );
    size_t specbin = glob_bin - ( spec * specsize_ );

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
      maj = fake_psf_fwhm_;
      min = fake_psf_fwhm_;
      ang = 0.0;

    } else {

      amp = resp_[ spec ].amp.Get( specbin, 0 );
      maj = resp_[ spec ].maj.Get( specbin, 0 );
      min = resp_[ spec ].min.Get( specbin, 0 );
      ang = resp_[ spec ].ang.Get( specbin, 0 );
      xcenter = resp_[ spec ].x.Get( specbin, 0 );
      ycenter = resp_[ spec ].y.Get( specbin, 0 );

    }

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

        AT.Update ( glob_bin, pix, vals.Get( p, 0 ) );

        ++p;
      }

    }

  }

  AT.StopAssembly();
  
  return;
}


void harp::psf_sandbox::fake_spec ( matrix_dist & data ) {

  size_t nbins = nspec_ * specsize_;

  data.ResizeTo ( nbins, 1 );

  double PI = std::atan2 ( 0.0, -1.0 );

  size_t half_peak = fake_peak_space_ >> 1;
  size_t half_obj = (size_t)(1.5 * (double)half_peak);
  size_t obj_space = 2 * half_obj;

  size_t bin = 0;

  for ( size_t i = 0; i < nspec_; ++i ) {
    for ( size_t j = 0; j < specsize_; ++j ) {

      double val = fake_background_ * sin ( 3.0 * (double)j / ( fake_peak_space_ * 2.0 * PI ) ) * sin ( 7.0 * (double)j / ( fake_peak_space_ * 2.0 * PI ) ) * sin ( 11.0 * (double)j / ( fake_peak_space_ * 2.0 * PI ) );

      val += 2.0 * fake_background_;

      if ( ( ( j + half_peak ) % fake_peak_space_ ) == 0 ) {
        val += fake_peak_amp_;
      }

      size_t temp = i % 25;
      temp += j + half_obj;

      //cerr << " temp = " << temp << " obj_space = " << obj_space << endl;

      if ( ( temp % obj_space ) == 0 ) {
        val += fake_peak_obj_;
      }

      data.Set ( bin, 0, val );

      ++bin;
    }
  }

  return;
}




