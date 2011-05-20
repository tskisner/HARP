// @COPYRIGHT@

#include <harp_internal.hpp>


using namespace std;
using namespace harp;

static const char * format_toy = "toy";

static const char * toy_psf_key_path = "path";

static const char * toy_psf_key_name = "PSFPARAM";

static const char * toy_psf_hdu_x = "X";
static const char * toy_psf_hdu_y = "Y";
static const char * toy_psf_hdu_lambda = "Wavelength";
static const char * toy_psf_hdu_amp = "Amplitude";
static const char * toy_psf_hdu_maj = "MajorAxis";
static const char * toy_psf_hdu_min = "MinorAxis";
static const char * toy_psf_hdu_ang = "Angle";

static const char * toy_psf_key_corr = "corr";


harp::psf_toy::psf_toy ( std::map < std::string, std::string > const & params ) : psf ( format_toy, params ) {
  
  map < std::string, std::string > :: const_iterator val;
  
  val = params.find( toy_psf_key_path );
  
  if ( val == params.end() ) {
    MOAT_THROW( "toy_psf: must specify path to PSF file" );
  }
    
  path_ = val->second;
  
  val = params.find( toy_psf_key_corr );
  
  if ( val == params.end() ) {
    MOAT_THROW( "toy_psf: must specify pixel space correlation length" );
  }
  
  pixcorr_ = atoi ( val->second.c_str() );
  
  
  fitsfile *fp;

  fits::open_read ( fp, path_ );
  
  int hdu = fits::img_seek ( fp, toy_psf_key_name, toy_psf_hdu_x );
  hdus_[ toy_psf_hdu_x ] = hdu;
  fits::img_dims ( fp, nspec_, nbins_ );
  
  size_t rows, cols;
  
  hdu = fits::img_seek ( fp, toy_psf_key_name, toy_psf_hdu_y );
  hdus_[ toy_psf_hdu_y ] = hdu;
  fits::img_dims ( fp, rows, cols );
  if ( ( rows != nspec_ ) || ( cols != nbins_ ) ) {
    MOAT_THROW( "toy_psf: PSF file must have identical dimensions for all HDUs" );
  }
  
  hdu = fits::img_seek ( fp, toy_psf_key_name, toy_psf_hdu_lambda );
  hdus_[ toy_psf_hdu_lambda ] = hdu;
  fits::img_dims ( fp, rows, cols );
  if ( ( rows != nspec_ ) || ( cols != nbins_ ) ) {
    MOAT_THROW( "toy_psf: PSF file must have identical dimensions for all HDUs" );
  }
  
  hdu = fits::img_seek ( fp, toy_psf_key_name, toy_psf_hdu_amp );
  hdus_[ toy_psf_hdu_amp ] = hdu;
  fits::img_dims ( fp, rows, cols );
  if ( ( rows != nspec_ ) || ( cols != nbins_ ) ) {
    MOAT_THROW( "toy_psf: PSF file must have identical dimensions for all HDUs" );
  }
  
  hdu = fits::img_seek ( fp, toy_psf_key_name, toy_psf_hdu_maj );
  hdus_[ toy_psf_hdu_maj ] = hdu;
  fits::img_dims ( fp, rows, cols );
  if ( ( rows != nspec_ ) || ( cols != nbins_ ) ) {
    MOAT_THROW( "toy_psf: PSF file must have identical dimensions for all HDUs" );
  }
  
  hdu = fits::img_seek ( fp, toy_psf_key_name, toy_psf_hdu_min );
  hdus_[ toy_psf_hdu_min ] = hdu;
  fits::img_dims ( fp, rows, cols );
  if ( ( rows != nspec_ ) || ( cols != nbins_ ) ) {
    MOAT_THROW( "toy_psf: PSF file must have identical dimensions for all HDUs" );
  }
  
  hdu = fits::img_seek ( fp, toy_psf_key_name, toy_psf_hdu_ang );
  hdus_[ toy_psf_hdu_ang ] = hdu;
  fits::img_dims ( fp, rows, cols );
  if ( ( rows != nspec_ ) || ( cols != nbins_ ) ) {
    MOAT_THROW( "toy_psf: PSF file must have identical dimensions for all HDUs" );
  }
  
  fits::close ( fp );
  
}


harp::psf_toy::~psf_toy ( ) {
  
  cleanup();
  
}


void harp::psf_toy::cache_spec ( size_t first, size_t last ) {
  
  fitsfile *fp = NULL;

  for ( size_t spec = first; spec <= last; ++spec ) {
    map < size_t, psf_toy_resp > :: iterator check;
    
    check = resp_.find ( spec );

    if ( check == resp_.end() ) {
      //cerr << "cache_spec reading spectrum " << spec << endl;
      if ( ! fp ) {
        //cerr << "cache_spec opening file " << path_ << endl;
        fits::open_read ( fp, path_ );
      }
      
      psf_toy_resp temp;
      resp_[ spec ] = temp;
      
      //cerr << "cache_spec resizing vectors to " << nbins_ << " elements" << endl;
      
      resp_[ spec ].x.resize ( nbins_ );
      resp_[ spec ].y.resize ( nbins_ );
      resp_[ spec ].lambda.resize ( nbins_ );
      resp_[ spec ].amp.resize ( nbins_ );
      resp_[ spec ].maj.resize ( nbins_ );
      resp_[ spec ].min.resize ( nbins_ );
      resp_[ spec ].ang.resize ( nbins_ );
      
      fits::img_seek ( fp, hdus_[ toy_psf_hdu_x ] );      
      fits::img_read_row_int ( fp, spec, resp_[ spec ].x );
      
      fits::img_seek ( fp, hdus_[ toy_psf_hdu_y ] );      
      fits::img_read_row_int ( fp, spec, resp_[ spec ].y );
      
      fits::img_seek ( fp, hdus_[ toy_psf_hdu_lambda ] );      
      fits::img_read_row ( fp, spec, resp_[ spec ].lambda );
      
      fits::img_seek ( fp, hdus_[ toy_psf_hdu_amp ] );      
      fits::img_read_row ( fp, spec, resp_[ spec ].amp );
      
      fits::img_seek ( fp, hdus_[ toy_psf_hdu_maj ] );      
      fits::img_read_row ( fp, spec, resp_[ spec ].maj );
      
      fits::img_seek ( fp, hdus_[ toy_psf_hdu_min ] );      
      fits::img_read_row ( fp, spec, resp_[ spec ].min );
      
      fits::img_seek ( fp, hdus_[ toy_psf_hdu_ang ] );      
      fits::img_read_row ( fp, spec, resp_[ spec ].ang );
      
      /*
      for ( size_t i = 0; i < nbins_; ++i ) {
        cout << "spec[" << spec << "](" << i << ") = " << resp_[ spec ].x[i] << " " << resp_[ spec ].y[i] << " " << resp_[ spec ].lambda[i] << " " << resp_[ spec ].amp[i] << " " << resp_[ spec ].maj[i] << " " << resp_[ spec ].min[i] << " " << resp_[ spec ].ang[i] << endl;
      }
      */
      
    }
  }
  
  if ( fp ) {
    fits::close ( fp );
  }
  
  return;
}


void harp::psf_toy::lambda ( size_t specnum, data_vec & data ) {
  //cerr << "lambda caching spectrum " << specnum << endl;
  cache_spec ( specnum, specnum );
  
  size_t bins = resp_[ specnum ].x.size();
  //cerr << "lambda found " << bins << " spectral bins" << endl;
  data.resize ( bins );
  
  for ( size_t i = 0; i < bins; ++i ) {
    data[ i ] = resp_[ specnum ].lambda[ i ];
  }
  
  return;
}


void harp::psf_toy::extent ( size_t firstspec, size_t lastspec, size_t firstbin, size_t lastbin, size_t & firstX, size_t & firstY, size_t & lastX, size_t & lastY ) {
  
  cache_spec ( firstspec, lastspec );
  
  int upleftX = resp_[ firstspec ].x[ firstbin ];
  int upleftY = resp_[ firstspec ].y[ firstbin ];
  
  int lowrightX = resp_[ lastspec ].x[ lastbin ];
  int lowrightY = resp_[ lastspec ].y[ lastbin ];
  
  upleftX -= (int)pixcorr_;
  upleftY += (int)pixcorr_;
  
  lowrightX += (int)pixcorr_;
  lowrightY -= (int)pixcorr_;
  
  if ( upleftX < 0 ) {
    upleftX = 0;
  }
  
  if ( lowrightY < 0 ) {
    lowrightY = 0;
  }
  
  firstX = upleftX;
  lastX = lowrightX;
  
  firstY = lowrightY;
  lastY = upleftY;
  
  return;
}


void harp::psf_toy::gauss_sample ( data_vec & vals, data_vec & xrel, data_vec & yrel, double amp, double maj, double min, double ang ) {
  
  //cerr << "gauss:  amp = " << amp << " maj = " << maj << " min = " << min << " ang = " << ang << endl;
  
  size_t nvals = xrel.size();
  
  double * buf = moat::double_alloc ( nvals );
  
  double cang = cos ( ang );
  double sang = sin ( ang );
  double s2ang = sin ( 2.0 * ang );
  
  double mjmj = 0.5 / ( maj * maj );
  double mnmn = 0.5 / ( min * min );
  
  double cxx = mjmj * cang * cang + mnmn * sang * sang;
  double cxy = 0.5 * mjmj * s2ang - 0.5 * mnmn * s2ang;
  double cyy  = mjmj * sang * sang + mnmn * cang * cang;

  for ( size_t i = 0; i < nvals; ++i ) {
    buf[i] = - ( cxx * xrel[i] * xrel[i] + 2.0 * cxy * xrel[i] * yrel[i] + cyy * yrel[i] * yrel[i] );
  }
  
  moat::sf::fast_exp ( nvals, buf, &(vals[0]) );
  
  for ( size_t i = 0; i < nvals; ++i ) {
    vals[i] *= amp;
  }

  return;
}


void harp::psf_toy::projection ( size_t firstspec, size_t lastspec, size_t firstbin, size_t lastbin, size_t firstX, size_t lastX, size_t firstY, size_t lastY, sparse_mat_view & data ) {
  
  size_t nx = lastX - firstX + 1;
  size_t ny = lastY - firstY + 1;
  
  size_t nbins = data.size2();
  
  size_t npix = data.size1();
  
  if ( ( npix != nx * ny ) || ( nbins != ( lastspec - firstspec + 1 ) * ( lastbin - firstbin + 1 ) ) ) {
    std::ostringstream o;
    o << "toy_psf: PSF projection ranges must match dimensions of projection data (" << npix << " x " << nbins << ")";
    MOAT_THROW( o.str().c_str() );
  }
  
  cache_spec ( firstspec, lastspec );
  
  vector < size_t > binlist ( nbins );
  
  size_t b = 0;
  for ( size_t i = firstspec; i <= lastspec; ++i ) {
    for ( size_t j = firstbin; j <= lastbin; ++j ) {
      binlist[b] = i * nbins_ + j;
      //cerr << "binlist[" << b << "] = " << binlist[b] << endl;
      ++b;
    }
  }

  #ifdef _OPENMP
  #pragma omp parallel for default(none) private(b) shared(nbins, binlist, firstX, firstY, lastX, lastY, data, cerr) schedule(static)
  #endif
  for ( b = 0; b < nbins; ++b ) {
    
    //cerr << "bin " << b << " of " << nbins << endl;
    
    size_t specbin = binlist[b] % nbins_;
    
    //cerr << "bin " << b << ": global = " << binlist[b] << " specbin = " << specbin << endl;
    
    psf_toy_resp & spec = resp_[ (size_t)(binlist[b] / nbins_) ];
    
    //cerr << "spec.x = " << spec.x[specbin] << " spec.y = " << spec.y[specbin] << endl;
    
    //cerr << "set spec" << endl;
    
    size_t startX;
    size_t stopX;
    size_t startY;
    size_t stopY;
    
    bool valid = true;
    
    if ( spec.x[specbin] > pixcorr_ ) {
      startX = spec.x[specbin] - pixcorr_;    
    } else {
      startX = 0;
    }
    if ( startX < firstX ) {
      startX = firstX;
    } else if ( startX > lastX ) {
      valid = false;
    }
    
    stopX = spec.x[specbin] + pixcorr_;
    if ( stopX > lastX ) {
      stopX = lastX;
    } else if ( stopX < firstX ) {
      valid = false;
    }
    
    if ( spec.y[specbin] > pixcorr_ ) {
      startY = spec.y[specbin] - pixcorr_;    
    } else {
      startY = 0;
    }
    if ( startY < firstY ) {
      startY = firstY;
    } else if ( startY > lastY ) {
      valid = false;
    }
    
    stopY = spec.y[specbin] + pixcorr_;
    if ( stopY > lastY ) {
      stopY = lastY;
    } else if ( stopY < firstY ) {
      valid = false;
    }
      
    //cerr << "pix range = [" << startX << ", " << stopX << "] and [" << startY << ", " << stopY << "]" << endl;
    
    if ( valid ) {
      
      size_t nvalid = ( stopY - startY + 1 ) * ( stopX - startX + 1);
      
      double amp = spec.amp[specbin];
      double maj = spec.maj[specbin];
      double min = spec.min[specbin];
      double ang = spec.ang[specbin];
      int xbin = spec.x[specbin];
      int ybin = spec.y[specbin];
      
      data_vec fxdist ( nvalid );
      data_vec fydist ( nvalid );
      
      int xdist;
      int ydist;
      
      size_t pix = 0;
      
      //cerr << "  computing distances" << endl;
      
      for ( size_t imgrow = startY; imgrow <= stopY; ++imgrow ) {
        
        for ( size_t imgcol = startX; imgcol <= stopX; ++imgcol ) {
        
          xdist = (int)imgcol - xbin;
          ydist = (int)imgrow - ybin;
        
          fxdist[pix] = (double)xdist;
          fydist[pix] = (double)ydist;
          
          ++pix;
        }
        
      }
      
      data_vec vals ( nvalid );
          
      //cerr << "  calling gauss_sample" << endl;
      
      gauss_sample ( vals, fxdist, fydist, amp, maj, min, ang );
      
      size_t datarow;
      
      pix = 0;
      
      #ifdef _OPENMP
      #pragma omp critical
      #endif
      {
        //cerr << "  incrementing matrix values" << endl;
        
        size_t rowoff;
        
        for ( size_t imgrow = startY; imgrow <= stopY; ++imgrow ) {
          rowoff = ( imgrow - firstY ) * ( lastX - firstX + 1 );
        
          for ( size_t imgcol = startX; imgcol <= stopX; ++imgcol ) {
          
            datarow = rowoff + ( imgcol - firstX );

            data ( datarow, b ) = vals[pix];
          
            ++pix;
          }
        }
      }
      
    }
    
  }
  
  
  return;
}




