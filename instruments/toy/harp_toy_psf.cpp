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


double harp::psf_toy::gauss_sample ( double xrel, double yrel, double amp, double maj, double min, double ang ) {
  double val;
  
  double cxx;
  double cxy;
  double cyy;
  double cang;
  double sang;
  double s2ang;
  
  cang = cos ( ang );
  sang = sin ( ang );
  s2ang = sin ( 2.0 * ang );

  cxx = ( cang * cang ) / ( 2.0 * maj * maj ) + ( sang * sang ) / ( 2.0 * min * min );

  cxy = - s2ang / ( 4.0 * maj * maj ) + s2ang / ( 4.0 * min * min );

  cyy = ( sang * sang ) / ( 2.0 * maj * maj ) + ( cang * cang ) / ( 2.0 * min * min );

  val = amp * exp ( - ( cxx * xrel * xrel + 2.0 * cxy * xrel * yrel + cyy * yrel * yrel ) );
  
  //if ( isnan ( val ) ) {
  //  cout << "xrel = " << xrel << " yrel = " << yrel << " cxx = " << cxx << " cxy = " << cxy << " cyy = " << cyy << endl;
  //  cout << "    amp = " << amp << " maj = " << maj << " min = " << min << " ang = " << ang << endl;
  //}
  
  return val;
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
  
  /*
  for ( size_t i = 0; i < npix; ++i ) {
    for ( size_t j = 0; j < nbins; ++j ) {
      data ( i, j ) = 0.0;
    }
  }
  */
  
  std::map < size_t, psf_toy_resp > :: iterator itspec;
  
  int64_t xdist;
  double fxdist;
  int64_t ydist;
  double fydist;
  
  double val;
  size_t datarow;
  size_t datacol = 0;
  
  for ( itspec = resp_.begin(); itspec != resp_.end(); ++itspec ) {
    
    size_t specindx = itspec->first;
    
    if ( ( specindx >= firstspec ) && ( specindx <= lastspec ) ) {
      
      // process this spectrum
      
      psf_toy_resp & spec = itspec->second;
      size_t specbins = spec.x.size();
    
      for ( size_t bin = 0; bin < specbins; ++bin ) {
      
        if ( ( bin >= firstbin ) && ( bin <= lastbin ) ) {
          
          // process this flux bin
          
          //cerr << "processing bin " << bin << " of spectrum " << specindx << endl;
    
          for ( size_t imgcol = firstX; imgcol <= lastX; ++imgcol ) {
    
            for ( size_t imgrow = firstY; imgrow <= lastY; ++imgrow ) {
              
              xdist = spec.x[bin] - (int)imgcol;
              ydist = spec.y[bin] - (int)imgrow;
              
              if ( ( labs ( (long)xdist ) < pixcorr_ ) && ( labs ( (long)ydist ) < pixcorr_ ) ) {
              
                fxdist = (double)xdist;
                fydist = (double)ydist;
              
                val = gauss_sample ( fxdist, fydist, spec.amp[bin], spec.maj[bin], spec.min[bin], spec.ang[bin] );
                
                datarow = ( imgrow - firstY ) * ( lastX - firstX + 1 ) + ( imgcol - firstX );
                
                //if ( isnan ( val ) ) {
                //  cout << "    NAN at spec " << specindx << ", bin " << bin << ": imgcol = " << imgcol << " imgrow = " << imgrow << endl;
                //}
              
                //cout << "data(" << datarow << "," << datacol << ") = " << data(datarow,datacol);
                data ( datarow, datacol ) += val;
                //cout << " + " << val << " ==> " << data(datarow, datacol) << endl;
              
              }
      
            }
      
          }
          
          ++datacol;
    
        }
      }
    }
  }
  
  
  return;
}




