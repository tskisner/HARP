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


harp::psf_toy::psf_toy ( std::map < std::string, std::string > const & params ) : psf ( format_toy, params ) {
  
  map < std::string, std::string > :: const_iterator val;
  
  val = params.find( toy_psf_key_path );
  
  if ( val == params.end() ) {
    MOAT_THROW( "toy_psf: must specify path to PSF file" );
  }
    
  path_ = val->second;
  
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
  
  // FIXME: this should be configurable at runtime
  pixcorr_ = (size_t) ( 5.0 * 4000.0 / (double) nspec_ ); 
  
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
      if ( ! fp ) {
        fits::open_read ( fp, path_ );
      }
      
      psf_toy_resp temp;
      resp_[ spec ] = temp;
      
      resp_[ spec ].x.resize ( nbins_ );
      resp_[ spec ].y.resize ( nbins_ );
      resp_[ spec ].lambda.resize ( nbins_ );
      resp_[ spec ].amp.resize ( nbins_ );
      resp_[ spec ].maj.resize ( nbins_ );
      resp_[ spec ].min.resize ( nbins_ );
      resp_[ spec ].ang.resize ( nbins_ );
      
      fits::img_seek ( fp, hdus_[ toy_psf_hdu_x ] );      
      fits::img_read_row_int ( fp, 1 + spec, resp_[ spec ].x );
      
      fits::img_seek ( fp, hdus_[ toy_psf_hdu_y ] );      
      fits::img_read_row_int ( fp, 1 + spec, resp_[ spec ].y );
      
      fits::img_seek ( fp, hdus_[ toy_psf_hdu_lambda ] );      
      fits::img_read_row ( fp, 1 + spec, resp_[ spec ].lambda );
      
      fits::img_seek ( fp, hdus_[ toy_psf_hdu_amp ] );      
      fits::img_read_row ( fp, 1 + spec, resp_[ spec ].amp );
      
      fits::img_seek ( fp, hdus_[ toy_psf_hdu_maj ] );      
      fits::img_read_row ( fp, 1 + spec, resp_[ spec ].maj );
      
      fits::img_seek ( fp, hdus_[ toy_psf_hdu_min ] );      
      fits::img_read_row ( fp, 1 + spec, resp_[ spec ].min );
      
      fits::img_seek ( fp, hdus_[ toy_psf_hdu_ang ] );      
      fits::img_read_row ( fp, 1 + spec, resp_[ spec ].ang );
      
    }
  }
  
  if ( fp ) {
    fits::close ( fp );
  }
  
  return;
}


void harp::psf_toy::lambda ( size_t specnum, data_vec & data ) {
  cerr << "lambda caching spectrum " << specnum << endl;
  cache_spec ( specnum, specnum );
  
  size_t bins = resp_[ specnum ].x.size();
  cerr << "lambda found " << bins << " spectral bins" << endl;
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


void harp::psf_toy::projection ( size_t firstX, size_t firstY, size_t lastX, size_t lastY, sparse_mat_view & data ) {
  
  size_t nx = lastX - firstX + 1;
  size_t ny = lastY - firstY + 1;
  
  size_t nbins = data.size1();
  
  size_t npix = data.size2();
  
  if ( npix != nx * ny ) {
    MOAT_THROW( "toy_psf: PSF projection ranges must match dimensions of projection data" );
  }
  
  
  
  for ( size_t imgcol = firstX; imgcol <= lastX; ++imgcol ) {
    
    for ( size_t imgrow = firstY; imgrow <= lastY; ++imgrow ) {
      
      
      
      
    }
    
  }
  
  
  return;
}




