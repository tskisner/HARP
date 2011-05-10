// @COPYRIGHT@

#include <harp_internal.hpp>


using namespace std;
using namespace harp;


static const char * toy_psf_key_path = "path";

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
  
  int hdu = fits::img_seek ( fp, toy_psf_hdu_x );
  hdus_[ toy_psf_hdu_x ] = hdu;
  fits::img_dims ( fp, nspec_, nbins_ );
  
  size_t rows, cols;
  
  hdu = fits::img_seek ( fp, toy_psf_hdu_y );
  hdus_[ toy_psf_hdu_y ] = hdu;
  fits::img_dims ( fp, rows, cols );
  if ( ( rows != nspec_ ) || ( cols != nbins_ ) ) {
    MOAT_THROW( "toy_psf: PSF file must have identical dimensions for all HDUs" );
  }
  
  hdu = fits::img_seek ( fp, toy_psf_hdu_lambda );
  hdus_[ toy_psf_hdu_lambda ] = hdu;
  fits::img_dims ( fp, rows, cols );
  if ( ( rows != nspec_ ) || ( cols != nbins_ ) ) {
    MOAT_THROW( "toy_psf: PSF file must have identical dimensions for all HDUs" );
  }
  
  hdu = fits::img_seek ( fp, toy_psf_hdu_amp );
  hdus_[ toy_psf_hdu_amp ] = hdu;
  fits::img_dims ( fp, rows, cols );
  if ( ( rows != nspec_ ) || ( cols != nbins_ ) ) {
    MOAT_THROW( "toy_psf: PSF file must have identical dimensions for all HDUs" );
  }
  
  hdu = fits::img_seek ( fp, toy_psf_hdu_maj );
  hdus_[ toy_psf_hdu_maj ] = hdu;
  fits::img_dims ( fp, rows, cols );
  if ( ( rows != nspec_ ) || ( cols != nbins_ ) ) {
    MOAT_THROW( "toy_psf: PSF file must have identical dimensions for all HDUs" );
  }
  
  hdu = fits::img_seek ( fp, toy_psf_hdu_min );
  hdus_[ toy_psf_hdu_min ] = hdu;
  fits::img_dims ( fp, rows, cols );
  if ( ( rows != nspec_ ) || ( cols != nbins_ ) ) {
    MOAT_THROW( "toy_psf: PSF file must have identical dimensions for all HDUs" );
  }
  
  hdu = fits::img_seek ( fp, toy_psf_hdu_ang );
  hdus_[ toy_psf_hdu_ang ] = hdu;
  fits::img_dims ( fp, rows, cols );
  if ( ( rows != nspec_ ) || ( cols != nbins_ ) ) {
    MOAT_THROW( "toy_psf: PSF file must have identical dimensions for all HDUs" );
  }
  
  fits::close ( fp );
  
  pixcorr_ = (size_t) ( 5.0 * 4000.0 / (double) nspec_ ); 
  
}


harp::psf_toy::~psf_toy ( ) {
  
  cleanup();
  
}


size_t nspec ( ) { return nspec_; }
size_t specsize ( size_t specnum );
void lambda ( size_t specnum, data_vec & data );
size_t extent ( size_t firstspec, size_t lastspec, size_t firstbin, size_t lastbin );
void projection ( size_t firstX, size_t firstY, size_t lastX, size_t lastY, sparse_mat_view & data );

