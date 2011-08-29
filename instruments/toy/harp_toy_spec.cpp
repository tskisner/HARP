// @COPYRIGHT@

#include <harp_internal.hpp>


using namespace std;
using namespace harp;

static const char * format_toy = "toy";

static const char * toy_spec_key_path = "path";
static const char * toy_spec_key_hdu = "hdu";
static const char * toy_spec_key_nspec = "nspec";
static const char * toy_spec_key_specsize = "specsize";


harp::spec_toy::spec_toy ( std::map < std::string, std::string > const & params ) : spec ( format_toy, params ) {
  
  map < std::string, std::string > :: const_iterator val;
  
  val = params.find( toy_spec_key_hdu );
  if ( val == params.end() ) {
    hdu_ = 1;
  } else {
    hdu_ = atoi ( val->second.c_str() );
  }
    
  val = params.find( toy_spec_key_path );
  
  if ( val == params.end() ) {
    
    // no path specified- must specify number of spectra
    
    val = params.find( toy_spec_key_nspec );
    if ( val == params.end() ) {
      MOAT_THROW( "toy_spec: must specify number of spectra if no path is given" );
    }
    
    nspec_ = (size_t) atoll ( val->second.c_str() );
    
    // no path specified- must specify spectrum size
    
    val = params.find( toy_spec_key_specsize );
    if ( val == params.end() ) {
      MOAT_THROW( "toy_spec: must specify spectrum size if no path is given" );
    }
    
    specsize_ = (size_t) atoll ( val->second.c_str() );
    
  } else {
    
    path_ = val->second;
    
    // read size from the FITS header
    
    fitsfile *fp;

    fits::open_read ( fp, path_ );

    fits::img_seek ( fp, hdu_ );
    
    size_t rows, cols;
    
    fits::img_dims ( fp, rows, cols );
    
    nspec_ = rows;
    specsize_ = cols;
    
    fits::close ( fp );
    
  }
  
  size_ = nspec_ * specsize_;
  
}


harp::spec_toy::~spec_toy ( ) {
  
  cleanup();
  
}


void harp::spec_toy::read ( vec_dense & data ) {
  
  fitsfile *fp;

  fits::open_read ( fp, path_ );

  fits::img_seek ( fp, hdu_ );
  
  fits::img_read ( fp, data );
  
  fits::close ( fp );
  
  return;
}


void harp::spec_toy::write ( string const & path, vec_dense & data ) {
  
  fitsfile *fp;

  fits::open_readwrite ( fp, path );
  
  int nh = fits::nhdus ( fp );

  if ( nh < hdu_ ) {
    while ( nh < hdu_ ) {
      fits::img_append ( fp, nspec_, specsize_ );
      ++nh;
    }
  } else {
    fits::img_seek ( fp, hdu_ );
  }
  
  fits::img_write ( fp, data );
  
  fits::close ( fp );
  
  return;
}


void harp::spec_toy::read_spectrum ( size_t spectrum, vec_dense & data ) {
  
  fitsfile *fp;

  fits::open_read ( fp, path_ );

  fits::img_seek ( fp, hdu_ );
  
  fits::img_read_row ( fp, spectrum + 1, data );
  
  fits::close ( fp );
  
  return;
}


void harp::spec_toy::write_spectrum ( size_t spectrum, string const & path, vec_dense & data ) {
  
  fitsfile *fp;

  fits::open_readwrite ( fp, path );
  
  int nh = fits::nhdus ( fp );

  if ( nh < hdu_ ) {
    while ( nh < hdu_ ) {
      fits::img_append ( fp, nspec_, specsize_ );
      ++nh;
    }
  } else {
    fits::img_seek ( fp, hdu_ );
  }

  fits::img_write_row ( fp, spectrum + 1, data );
  
  fits::close ( fp );
  
  return;
}



