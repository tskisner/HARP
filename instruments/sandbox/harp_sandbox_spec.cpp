// @COPYRIGHT@

#include <harp_internal.hpp>


using namespace std;
using namespace harp;

static const char * format_sandbox = "sandbox";

static const char * sandbox_spec_key_path = "path";
static const char * sandbox_spec_key_hdu = "hdu";
static const char * sandbox_spec_key_nspec = "nspec";
static const char * sandbox_spec_key_specsize = "specsize";


harp::spec_sandbox::spec_sandbox ( boost::ptree const & params ) : spec ( format_sandbox, params ) {
  
  hdu_ = props.get ( sandbox_spec_key_hdu, 1 );

  path_ = props.get ( sandbox_spec_key_path, "" );

  if ( path_ == "" ) {

    // no path specified- must specify spectra size

    nspec_ = props.get < size_t > ( sandbox_spec_key_nspec );

    specsize_ = props.get < size_t > ( sandbox_spec_key_specsize );

  } else {
    
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


harp::spec_sandbox::~spec_sandbox ( ) {
  
  cleanup();
  
}


boost::ptree harp::spec_sandbox::serialize ( ) {
  boost::ptree ret;

  ret.put ( "format", spec::format() );

  if ( hdu_ != 1 ) {
    ret.put ( sandbox_spec_key_hdu, hdu_ );
  }

  if ( path_ == "" ) {
    ret.put ( sandbox_spec_key_nspec, nspec_ );
    ret.put ( sandbox_spec_key_specsize, specsize_ );
  } else {
    ret.put ( sandbox_spec_key_path, path_ );
  }

  return ret;
}


void harp::spec_sandbox::read ( vec_dense & data ) {
  
  fitsfile *fp;

  fits::open_read ( fp, path_ );

  fits::img_seek ( fp, hdu_ );
  
  fits::img_read ( fp, data );
  
  fits::close ( fp );
  
  return;
}


void harp::spec_sandbox::write ( string const & path, vec_dense & data ) {
  
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


void harp::spec_sandbox::read_spectrum ( size_t spectrum, vec_dense & data ) {
  
  fitsfile *fp;

  fits::open_read ( fp, path_ );

  fits::img_seek ( fp, hdu_ );
  
  fits::img_read_row ( fp, spectrum + 1, data );
  
  fits::close ( fp );
  
  return;
}


void harp::spec_sandbox::write_spectrum ( size_t spectrum, string const & path, vec_dense & data ) {
  
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



