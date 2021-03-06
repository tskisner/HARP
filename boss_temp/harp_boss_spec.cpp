/*
  High Performance Astrophysical Reconstruction and Processing (HARP)
  (c) 2014-2015, The Regents of the University of California, 
  through Lawrence Berkeley National Laboratory.  See top
  level LICENSE file for details.
*/

#include <harp_internal.hpp>


using namespace std;
using namespace harp;

static const char * format_boss = "boss";

static const char * boss_spec_key_path = "path";
static const char * boss_spec_key_hdu = "hdu";
static const char * boss_spec_key_nspec = "nspec";
static const char * boss_spec_key_specsize = "specsize";


harp::spec_boss::spec_boss ( std::map < std::string, std::string > const & params ) : spec ( format_boss, params ) {
  
  map < std::string, std::string > :: const_iterator val;
  
  val = params.find( boss_spec_key_hdu );
  if ( val == params.end() ) {
    hdu_ = 1;
  } else {
    hdu_ = atoi ( val->second.c_str() );
  }
    
  val = params.find( boss_spec_key_path );
  
  if ( val == params.end() ) {
    
    // no path specified- must specify number of spectra
    
    val = params.find( boss_spec_key_nspec );
    if ( val == params.end() ) {
      MOAT_THROW( "boss_spec: must specify number of spectra if no path is given" );
    }
    
    nspec_ = (size_t) atoll ( val->second.c_str() );
    
    // no path specified- must specify spectrum size
    
    val = params.find( boss_spec_key_specsize );
    if ( val == params.end() ) {
      MOAT_THROW( "boss_spec: must specify spectrum size if no path is given" );
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


harp::spec_boss::~spec_boss ( ) {
  
  cleanup();
  
}


void harp::spec_boss::read ( harp::data_vec_view & data ) {
  
  fitsfile *fp;

  fits::open_read ( fp, path_ );

  fits::img_seek ( fp, hdu_ );
  
  fits::img_read ( fp, data );
  
  fits::close ( fp );
  
  return;
}


void harp::spec_boss::write ( string const & path, harp::data_vec_view & data ) {
  
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


void harp::spec_boss::read_spectrum ( size_t spectrum, harp::data_vec_view & data ) {
  
  fitsfile *fp;

  fits::open_read ( fp, path_ );

  fits::img_seek ( fp, hdu_ );
  
  fits::img_read_row ( fp, spectrum + 1, data );
  
  fits::close ( fp );
  
  return;
}


void harp::spec_boss::write_spectrum ( size_t spectrum, string const & path, harp::data_vec_view & data ) {
  
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



