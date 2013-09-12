// @COPYRIGHT@

#include <harp_internal.hpp>


using namespace std;
using namespace harp;

static const char * image_fits_key_path = "path";
static const char * image_fits_key_signal = "signal";
static const char * image_fits_key_sky = "sky";
static const char * image_fits_key_noise = "noise";
static const char * image_fits_key_rows = "rows";
static const char * image_fits_key_cols = "cols";


harp::image_fits::image_fits ( boost::property_tree::ptree const & props ) : image ( props ) {

  sighdu_ = props.get ( image_fits_key_signal, 1 );

  nsehdu_ = props.get ( image_fits_key_noise, 2 );

  skyhdu_ = props.get ( image_fits_key_sky, 3 );

  path_ = props.get ( image_fits_key_path, "" );

  if ( path_ == "" ) {

    // no path specified- must specify rows / cols

    rows_ = props.get < size_t > ( image_fits_key_rows );

    cols_ = props.get < size_t > ( image_fits_key_cols );

  } else {
    
    // read rows / cols from the FITS header

    fitsfile *fp;
    
    fits::open_read ( fp, path_ );

    fits::img_seek ( fp, sighdu_ );
    
    fits::img_dims ( fp, rows_, cols_ );

    size_t rowcheck;
    size_t colcheck;

    fits::img_seek ( fp, nsehdu_ );
    
    fits::img_dims ( fp, rowcheck, colcheck );

    if ( ( rowcheck != rows_ ) || ( colcheck != cols_ ) ) {
      HARP_THROW( "noise variance dimensions do not match data dimensions" );
    }

    // read the sky information

    size_t nspec;

    vector < string > skycolnames(1);
    skycolnames[0] = "OBJTYPE";

    fits::bin_seek ( fp, skyhdu_ );

    fits::bin_info ( fp, nspec, skycolnames );

    vector < int > skycols;
    skycols = fits::bin_columns ( fp, skycolnames );

    vector < string > objnames;
    fits::bin_read_column_strings ( fp, 0, nspec - 1, cols[0], objnames );

    fits::close ( fp );

    sky_.resize ( nspec );

    for ( size_t i = 0; i < nspec; ++i ) {
      if ( objnames[i] == "SKY" ) {
        sky_[i] = true;
      } else {
        sky_[i] = false;
      }
    }
    
  }  
  
}


harp::image_fits::~image_fits ( ) {

}


void read ( vector_double & data, vector_double & invvar, std::vector < bool > & sky ) {

  fitsfile *fp;

  fits::open_read ( fp, path_ );

  fits::img_seek ( fp, sighdu_ );
    
  fits::img_read ( fp, data );

  fits::img_seek ( fp, nsehdu_ );
    
  fits::img_read ( fp, invvar );

  sky = sky_;
    
  fits::close ( fp );

  return;
}


void write ( std::string const & path, vector_double & data, vector_double & invvar, std::vector < bool > & sky ) {

  fitsfile *fp;
    
  fits::open_readwrite ( fp, path );
  
  int nh = fits::nhdus ( fp );

  if ( nh < sighdu_ ) {
    while ( nh < sighdu_ ) {
      fits::img_append ( fp, rows_, cols_ );
      ++nh;
    }
  } else {
    fits::img_seek ( fp, sighdu_ );
  }
  
  fits::img_write ( fp, data );

  if ( nh < nsehdu_ ) {
    while ( nh < nsehdu_ ) {
      fits::img_append ( fp, rows_, cols_ );
      ++nh;
    }
  } else {
    fits::img_seek ( fp, nsehdu_ );
  }
  
  fits::img_write ( fp, invvar );

  vector < string > colnames ( sky_.size() );
  vector < string > coltypes ( sky_.size() );
  vector < string > colunits ( sky_.size() );

  colnames[0] = "OBJTYPE";
  coltypes[0] = "6A";
  colunits[0] = "None";

  colnames[0] = "Z";
  coltypes[0] = "1E";
  colunits[0] = "None";

  colnames[0] = "O2FLUX";
  coltypes[0] = "1E";
  colunits[0] = "None";

  if ( nh < skyhdu_ ) {
    while ( nh < skyhdu_ ) {
      bin_create ( fp, string("TARGETINFO"), sky_.size(), colnames, coltypes, colunits );
      ++nh;
    }
  } else {
    fits::bin_seek ( fp, skyhdu_ );
  }

  vector < string > objnames ( sky_.size() );

  for ( size_t i = 0; i < sky_.size(); ++i ) {
    if ( sky_[i] ) {
      objnames[i] == "SKY";
    } else {
      objnames[i] = "Unknown";
    }
  }

  fits::bin_write_column_strings ( fp, 0, sky.size() - 1, 1, objnames );
  
  fits::close ( fp );

  return;
}

