// @COPYRIGHT@

#include <harp_internal.hpp>


using namespace std;
using namespace harp;


void harp::specter_read_sky ( fitsfile * fp, std::vector < bool > & sky ) {

  size_t nspec;

  vector < string > colnames;
  fits::bin_info ( fp, nspec, colnames );

  if ( colnames[0] != "OBJTYPE" ) {
    HARP_THROW( "specter format object table, expected first column to be \"OBJTYPE\"" );
  }
  if ( colnames[1] != "Z" ) {
    HARP_THROW( "specter format object table, expected second column to be \"Z\"" );
  }
  if ( colnames[2] != "O2FLUX" ) {
    HARP_THROW( "specter format object table, expected third column to be \"O2FLUX\"" );
  }

  vector < string > skycolnames(1);
  skycolnames[0] = "OBJTYPE";

  vector < int > skycols;
  skycols = fits::bin_columns ( fp, skycolnames );

  vector < string > objnames;
  fits::bin_read_column_strings ( fp, 0, nspec - 1, skycols[0], objnames );

  fits::close ( fp );

  sky.resize ( nspec );

  for ( size_t i = 0; i < nspec; ++i ) {
    if ( objnames[i] == "SKY" ) {
      sky[i] = true;
    } else {
      sky[i] = false;
    }
  }

  return;
}


void harp::specter_write_sky ( fitsfile * fp, std::vector < bool > const & sky ) {

  vector < string > colnames ( sky.size() );
  vector < string > coltypes ( sky.size() );
  vector < string > colunits ( sky.size() );

  colnames[0] = "OBJTYPE";
  coltypes[0] = "6A";
  colunits[0] = "None";

  colnames[0] = "Z";
  coltypes[0] = "1E";
  colunits[0] = "None";

  colnames[0] = "O2FLUX";
  coltypes[0] = "1E";
  colunits[0] = "None";

  fits::bin_create ( fp, string("TARGETINFO"), sky.size(), colnames, coltypes, colunits );

  vector < string > objnames ( sky.size() );

  for ( size_t i = 0; i < sky.size(); ++i ) {
    if ( sky[i] ) {
      objnames[i] == "SKY";
    } else {
      objnames[i] = "Unknown";
    }
  }

  fits::bin_write_column_strings ( fp, 0, sky.size() - 1, 1, objnames );

  return;
}



static const char * spec_specter_key_path = "path";
static const char * spec_specter_key_objonly = "objonly";
static const char * spec_specter_key_nspec = "nspec";
static const char * spec_specter_key_nlambda = "nlambda";


harp::spec_specter::spec_specter ( boost::property_tree::ptree const & props ) : spec ( props ) {

  path_ = props.get < string > ( spec_specter_key_path, "" );

  boost::optional < string > objval = props.get_optional < string > ( spec_specter_key_objonly );
  string obj = boost::get_optional_value_or ( objval, "FALSE" );
  objonly_ = ( obj == "TRUE" );

  if ( path_ == "" ) {

    nspec_ = props.get < size_t > ( spec_specter_key_nspec );

    nlambda_ = props.get < size_t > ( spec_specter_key_nlambda );

    spechdu_ = 1;

    lambdahdu_ = 2;

    targethdu_ = 3;

  } else {

    fitsfile * fp;

    fits::open_read ( fp, path_ );

    if ( objonly_ ) {
      spechdu_ = fits::img_seek ( fp, "EXTNAME", "OBJPHOT" );
    } else {
      spechdu_ = fits::img_seek ( fp, "EXTNAME", "FLUX" );
    }

    fits::img_dims ( fp, nspec_, nlambda_ );

    lambdahdu_ = fits::img_seek ( fp, "EXTNAME", "WAVELENGTH" );
    targethdu_ = fits::bin_seek ( fp, "EXTNAME", "TARGETINFO" );

    fits::close ( fp );

  }

  nglobal_ = nspec_ * nlambda_;
  
}


harp::spec_specter::~spec_specter ( ) {
  
}


void harp::spec_specter::read ( vector_double & data, vector_double & lambda, std::vector < bool > & sky ) {

  data.resize ( nglobal_ );
  data.clear();

  lambda.resize ( nlambda_ );
  sky.resize ( nspec_ );

  fitsfile * fp;

  fits::open_read ( fp, path_ );

  // read the spectral data

  fits::img_seek ( fp, spechdu_ );   
  fits::img_read ( fp, data );

  // read the wavelength vector

  fits::img_seek ( fp, lambdahdu_ );
  fits::img_read ( fp, lambda );

  // read the sky flag

  fits::bin_seek ( fp, targethdu_ );
  specter_read_sky ( fp, sky );

  fits::close ( fp );

  return;
}


void harp::spec_specter::write ( std::string const & path, vector_double & data, vector_double const & lambda, std::vector < bool > const & sky ) {

  fitsfile * fp;
    
  fits::create ( fp, path );

  fits::img_append < double > ( fp, nspec_, nlambda_ );
  fits::write_key ( fp, "EXTNAME", "FLUX", "" );
  fits::img_write ( fp, data );

  fits::img_append < double > ( fp, 1, nlambda_ );
  fits::write_key ( fp, "EXTNAME", "WAVELENGTH", "" );
  fits::img_write ( fp, lambda );

  specter_write_sky ( fp, sky );

  fits::close ( fp );

  return;
}



