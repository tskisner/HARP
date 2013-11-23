// @COPYRIGHT@

#include <harp_internal.hpp>


using namespace std;
using namespace harp;


void harp::specter_read_targets ( fitsfile * fp, std::vector < target > & target_list ) {

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

  target_list.clear();

  // FIXME: modify this to read all the object types supported by specter and
  // map them to HARP types.

  for ( size_t i = 0; i < nspec; ++i ) {
    if ( objnames[i] == "SKY" ) {
      target_list.push_back ( target ( TARGET_SKY, "sky" ) );
    } else {
      target_list.push_back ( target ( TARGET_UNKNOWN, objnames[i] ) );
    }
  }

  return;
}


void harp::specter_write_targets ( fitsfile * fp, std::vector < target > const & target_list ) {

  vector < string > colnames ( 3 );
  vector < string > coltypes ( 3 );
  vector < string > colunits ( 3 );

  colnames[0] = "OBJTYPE";
  coltypes[0] = "6A";
  colunits[0] = "None";

  colnames[1] = "Z";
  coltypes[1] = "1E";
  colunits[1] = "None";

  colnames[2] = "O2FLUX";
  coltypes[2] = "1E";
  colunits[2] = "None";

  fits::bin_create ( fp, string("TARGETINFO"), target_list.size(), colnames, coltypes, colunits );

  vector < string > objnames ( target_list.size() );

  for ( size_t i = 0; i < target_list.size(); ++i ) {
    if ( target_list[i].type() == TARGET_SKY ) {
      objnames[i] = "SKY";
    } else {
      objnames[i] = "NA";
    }
  }

  fits::bin_write_column_strings ( fp, 0, target_list.size() - 1, 0, objnames );

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


void harp::spec_specter::values ( vector_double & data ) const {

  data.resize ( nglobal_ );
  data.clear();

  fitsfile * fp;

  fits::open_read ( fp, path_ );

  // read the spectral data

  fits::img_seek ( fp, spechdu_ );   
  fits::img_read ( fp, data );

  fits::close ( fp );

  return;
}


void harp::spec_specter::lambda ( vector_double & lambda ) const {

  lambda.resize ( nlambda_ );

  fitsfile * fp;

  fits::open_read ( fp, path_ );

  // read the wavelength vector

  fits::img_seek ( fp, lambdahdu_ );
  fits::img_read ( fp, lambda );

  fits::close ( fp );

  return;
}


void harp::spec_specter::targets ( std::vector < target > & target_list ) const {

  fitsfile * fp;

  fits::open_read ( fp, path_ );

  // read the target list

  fits::bin_seek ( fp, targethdu_ );
  specter_read_targets ( fp, target_list );

  fits::close ( fp );

  return;
}


void harp::spec_specter::write ( std::string const & path, vector_double & data, vector_double const & lambda, std::vector < target > const & target_list ) {

  fitsfile * fp;
    
  fits::create ( fp, path );

  fits::img_append < double > ( fp, nspec_, nlambda_ );
  fits::write_key ( fp, "EXTNAME", "FLUX", "" );
  fits::img_write ( fp, data );

  fits::img_append < double > ( fp, 1, nlambda_ );
  fits::write_key ( fp, "EXTNAME", "WAVELENGTH", "" );
  fits::img_write ( fp, lambda );

  specter_write_targets ( fp, target_list );

  fits::close ( fp );

  return;
}


void harp::spec_specter::write ( std::string const & path, matrix_double & data, vector_double & lambda, std::vector < target > const & target_list ) {

  size_t nelem = nspec_ * nlambda_;

  if ( ( nlambda_ != data.size2() ) || ( nspec_ != data.size1() ) ) {
    HARP_THROW( "data size does not match spec dimensions" );
  }

  vector_double tempdata ( nelem );

  for ( size_t i = 0; i < nspec_; ++i ) {
    for ( size_t j = 0; j < nlambda_; ++j ) {
      tempdata[ i * nlambda_ + j ] = data( i, j );
    }
  }

  write ( path, tempdata, lambda, target_list );

  return;
}



