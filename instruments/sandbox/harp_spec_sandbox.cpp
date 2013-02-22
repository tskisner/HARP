// @COPYRIGHT@

#include <harp_internal.hpp>


using namespace std;
using namespace harp;

static const char * format_sandbox = "sandbox";

static const char * sandbox_spec_key_nspec = "nspec";
static const char * sandbox_spec_key_nlambda = "nlambda";
static const char * sandbox_spec_key_back = "back";
static const char * sandbox_spec_key_atm = "atm";
static const char * sandbox_spec_key_obj = "obj";
static const char * sandbox_spec_key_atmspace = "atmspace";
static const char * sandbox_spec_key_skymod = "skymod";


harp::spec_sandbox::spec_sandbox ( boost::property_tree::ptree const & props ) : spec ( format_sandbox, props ) {
  
  nspec_ = props.get < size_t > ( sandbox_spec_key_nspec );

  nlambda_ = props.get < size_t > ( sandbox_spec_key_nlambda );

  background_ = props.get < double > ( sandbox_spec_key_back );

  atmpeak_ = props.get < double > ( sandbox_spec_key_atm );
  
  objpeak_ = props.get < double > ( sandbox_spec_key_obj );
  
  atmspace_ = props.get < size_t > ( sandbox_spec_key_atmspace );

  objspace_ = (size_t)( (double)atmspace_ * 1.625 );

  skymod_ = props.get < size_t > ( sandbox_spec_key_skymod );
  
  size_ = nspec_ * nlambda_;
  
}


harp::spec_sandbox::~spec_sandbox ( ) {
  
  cleanup();
  
}


boost::property_tree::ptree harp::spec_sandbox::serialize ( ) {
  boost::property_tree::ptree ret;

  ret.put ( "format", spec::format() );


  return ret;
}


void harp::spec_sandbox::read ( matrix_dist & data, std::vector < double > & lambda, std::vector < bool > & sky ) {

  double PI = std::atan2 ( 0.0, -1.0 );

  lambda.resize ( nlambda_ );
  sky.resize ( nspec_ );

  for ( size_t j = 0; j < nlambda_; ++j ) {
    lambda[j] = 0.0;
  }

  size_t bin = 0;

  for ( size_t i = 0; i < nspec_; ++i ) {

    if ( i % skymod_ == 0 ) {
      sky[ i ] = true;
    } else {
      sky[ i ] = false;
    }

    size_t objoff = 2 * i;

    for ( size_t j = 0; j < nlambda_; ++j ) {

      double val = background_ * sin ( 3.0 * (double)j / ( (double)atmspace_ * 2.0 * PI ) ) * sin ( 7.0 * (double)j / ( (double)atmspace_ * 2.0 * PI ) ) * sin ( 11.0 * (double)j / ( (double)atmspace_ * 2.0 * PI ) );

      val += 2.0 * background_;

      if ( ( ( j + objoff ) % objspace_ ) == 0 ) {
        val += objpeak_;
      }

      data.Set ( bin, 0, val );

      ++bin;
    }
  }

  
  return;
}


void harp::spec_sandbox::write ( std::string const & path, matrix_dist & data, std::vector < double > const & lambda, std::vector < bool > const & sky ) {

  HARP_THROW( "sandbox spec format does not support writing" );
  
  return;
}



