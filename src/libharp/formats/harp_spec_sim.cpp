// @COPYRIGHT@

#include <harp_internal.hpp>

using namespace std;
using namespace harp;


static const char * spec_sim_key_nspec = "nspec";
static const char * spec_sim_key_nlambda = "nlambda";
static const char * spec_sim_key_firstlambda = "first_lambda";
static const char * spec_sim_key_lastlambda = "last_lambda";
static const char * spec_sim_key_back = "back";
static const char * spec_sim_key_atm = "atm";
static const char * spec_sim_key_obj = "obj";
static const char * spec_sim_key_atmspace = "atmspace";
static const char * spec_sim_key_skymod = "skymod";


harp::spec_sim::spec_sim ( boost::property_tree::ptree const & props ) : spec ( props ) {
  
  nspec_ = props.get < size_t > ( spec_sim_key_nspec );

  nlambda_ = props.get < size_t > ( spec_sim_key_nlambda );

  first_lambda_ = props.get < double > ( spec_sim_key_firstlambda );

  last_lambda_ = props.get < double > ( spec_sim_key_lastlambda );  

  background_ = props.get < double > ( spec_sim_key_back );

  atmpeak_ = props.get < double > ( spec_sim_key_atm );
  
  objpeak_ = props.get < double > ( spec_sim_key_obj );
  
  atmspace_ = props.get < size_t > ( spec_sim_key_atmspace );

  objspace_ = (size_t)( (double)atmspace_ * 1.625 );

  skymod_ = props.get < size_t > ( spec_sim_key_skymod );
  
  size_ = nspec_ * nlambda_;
  
}


harp::spec_sim::~spec_sim ( ) {
  
  cleanup();
  
}


boost::property_tree::ptree harp::spec_sim::serialize ( ) {
  boost::property_tree::ptree ret;

  ret.put ( "format", spec::format() );


  return ret;
}


void harp::spec_sim::read ( matrix_dist & data, std::vector < double > & lambda, std::vector < bool > & sky ) {

  double PI = std::atan2 ( 0.0, -1.0 );

  lambda.resize ( nlambda_ );
  sky.resize ( nspec_ );

  double incr = (last_lambda_ - first_lambda_) / (double)( nlambda_ - 1 );

  for ( size_t j = 0; j < nlambda_; ++j ) {
    lambda[j] = first_lambda_ + incr * (double)j;
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

      if ( ( j % atmspace_ ) == 0 ) {
        val += atmpeak_;
      }

      if ( ( ( j + objoff ) % objspace_ ) == 0 ) {
        val += objpeak_;
      }

      data.Set ( bin, 0, val );

      ++bin;
    }
  }

  
  return;
}


void harp::spec_sim::write ( std::string const & path, matrix_dist & data, std::vector < double > const & lambda, std::vector < bool > const & sky ) {

  HARP_THROW( "sim spec format does not support writing" );
  
  return;
}



