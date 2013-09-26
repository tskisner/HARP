// @COPYRIGHT@

#include <harp_internal.hpp>

using namespace std;
using namespace harp;


static const char * spec_sim_key_nspec = "nspec";
static const char * spec_sim_key_nlambda = "lambda_n";
static const char * spec_sim_key_firstlambda = "lambda_start";
static const char * spec_sim_key_lastlambda = "lambda_stop";
static const char * spec_sim_key_back = "back";
static const char * spec_sim_key_atm = "atm";
static const char * spec_sim_key_obj = "obj";
static const char * spec_sim_key_atmspace = "atmspace";
static const char * spec_sim_key_skymod = "skymod";


harp::spec_sim::spec_sim ( boost::property_tree::ptree const & props ) : spec ( props ) {
  
  nspec_ = props.get < size_t > ( spec_sim_key_nspec );

  nlambda_ = props.get < size_t > ( spec_sim_key_nlambda, 4697 );

  first_lambda_ = props.get < double > ( spec_sim_key_firstlambda, 7460.0 );

  last_lambda_ = props.get < double > ( spec_sim_key_lastlambda, 9808.0 );  

  background_ = props.get < double > ( spec_sim_key_back );

  atmpeak_ = props.get < double > ( spec_sim_key_atm );
  
  objpeak_ = props.get < double > ( spec_sim_key_obj );
  
  atmspace_ = props.get < size_t > ( spec_sim_key_atmspace );

  objspace_ = (size_t)( (double)atmspace_ * 1.625 );

  skymod_ = props.get < size_t > ( spec_sim_key_skymod );
  
  size_ = nspec_ * nlambda_;
  
}


harp::spec_sim::~spec_sim ( ) {
  
}


void harp::spec_sim::read ( vector_double & data, vector_double & lambda, std::vector < bool > & sky ) {

  double PI = std::atan2 ( 0.0, -1.0 );

  data.resize ( size_ );
  lambda.resize ( nlambda_ );
  sky.resize ( nspec_ );

  double incr = (last_lambda_ - first_lambda_) / (double)( nlambda_ - 1 );

  for ( size_t j = 0; j < nlambda_; ++j ) {
    lambda[j] = first_lambda_ + incr * (double)j;
  }

  size_t bin = 0;

  size_t halfspace = (size_t)( atmspace_ / 2 );

  for ( size_t i = 0; i < nspec_; ++i ) {

    if ( i % skymod_ == 0 ) {
      sky[ i ] = true;
    } else {
      sky[ i ] = false;
    }

    size_t objoff = 2 * i;

    for ( size_t j = 0; j < nlambda_; ++j ) {

      double val = 0.0;

      val += background_ * sin ( 3.0 * (double)j / ( (double)atmspace_ * 2.0 * PI ) ) * sin ( 7.0 * (double)j / ( (double)atmspace_ * 2.0 * PI ) ) * sin ( 11.0 * (double)j / ( (double)atmspace_ * 2.0 * PI ) );

      val += 2.0 * background_;

      if ( ( ( j + halfspace ) % atmspace_ ) == 0 ) {
        val += atmpeak_;
      }

      if ( ! sky[i] ) {
        if ( ( ( j + objoff ) % objspace_ ) == 0 ) {
          val += objpeak_;
        }
      }

      data[ bin ] = val;

      ++bin;
    }
  }

  
  return;
}


void harp::spec_sim::write ( std::string const & path, vector_double & data, vector_double const & lambda, std::vector < bool > const & sky ) {

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


void harp::spec_sim::sky_truth ( vector_double & data ) {

  double PI = std::atan2 ( 0.0, -1.0 );

  size_t halfspace = (size_t)( atmspace_ / 2 );

  size_t nreduced = 0;

  for ( size_t i = 0; i < nspec_; ++i ) {
    if ( i % skymod_ != 0 ) {
      ++nreduced;
    }
  }

  ++nreduced;

  size_t nbins = nreduced * nlambda_;

  data.resize ( nbins );

  size_t bin = 0;

  for ( size_t i = 0; i < nspec_; ++i ) {

    if ( i % skymod_ != 0 ) {

      size_t objoff = 2 * i;

      for ( size_t j = 0; j < nlambda_; ++j ) {

        double val = 0.0;

        if ( ( ( j + objoff ) % objspace_ ) == 0 ) {
          val += objpeak_;
        }

        data[ bin ] = val;

        ++bin;
      }

    }

  }

  for ( size_t j = 0; j < nlambda_; ++j ) {

    double val = 0.0;

    val += background_ * sin ( 3.0 * (double)j / ( (double)atmspace_ * 2.0 * PI ) ) * sin ( 7.0 * (double)j / ( (double)atmspace_ * 2.0 * PI ) ) * sin ( 11.0 * (double)j / ( (double)atmspace_ * 2.0 * PI ) );

    val += 2.0 * background_;

    if ( ( ( j + halfspace ) % atmspace_ ) == 0 ) {
      val += atmpeak_;
    }

    data[ bin ] = val;

    ++bin;
  }

  return;
}



