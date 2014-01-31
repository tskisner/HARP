// @COPYRIGHT@

#include <harp_internal.hpp>

#include <harp/plugin.hpp>
#ifdef HAVE_BOOST_MPI_HPP
#include <harp/plugin_mpi.hpp>
#endif

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


boost::property_tree::ptree harp::spec_sim::metadata ( ) const {

  return boost::property_tree::ptree();
}


void harp::spec_sim::values ( vector_double & data ) const {

  double PI = std::atan2 ( 0.0, -1.0 );

  data.resize ( size_ );

  size_t bin = 0;

  size_t halfspace = (size_t)( atmspace_ / 2 );

  for ( size_t i = 0; i < nspec_; ++i ) {

    bool dosky = ( i % skymod_ == 0 ) ? true : false;

    size_t objoff = 2 * i;

    for ( size_t j = 0; j < nlambda_; ++j ) {

      double val = 0.0;

      val += background_ * sin ( 3.0 * (double)j / ( (double)atmspace_ * 2.0 * PI ) ) * sin ( 7.0 * (double)j / ( (double)atmspace_ * 2.0 * PI ) ) * sin ( 11.0 * (double)j / ( (double)atmspace_ * 2.0 * PI ) );

      val += 2.0 * background_;

      if ( ( ( j + halfspace ) % atmspace_ ) == 0 ) {
        val += atmpeak_;
      }

      if ( ! dosky ) {
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


void harp::spec_sim::lambda ( vector_double & lambda_vals ) const {

  lambda_vals.resize ( nlambda_ );

  double incr = (last_lambda_ - first_lambda_) / (double)( nlambda_ - 1 );

  for ( size_t j = 0; j < nlambda_; ++j ) {
    lambda_vals[j] = first_lambda_ + incr * (double)j;
  }
  
  return;
}


void harp::spec_sim::targets ( std::vector < target > & target_list ) const {

  target_list.clear();

  for ( size_t i = 0; i < nspec_; ++i ) {

    if ( i % skymod_ == 0 ) {
      target_list.push_back ( target ( TARGET_SKY, "sky" ) );
    } else {
      target_list.push_back ( target ( TARGET_UNKNOWN, "object" ) );
    }

  }
  
  return;
}


void harp::spec_sim::sky_truth ( vector_double & data ) const {

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


BOOST_CLASS_EXPORT(harp::spec_sim)

