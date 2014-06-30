// @COPYRIGHT@


#include <harp/static_plugins.hpp>

using namespace std;
using namespace harp;


static const char * targets_sim_key_nobj = "nobj";
static const char * targets_sim_key_skymod = "skymod";

harp::targets_sim::targets_sim ( boost::property_tree::ptree const & props ) : targets ( "sim", props ) {

  nobjects_ = props.get < size_t > ( targets_sim_key_nobj );

  skymod_ = props.get < size_t > ( targets_sim_key_skymod );

  objects_.clear();

  for ( size_t i = 0; i < nrows; ++i ) {
    object_type type;
    if ( i % skymod_ == 0 ) {
      type = OBJECT_SKY;
    } else {
      type = OBJECT_UNKNOWN;
    }
    objects_.push_back ( object_p ( new object ( type, "" ) ) );
  }
  
}


harp::targets_sim::~targets_sim ( ) {
  
}


BOOST_CLASS_EXPORT(harp::targets_sim)

targets * harp::targets_sim_create ( boost::property_tree::ptree const & props ) {
  return new targets_sim ( props );
}

