/*
  High Performance Astrophysical Reconstruction and Processing (HARP)
  (c) 2014-2015, The Regents of the University of California, 
  through Lawrence Berkeley National Laboratory.  See top
  level LICENSE file for details.
*/


#include <harp/static_plugins.hpp>

using namespace std;
using namespace harp;


static const char * targets_sim_key_nobj = "nobj";
static const char * targets_sim_key_skymod = "skymod";


harp::targets_sim::targets_sim ( ) : targets () {
  nobjects_ = 0;
  objects_.clear();
  skymod_ = 0;
}


harp::targets_sim::targets_sim ( boost::property_tree::ptree const & props ) : targets ( "sim", props ) {

  nobjects_ = props.get < size_t > ( targets_sim_key_nobj );

  skymod_ = props.get < size_t > ( targets_sim_key_skymod );

  objects_.clear();

  for ( size_t i = 0; i < nobjects_; ++i ) {
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


size_t harp::targets_sim::n_objects ( ) const {
  return nobjects_;
}


std::vector < object_p > harp::targets_sim::objects ( ) const {
  return objects_;
}


BOOST_CLASS_EXPORT(harp::targets_sim)

targets * harp::targets_sim_create ( boost::property_tree::ptree const & props ) {
  return new targets_sim ( props );
}

