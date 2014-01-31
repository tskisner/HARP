// @COPYRIGHT@

#include <harp_data_internal.hpp>


using namespace std;
using namespace harp;


harp::image::image ( std::string const & name, boost::property_tree::ptree const & props ) {
  props_ = props;
  name_ = name;
}


string harp::image::name ( ) const {
  return name_;
}


