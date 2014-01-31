// @COPYRIGHT@

#include <harp_data_internal.hpp>


using namespace std;
using namespace harp;


harp::spec::spec ( std::string const & name, boost::property_tree::ptree const & props ) {
  props_ = props;
  name_ = name;
}


string harp::spec::name ( ) const {
  return name_;
}

