/*
  High Performance Astrophysical Reconstruction and Processing (HARP)
  (c) 2014-2015, The Regents of the University of California, 
  through Lawrence Berkeley National Laboratory.  See top
  level LICENSE file for details.
*/

#include <harp_data_internal.hpp>


using namespace std;
using namespace harp;


harp::image::image ( std::string const & type, boost::property_tree::ptree const & props ) {
  props_ = props;
  type_ = type;
}


string harp::image::type ( ) const {
  return type_;
}


