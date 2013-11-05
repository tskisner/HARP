// @COPYRIGHT@

#ifndef HARP_HPP
#define HARP_HPP

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree_serialization.hpp>

#include <boost/smart_ptr/shared_ptr.hpp>
#include <boost/smart_ptr/weak_ptr.hpp>
#include <boost/enable_shared_from_this.hpp>

#include <boost/serialization/string.hpp>
#include <boost/serialization/version.hpp>
#include <boost/serialization/split_member.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/nvp.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/export.hpp>
#include <boost/serialization/assume_abstract.hpp>

#include <boost/serialization/list.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/map.hpp>

#include <boost/foreach.hpp>

#include <harp/math.hpp>


namespace harp {

  // utilities

  std::string ptree_quote ( const std::string & s );
  void ptree_print ( const boost::property_tree::ptree & pt, int level );
  void ptree_print ( const boost::property_tree::ptree & pt );
  
}

#include <harp/fits.hpp>

#include <harp/image.hpp>
#include <harp/spec.hpp>
#include <harp/psf.hpp>

#include <harp/formats_image.hpp>
#include <harp/formats_spec.hpp>
#include <harp/formats_psf.hpp>

#include <harp/extraction.hpp>

#include <harp/sky.hpp>

#endif



