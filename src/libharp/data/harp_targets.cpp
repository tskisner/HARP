// @COPYRIGHT@

#include <harp_data_internal.hpp>


using namespace std;
using namespace harp;

object_type harp::object_str2type ( std::string const & in ) {
  object_type ret;
  if ( in.compare ( "CALIB" ) == 0 ) {
    ret = OBJECT_CALIB;
  } else if ( in.compare ( "SKY" ) == 0 ) {
    ret = OBJECT_SKY;
  } else if ( in.compare ( "LRG" ) == 0 ) {
    ret = OBJECT_LRG;
  } else if ( in.compare ( "ELG" ) == 0 ) {
    ret = OBJECT_ELG;
  } else if ( in.compare ( "QSO" ) == 0 ) {
    ret = OBJECT_QSO;
  } else if ( in.compare ( "STAR" ) == 0 ) {
    ret = OBJECT_STAR;
  } else {
  	ret = OBJECT_UNKNOWN;
  }
  return ret;
}

std::string harp::object_type2str ( object_type const & in ) {
  string ret;
  switch ( in ) {
    case OBJECT_CALIB :
      ret = "CALIB";
      break;
    case OBJECT_SKY :
      ret = "SKY";
      break;
    case OBJECT_LRG :
      ret = "LRG";
      break;
    case OBJECT_ELG :
      ret = "ELG";
      break;
    case OBJECT_QSO :
      ret = "QSO";
      break;
    case OBJECT_STAR :
      ret = "STAR";
      break;
    default :
      ret = "UNKNOWN";
  }
  return ret;
}


harp::targets::targets ( std::string const & type, boost::property_tree::ptree const & props ) {
  props_ = props;
  type_ = type;
}


string harp::targets::type ( ) const {
  return type_;
}


