/*
  High Performance Astrophysical Reconstruction and Processing (HARP)
  (c) 2014-2015, The Regents of the University of California, 
  through Lawrence Berkeley National Laboratory.  See top
  level LICENSE file for details.
*/

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


harp::object::object ( ) {
  type_ = OBJECT_UNKNOWN;
  name_ = "";
}


harp::object::object ( object_type type, std::string name ) {
  type_ = type;
  name_ = name;
}

harp::object::~object ( ) {
}


object_type harp::object::type ( ) const {
  return type_;
}


std::string harp::object::name ( ) const {
  return name_;
}


void harp::object::set_type ( object_type t ) {
  type_ = t;
  return;
}


void harp::object::set_name ( std::string const & n ) {
  name_ = n;
  return;
}


bool harp::object::is_sky ( ) const {
  return ( type_ == OBJECT_SKY );
}


harp::targets::targets ( ) {
  type_ = "";
}


harp::targets::targets ( std::string const & type, boost::property_tree::ptree const & props ) {
  props_ = props;
  type_ = type;
}


harp::targets::~targets ( ) {
}


size_t harp::targets::n_objects ( ) const {
  HARP_THROW( "fell through to virtual method" );
  return 0;
}


std::vector < object_p > harp::targets::objects ( ) const {
  HARP_THROW( "fell through to virtual method" );
  return std::vector < object_p > ();
}


string harp::targets::type ( ) const {
  return type_;
}


