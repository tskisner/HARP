#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdio>

#include <harp_test.hpp>

#include <boost/random.hpp>

extern "C" {
  #include <unistd.h>
  #include <sys/stat.h>
}

using namespace std;
using namespace harp;



harp::data_base * harp::data_base::create ( boost::property_tree::ptree const & props ) {

  string format = props.get < string > ( "format" );

  if ( format == "one" ) {
    return static_cast < harp::data_base * > ( new harp::data_one ( props ) );
  }

  std::ostringstream o;
  o << "Cannot create data of unknown format (" << format << ")";
  HARP_THROW( o.str().c_str() );

  return NULL;

}


class fake {

  friend class boost::serialization::access;

  public :
    double d_blah;
    size_t s_blah;

  private :

    template < class Archive >
    void serialize ( Archive & ar, const unsigned int version ) {
      ar & BOOST_SERIALIZATION_NVP(d_blah);
      ar & BOOST_SERIALIZATION_NVP(s_blah);
      return;
    }
};



BOOST_CLASS_EXPORT(fake)
BOOST_CLASS_EXPORT(data_base)
BOOST_CLASS_EXPORT(data_one)


void harp::test_serialize ( string const & datadir ) {

  cerr << "Testing general serialization..." << endl;

  // first, test simple serialization of a class.

  fake one;
  one.d_blah = 1.0;
  one.s_blah = 1;

  string path = datadir + "/serialize_fake.xml.out";

  {
    ofstream ofs ( path.c_str() );    
    boost::archive::xml_oarchive oa ( ofs );
    oa << BOOST_SERIALIZATION_NVP(one);
  }

  {
    ifstream ifs ( path.c_str() );
    boost::archive::xml_iarchive ia ( ifs );
    ia >> BOOST_SERIALIZATION_NVP(one);
  }

  data_one d_one;

  path = datadir + "/serialize_data-one_empty.xml.out";

  {
    ofstream ofs ( path.c_str() );    
    boost::archive::xml_oarchive oa ( ofs );
    oa << BOOST_SERIALIZATION_NVP(d_one);
  }

  {
    ifstream ifs ( path.c_str() );
    boost::archive::xml_iarchive ia ( ifs );
    ia >> BOOST_SERIALIZATION_NVP(d_one);
  }

  boost::property_tree::ptree props;

  props.put ( "format", "one" );
  props.put ( "sizet", 2 );
  props.put ( "double", 2.0 );

  data_one d_two ( props );

  path = datadir + "/serialize_data-one_props.xml.out";

  {
    ofstream ofs ( path.c_str() );    
    boost::archive::xml_oarchive oa ( ofs );
    oa << BOOST_SERIALIZATION_NVP(d_two);
  }

  {
    ifstream ifs ( path.c_str() );
    boost::archive::xml_iarchive ia ( ifs );
    ia >> BOOST_SERIALIZATION_NVP(d_two);
  }

  data_base_p dp ( data_base::create ( props ) );

  path = datadir + "/serialize_data_shptr.xml.out";

  {
    ofstream ofs ( path.c_str() );    
    boost::archive::xml_oarchive oa ( ofs );
    oa << BOOST_SERIALIZATION_NVP(dp);
  }

  {
    ifstream ifs ( path.c_str() );
    boost::archive::xml_iarchive ia ( ifs );
    ia >> BOOST_SERIALIZATION_NVP(dp);
  }
  

  cerr << "  (PASSED)" << endl;

  return;
}



