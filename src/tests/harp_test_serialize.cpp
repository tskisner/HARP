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


void harp::test_serialize ( string const & datadir ) {

  cerr << "Testing general serialization..." << endl;

  // first, test simple serialization of a class.

  data_one one;

  string path = datadir + "/serialize_empty.xml.out";

  {
    ofstream ofs ( path.c_str() );
    boost::archive::binary_oarchive oa ( ofs );
    oa << one;
  }
  {
    ifstream ifs ( path.c_str() );
    boost::archive::binary_iarchive ia ( ifs );
    ia >> one;
  }


  cerr << "  (PASSED)" << endl;

  return;
}



