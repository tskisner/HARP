// @COPYRIGHT@

// Test that we can instantiate our plugin.

#include <iostream>
#include <cstdio>
#include <string>

#include <harp.hpp>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree_serialization.hpp>


using namespace std;
using namespace harp;


int main ( int argc, char *argv[] ) {

  boost::property_tree::ptree params;

  cout << endl;
  cout << "Testing plugin registration..." << endl;

  bool reg_debug = true;
  plugin_registry & reg = plugin_registry::get( reg_debug );

  cout << "Testing instantiation of example spec type" << endl;
  spec_p testspec ( reg.create_spec ( "example", params ) );

  cout << "Testing instantiation of example psf type" << endl;
  psf_p testpsf ( reg.create_psf ( "example", params ) );

  cout << "Testing instantiation of example image type" << endl;
  image_p testimage ( reg.create_image ( "example", params ) );

  cout << endl;
  cout << "PASSED" << endl << endl;

  return 0;
}



