
#include <iostream>

#include <harp_test.hpp>

using namespace std;
using namespace harp;


int main ( int argc, char *argv[] ) {

  string datadir = "testdata";

  // run built-in tests
  
  cerr << endl;
  

  // run format tests
  
# include "harp_testcommands.cpp"

  cerr << endl;
  
  return 0;
}
