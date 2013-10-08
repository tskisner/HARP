#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include <harp_test.hpp>


extern "C" {
#include <unistd.h>
}

using namespace std;
using namespace harp;


void harp::test_spec_sim ( string const & datadir ) {

  cerr << "Testing simulated spec generation..." << endl;

  size_t nlambda = 50;
  size_t nspec = 100;
  size_t first_lambda = 8000.0;
  size_t last_lambda = 8001.1;

  size_t nbins = nspec * nlambda;

  boost::property_tree::ptree spec_props;
  spec_props.clear();
  spec_props.put ( "format", "sim" );
  spec_props.put ( "nspec", nspec );
  spec_props.put ( "lambda_n", nlambda );
  spec_props.put ( "lambda_start", first_lambda );
  spec_props.put ( "lambda_stop", last_lambda );
  spec_props.put ( "back", 10.0 );
  spec_props.put ( "atm", 500.0 );
  spec_props.put ( "obj", 80.0 );
  spec_props.put ( "atmspace", 12 );
  spec_props.put ( "skymod", 25 );
  spec_p testspec ( spec::create ( spec_props ) );

  // immediately serialize and restore, so that any issues with that process will impact the code that follows

  string serialpath = datadir + "/test_specsim_serialize.xml.out";
  {
    ofstream ofs ( serialpath.c_str() );
    boost::archive::xml_oarchive oa ( ofs );
    oa << BOOST_SERIALIZATION_NVP(testspec);
  }
  {
    ifstream ifs ( serialpath.c_str() );
    boost::archive::xml_iarchive ia ( ifs );
    ia >> BOOST_SERIALIZATION_NVP(testspec);
  }

  if ( nspec != testspec->n_spec() ) {
    cerr << "FAIL:  simulated spec nspec (" << testspec->n_spec() << ") does not match input (" << nspec << ")" << endl;
    exit(1);
  }

  if ( nlambda != testspec->n_lambda() ) {
    cerr << "FAIL:  simulated spec nlambda (" << testspec->n_lambda() << ") does not match input (" << nlambda << ")" << endl;
    exit(1);
  }

  vector_double data;
  vector_double lambda;
  vector < bool > sky;

  testspec->read( data, lambda, sky );

  string outfile = datadir + "/spec_sim.fits.out";

  testspec->write ( outfile, data, lambda, sky );

  spec_props.clear();
  spec_props.put ( "format", "specter" );
  spec_props.put ( "path", outfile );

  spec_p checkspec ( spec::create ( spec_props ) );

  if ( nspec != checkspec->n_spec() ) {
    cerr << "FAIL:  simulated spec written nspec (" << checkspec->n_spec() << ") does not match input (" << nspec << ")" << endl;
    exit(1);
  }

  if ( nlambda != checkspec->n_lambda() ) {
    cerr << "FAIL:  simulated spec written nlambda (" << checkspec->n_lambda() << ") does not match input (" << nlambda << ")" << endl;
    exit(1);
  }

  vector_double check_data;
  vector_double check_lambda;
  vector < bool > check_sky;

  checkspec->read( check_data, check_lambda, check_sky );

  for ( size_t i = 0; i < sky.size(); ++i ) {
    if ( check_sky[i] != sky[i] ) {
      cerr << "FAIL:  written sky element " << i << " does not match original value" << endl;
      exit(1);
    }
  }

  for ( size_t i = 0; i < lambda.size(); ++i ) {
    if ( fabs ( check_lambda[i] - lambda[i] ) / lambda[i] > 1.0e-5 ) {
      cerr << "FAIL:  written lambda element " << i << " does not match original value" << endl;
      exit(1);
    }
  }

  for ( size_t i = 0; i < nspec; ++i ) {
    for ( size_t j = 0; j < nlambda; ++j ) {
      double orig = data[ i * nlambda + j ];
      double check = check_data[ i * nlambda + j ];

      if ( fabs ( check - orig ) / orig > 1.0e-5 ) {
        cerr << "FAIL:  written spec " << i << ", bin " << j << " does not match original value" << endl;
        exit(1);
      }
    }
  }

  cerr << "  (PASSED)" << endl;
     
  return;
}



