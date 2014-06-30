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


void harp::test_spec_simspecter ( string const & datadir ) {

  plugin_registry & reg = plugin_registry::get();

  cout << "Testing simulated spec / targets generation..." << endl;

  size_t nlambda = 50;
  size_t nspec = 25;
  size_t skymod = 25;
  double first_lambda = 8000.0;
  double last_lambda = 8001.1;

  size_t nbins = nspec * nlambda;

  boost::property_tree::ptree spec_props;
  spec_props.clear();
  spec_props.put ( "nspec", nspec );
  spec_props.put ( "lambda_n", nlambda );
  spec_props.put ( "lambda_start", first_lambda );
  spec_props.put ( "lambda_stop", last_lambda );
  spec_props.put ( "back", 10.0 );
  spec_props.put ( "atm", 500.0 );
  spec_props.put ( "obj", 80.0 );
  spec_props.put ( "atmspace", 12 );
  spec_props.put ( "skymod", skymod );
  spec_p testspec ( reg.create_spec ( "sim", spec_props ) );

  // immediately serialize and restore, so that any issues with that process will impact the code that follows

  string serialpath = datadir + "/spec_sim_serialize.xml.out";
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

  testspec->values ( data );
  testspec->lambda ( lambda );

  vector_double invvar ( data );
  invvar.clear();

  // create targets

  boost::property_tree::ptree target_props;
  target_props.clear();
  target_props.put ( "nobj", nspec );
  target_props.put ( "skymod", skymod );
  targets_p testtargets ( reg.create_targets ( "sim", target_props ) );

  serialpath = datadir + "/targets_sim_serialize.xml.out";
  {
    ofstream ofs ( serialpath.c_str() );
    boost::archive::xml_oarchive oa ( ofs );
    oa << BOOST_SERIALIZATION_NVP(testtargets);
  }
  {
    ifstream ifs ( serialpath.c_str() );
    boost::archive::xml_iarchive ia ( ifs );
    ia >> BOOST_SERIALIZATION_NVP(testtargets);
  }

  if ( nspec != testtargets->n_objects() ) {
    cerr << "FAIL:  simulated targets n_objects (" << testtargets->n_objects() << ") does not match input (" << nspec << ")" << endl;
    exit(1);
  }

  vector < object_p > obj_list = testtargets->objects();

  for ( size_t i = 0; i < nspec; ++i ) {
    if ( i % skymod == 0 ) {
      if ( obj_list[i]->type() != OBJECT_SKY ) {
        cerr << "FAIL:  simulated targets object[" << i << "] not SKY" << endl;
        exit(1);
      }
    } else {
      if ( obj_list[i]->type() != OBJECT_UNKNOWN ) {
        cerr << "FAIL:  simulated targets object[" << i << "] not UNKNOWN" << endl;
        exit(1);
      }
    }
  }

  // we now write out the data to FITS format files

  string outfile = datadir + "/spec_sim.fits.out";

  spec_fits::write ( outfile, data, invvar, lambda );

  outfile = datadir + "/targets_sim.fits.out";

  targets_fits::write ( outfile, 2, obj_list );

  cout << "  (PASSED)" << endl;


  // now read it back in and check values

  cout << "Testing spec / target FITS operations..." << endl;

  outfile = datadir + "/spec_sim.fits.out";
  spec_props.clear();
  spec_props.put ( "path", outfile );

  spec_p checkspec ( reg.create_spec ( "fits", spec_props ) );

  serialpath = datadir + "/spec_fits_serialize.xml.out";
  {
    ofstream ofs ( serialpath.c_str() );
    boost::archive::xml_oarchive oa ( ofs );
    oa << BOOST_SERIALIZATION_NVP(checkspec);
  }
  {
    ifstream ifs ( serialpath.c_str() );
    boost::archive::xml_iarchive ia ( ifs );
    ia >> BOOST_SERIALIZATION_NVP(checkspec);
  }

  if ( nspec != checkspec->n_spec() ) {
    cerr << "FAIL:  fits written nspec (" << checkspec->n_spec() << ") does not match input (" << nspec << ")" << endl;
    exit(1);
  }

  if ( nlambda != checkspec->n_lambda() ) {
    cerr << "FAIL:  fits written nlambda (" << checkspec->n_lambda() << ") does not match input (" << nlambda << ")" << endl;
    exit(1);
  }

  vector_double check_data;
  vector_double check_lambda;

  checkspec->values ( check_data );
  checkspec->lambda ( check_lambda );

  for ( size_t i = 0; i < lambda.size(); ++i ) {
    if ( fabs ( ( check_lambda[i] - lambda[i] ) / lambda[i] ) > 1.0e-5 ) {
      cerr << "FAIL:  fits written lambda element " << i << " does not match original value" << endl;
      exit(1);
    }
  }

  for ( size_t i = 0; i < nspec; ++i ) {
    for ( size_t j = 0; j < nlambda; ++j ) {
      double orig = data[ i * nlambda + j ];
      double check = check_data[ i * nlambda + j ];

      if ( fabs ( ( check - orig ) / orig ) > 1.0e-5 ) {
        cerr << "FAIL:  fits written spec " << i << ", bin " << j << " does not match original value" << endl;
        exit(1);
      }
    }
  }

  outfile = datadir + "/targets_sim.fits.out";

  target_props.clear();
  target_props.put ( "path", outfile );
  target_props.put ( "hdu", "2" );
  targets_p checktargets ( reg.create_targets ( "fits", target_props ) );

  if ( nspec != checktargets->n_objects() ) {
    cerr << "FAIL:  fits targets n_objects (" << checktargets->n_objects() << ") does not match input (" << nspec << ")" << endl;
    exit(1);
  }

  obj_list = checktargets->objects();

  for ( size_t i = 0; i < nspec; ++i ) {
    if ( i % skymod == 0 ) {
      if ( obj_list[i]->type() != OBJECT_SKY ) {
        cerr << "FAIL:  fits targets object[" << i << "] not SKY" << endl;
        exit(1);
      }
    } else {
      if ( obj_list[i]->type() != OBJECT_UNKNOWN ) {
        cerr << "FAIL:  fits targets object[" << i << "] not UNKNOWN" << endl;
        exit(1);
      }
    }
  }

  cout << "  (PASSED)" << endl;
     
  return;
}



