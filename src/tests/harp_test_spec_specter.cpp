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


void harp::test_spec_specter ( string const & datadir ) {

  string filepath = datadir + "/test_specter_20120212.fits";
  string outpath = datadir + "/test_specter_20120212_check.fits.out";

  int statret;
  struct stat statbuf;
  statret = stat ( filepath.c_str(), &statbuf );
    
  if ( statret == 0 ) {

    cerr << "Testing spec_specter operations..." << endl;

    int ret = remove ( outpath.c_str() );

    boost::property_tree::ptree props;
    props.put ( "format", "specter" );
    props.put ( "path", filepath );

    spec_p testspec ( spec::create ( props ) );

    // immediately serialize and restore, so that any issues with that process will impact the code that follows

    string serialpath = datadir + "/test_specter_serialize.xml.out";
    {
      ofstream ofs ( serialpath.c_str() );
      boost::archive::binary_oarchive oa ( ofs );
      oa << testspec;
    }
    {
      ifstream ifs ( serialpath.c_str() );
      boost::archive::binary_iarchive ia ( ifs );
      ia >> testspec;
    }

    size_t nspec = testspec->n_spec();
    if ( nspec != 500 ) {
      cerr << "FAIL:  number of spectra (" << nspec << ") is not 500" << endl;
      exit(1);
    }

    size_t nlambda = testspec->n_lambda();
    if ( nlambda != 4697 ) {
      cerr << "FAIL:  number of wavelength points (" << nlambda << ") is not 4697" << endl;
      exit(1);
    }

    size_t nglobal = nspec * nlambda;

    vector_double specdata;
    vector_double lambda;
    vector < bool > sky;

    testspec->read ( specdata, lambda, sky );

    testspec->write ( outpath, specdata, lambda, sky );

    vector_double check_specdata;
    vector_double check_lambda;
    vector < bool > check_sky;

    props.clear();
    props.put ( "format", "specter" );
    props.put ( "path", outpath );
    spec_p checkspec ( spec::create ( props ) );

    if ( checkspec->n_spec() != nspec ) {
      cerr << "FAIL:  out file has " << checkspec->n_spec() << " spectra instead of " << nspec << endl;
      exit(1);
    }

    if ( checkspec->n_lambda() != nlambda ) {
      cerr << "FAIL:  out file has " << checkspec->n_lambda() << " wavelength points instead of " << nlambda << endl;
      exit(1);
    }

    checkspec->read ( check_specdata, check_lambda, check_sky );

    for ( size_t i = 0; i < sky.size(); ++i ) {
      if ( check_sky[i] != sky[i] ) {
        cerr << "FAIL:  sky element " << i << " does not match original value" << endl;
        exit(1);
      }
    }

    for ( size_t i = 0; i < lambda.size(); ++i ) {
      if ( fabs ( check_lambda[i] - lambda[i] ) / lambda[i] > 1.0e-5 ) {
        cerr << "FAIL:  lambda element " << i << " does not match original value" << endl;
        exit(1);
      }
    }

    for ( size_t i = 0; i < nspec; ++i ) {
      for ( size_t j = 0; j < nlambda; ++j ) {
        double orig = specdata[ i * nlambda + j ];
        double check = check_specdata[ i * nlambda + j ];

        if ( fabs ( check - orig ) / orig > 1.0e-5 ) {
          cerr << "FAIL:  spec " << i << ", bin " << j << " does not match original value" << endl;
          exit(1);
        }
      }
    }

    cerr << "  (PASSED)" << endl;

  } else {
    cerr << "Skipping spec_specter tests (file " << filepath << " not found)" << endl;
  }
     
  return;
}



