#include <iostream>
#include <fstream>
#include <sstream>

#include <harp_test.hpp>

extern "C" {
  #include <unistd.h>
  #include <sys/stat.h>
}

using namespace std;
using namespace harp;


void harp::test_psf_gauss ( string const & datadir ) {

  plugin_registry & reg = plugin_registry::get();

  string specpath = datadir + "/spec_sim.fits.out";
  string inpath = datadir + "/psf_gauss_sim.fits.out";
  string outpath = datadir + "/psf_gauss.fits.out";

  cout << "Testing elliptical gaussian PSF..." << endl;

  // instantiate a spec to check the values against

  boost::property_tree::ptree spec_props;
  spec_props.clear();
  spec_props.put ( "path", specpath );

  spec_p checkspec ( reg.create_spec ( "specter", spec_props ) );

  size_t spec_nspec = checkspec->n_spec();
  size_t spec_nlambda = checkspec->n_lambda();

  vector_double check_data;
  vector_double check_lambda;
  vector < target > check_target_list;

  checkspec->values ( check_data );
  checkspec->lambda ( check_lambda );
  checkspec->targets ( check_target_list );

  // create the psf

  boost::property_tree::ptree gauss_props;
  gauss_props.put ( "path", inpath );
  gauss_props.put ( "corr", 10 );
  gauss_props.put ( "imgrows", 70 );
  gauss_props.put ( "imgcols", 853 );

  psf_gauss gauss_psf ( gauss_props );

  // immediately serialize and restore, so that any issues with that process will impact the code that follows

  string serialpath = datadir + "/psf_gauss_serialize.xml.out";
  {
    ofstream ofs ( serialpath.c_str() );
    boost::archive::xml_oarchive oa ( ofs );
    oa << BOOST_SERIALIZATION_NVP(gauss_psf);
  }
  {
    ifstream ifs ( serialpath.c_str() );
    boost::archive::xml_iarchive ia ( ifs );
    ia >> BOOST_SERIALIZATION_NVP(gauss_psf);
  }

  size_t psf_nspec = gauss_psf.n_spec();
  size_t psf_nlambda = gauss_psf.n_lambda();

  if ( psf_nspec != 25 ) {
    cerr << "FAIL:  gauss psf nspec (" << psf_nspec << ") does not match input (" << 25 << ")" << endl;
    exit(1);
  }

  if ( psf_nlambda != 50 ) {
    cerr << "FAIL:  gauss psf nlambda (" << psf_nlambda << ") does not match input (" << 50 << ")" << endl;
    exit(1);
  }

  size_t gauss_rows = gauss_psf.img_rows();
  size_t gauss_cols = gauss_psf.img_cols();

  // pixel dimensions should match the default spectral spacing in the psf class

  if ( gauss_cols != 853 ) {
    cerr << "FAIL:  gauss psf image cols (" << gauss_cols << ") is not 853" << endl;
    exit(1);
  }

  if ( gauss_rows != 70 ) {
    cerr << "FAIL:  gauss psf image cols (" << gauss_rows << ") is not 70" << endl;
    exit(1);
  }

  // check wavelength solution

  vector_double gauss_lambda = gauss_psf.lambda();

  for ( size_t i = 0; i < gauss_lambda.size(); ++i ) {
    if ( fabs ( ( check_lambda[i] - gauss_lambda[i] ) / check_lambda[i] ) > 1.0e-5 ) {
      cerr << "FAIL:  PSF lambda element " << i << " does not match original value" << endl;
      exit(1);
    }
  }

  // write out

  gauss_psf.write ( outpath );

  cout << "  (PASSED)" << endl;

  return;
}


