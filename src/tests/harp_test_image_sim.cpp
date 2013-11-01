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


void harp::test_image_sim ( string const & datadir ) {

  string specpath = datadir + "/spec_sim.fits.out";
  string psfpath = datadir + "/psf_gauss_sim.fits.out";
  string imgpath = datadir + "/image_sim.fits.out";

  cerr << "Testing simulated image..." << endl;

  // instantiate the spec

  boost::property_tree::ptree spec_props;
  spec_props.clear();
  spec_props.put ( "format", "specter" );
  spec_props.put ( "path", specpath );

  spec_p checkspec ( spec::create ( spec_props ) );

  size_t spec_nspec = checkspec->n_spec();
  size_t spec_nlambda = checkspec->n_lambda();

  vector_double check_data;
  vector_double check_lambda;
  vector < bool > check_sky;

  checkspec->read( check_data, check_lambda, check_sky );

  // instantiate the PSF

  boost::property_tree::ptree gauss_props;
  gauss_props.put ( "format", "gauss_sim" );
  gauss_props.put_child ( "lambda_spec", spec_props );
  gauss_props.put ( "bundle_size", 25 );
  gauss_props.put ( "nbundle", 4 );

  psf_p gauss_psf ( psf::create ( gauss_props ) );

  size_t psf_nspec = gauss_psf->n_spec();
  size_t psf_nlambda = gauss_psf->n_lambda();
  size_t psf_imgrows = gauss_psf->img_rows();
  size_t psf_imgcols = gauss_psf->img_cols();

  // instantiate image

  boost::property_tree::ptree img_props;
  img_props.put ( "format", "sim" );
  img_props.put_child ( "spec", spec_props );
  img_props.put_child ( "psf", gauss_props );

  cerr << "creating sim image" << endl;

  image_p img ( image::create ( img_props ) );

  // immediately serialize and restore, so that any issues with that process will impact the code that follows

  string serialpath = datadir + "/image_sim_serialize.xml.out";
  {
    ofstream ofs ( serialpath.c_str() );
    boost::archive::xml_oarchive oa ( ofs );
    oa << BOOST_SERIALIZATION_NVP(img);
  }
  {
    ifstream ifs ( serialpath.c_str() );
    boost::archive::xml_iarchive ia ( ifs );
    ia >> BOOST_SERIALIZATION_NVP(img);
  }

  size_t img_nrows = img->n_rows();
  size_t img_ncols = img->n_cols();

  // pixel dimensions should match the psf class

  if ( img_nrows != psf_imgrows ) {
    cerr << "FAIL:  simulated image rows (" << img_nrows << ") does not equal psf rows (" << psf_imgrows << ")" << endl;
    exit(1);
  }

  if ( img_ncols != psf_imgcols ) {
    cerr << "FAIL:  simulated image cols (" << img_ncols << ") does not equal psf cols (" << psf_imgcols << ")" << endl;
    exit(1);
  }

  // read / write test

  vector_double img_data;
  vector_double img_inv;
  vector < bool > sky;

  img->read ( img_data, img_inv, sky );

  cerr << "sim image len = " << img_data.size() << endl;

  for ( size_t i = 0; i < sky.size(); ++i ) {
    if ( sky[i] != check_sky[i] ) {
      cerr << "FAIL:  image sky element " << i << " does not match original value" << endl;
      exit(1);
    }
  }

  img->write ( imgpath, img_data, img_inv, sky );

  cerr << "  (PASSED)" << endl;

  return;
}



