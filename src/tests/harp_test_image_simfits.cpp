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


void harp::test_image_simfits ( string const & datadir ) {

  plugin_registry & reg = plugin_registry::get();

  string specpath = datadir + "/spec_sim.fits.out";
  string psfpath = datadir + "/psf_gauss_sim.fits.out";
  string imgpath = datadir + "/image_sim.fits.out";

  cout << "Testing simulated image..." << endl;

  // instantiate the spec

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

  // instantiate the PSF

  boost::property_tree::ptree gauss_props;
  gauss_props.put ( "lambda_spec_type", "specter" );
  gauss_props.put_child ( "lambda_spec", spec_props );
  gauss_props.put ( "bundle_size", 25 );
  gauss_props.put ( "nbundle", 1 );

  psf_p gauss_psf ( reg.create_psf ( "gauss_sim", gauss_props ) );

  size_t psf_nspec = gauss_psf->n_spec();
  size_t psf_nlambda = gauss_psf->n_lambda();
  size_t psf_imgrows = gauss_psf->img_rows();
  size_t psf_imgcols = gauss_psf->img_cols();

  // instantiate image

  boost::property_tree::ptree img_props;
  img_props.put ( "spec_type", "specter" );
  img_props.put ( "psf_type", "gauss_sim" );
  img_props.put_child ( "spec", spec_props );
  img_props.put_child ( "psf", gauss_props );

  image_p img ( reg.create_image ( "sim", img_props ) );

  // write out a json file for testing with external tools

  boost::property_tree::ptree full_props;
  full_props.put_child ( "psf", gauss_props );
  full_props.put_child ( "image", img_props );
  full_props.put_child ( "truth", spec_props );

  string jsonpath = datadir + "/test_sim.json.out";
  boost::property_tree::json_parser::write_json ( jsonpath, full_props );

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

  // read test

  vector_double img_data;
  vector_double img_inv;

  img->values ( img_data );
  img->inv_variance ( img_inv );

  cout << "  (PASSED)" << endl;


  cout << "Testing FITS image..." << endl;

  img_props.clear();
  img_props.put ( "rows", psf_imgrows );
  img_props.put ( "cols", psf_imgcols );

  image_fits outimg ( img_props );

  // immediately serialize and restore, so that any issues with that process will impact the code that follows

  serialpath = datadir + "/image_fits_serialize.xml.out";
  {
    ofstream ofs ( serialpath.c_str() );
    boost::archive::xml_oarchive oa ( ofs );
    oa << BOOST_SERIALIZATION_NVP(outimg);
  }
  {
    ifstream ifs ( serialpath.c_str() );
    boost::archive::xml_iarchive ia ( ifs );
    ia >> BOOST_SERIALIZATION_NVP(outimg);
  }

  img_nrows = outimg.n_rows();
  img_ncols = outimg.n_cols();

  // pixel dimensions should match the psf class

  if ( img_nrows != psf_imgrows ) {
    cerr << "FAIL:  FITS image rows (" << img_nrows << ") does not equal psf rows (" << psf_imgrows << ")" << endl;
    exit(1);
  }

  if ( img_ncols != psf_imgcols ) {
    cerr << "FAIL:  FITS image cols (" << img_ncols << ") does not equal psf cols (" << psf_imgcols << ")" << endl;
    exit(1);
  }

  // write test

  outimg.write ( imgpath, img_data, img_inv );

  // read test

  img_props.clear();
  img_props.put ( "path", imgpath );

  image_p checkimg ( reg.create_image ( "fits", img_props ) );

  vector_double check_img_data;
  vector_double check_img_inv;

  checkimg->values ( check_img_data );
  checkimg->inv_variance ( check_img_inv );

  cout << "  (PASSED)" << endl;

  return;
}



