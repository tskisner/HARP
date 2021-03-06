#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include <harp_test.hpp>


using namespace std;
using namespace harp;


void harp::test_extract ( string const & datadir ) {

  plugin_registry & reg = plugin_registry::get();

  cout << "Testing extraction spectral sub/accum functions..." << endl;

  size_t nspec = 5;
  size_t nlambda = 200;
  size_t chunk_spec = 5;
  size_t overlap_spec = 0;
  size_t chunk_lambda = 20;
  size_t overlap_lambda = 10;

  size_t nbin = nspec * nlambda;

  vector_double full_data ( nbin );
  full_data.clear();

  vector_double check_data ( nbin );
  check_data.clear();

  for ( size_t i = 0; i < nbin; ++i ) {
    full_data[i] = (double)i;
  }
  
  size_t procs = 10;

  spec_slice_p slice ( new spec_slice ( procs, 0, overlap_lambda, nspec, nlambda-2*overlap_lambda, chunk_spec, chunk_lambda, overlap_spec, overlap_lambda ) );

  spec_slice_region full_region = slice->full_region();

  for ( size_t p = 0; p < procs; ++p ) {
    std::vector < spec_slice_region > procslice = slice->regions ( p );

    size_t indx = 0;

    for ( std::vector < spec_slice_region > :: const_iterator sit = procslice.begin(); sit != procslice.end(); ++sit ) {

      size_t sub_nbin = sit->n_spec * sit->n_lambda;

      vector_double sub_data ( sub_nbin );
      sub_data.clear();

      sub_spec ( full_region, (*sit), full_data, false, sub_data );

      for ( size_t i = 0; i < sub_nbin; ++i ) {
        size_t spec_offset = (size_t)( i / sit->n_lambda );
        size_t cur_spec = sit->first_spec + spec_offset;
        size_t cur_lambda = sit->first_lambda + i - ( spec_offset * sit->n_lambda );
        size_t cur_bin = cur_spec * nlambda + cur_lambda;
        if ( (size_t)sub_data[i] != cur_bin ) {
          cerr << "FAIL on sub_spec proc " << p << ", region " << indx << ", bin " << i << ", " << (size_t)sub_data[i] << " != " << cur_bin << endl;
          exit(1);
        }
      }

      accum_spec ( (*sit), full_region, sub_data, true, check_data );

      for ( size_t i = 0; i < sub_nbin; ++i ) {
        size_t spec_offset = (size_t)( i / sit->n_lambda );
        size_t cur_spec = sit->first_spec + spec_offset;
        size_t cur_lambda = sit->first_lambda + i - ( spec_offset * sit->n_lambda );
        size_t cur_bin = cur_spec * nlambda + cur_lambda;
        if ( ( cur_spec >= sit->first_good_spec ) && ( cur_spec < sit->first_good_spec + sit->n_good_spec ) && ( cur_lambda >= sit->first_good_lambda ) && ( cur_lambda < sit->first_good_lambda + sit->n_good_lambda ) ) {
          if ( (size_t)check_data[ cur_bin ] != cur_bin ) {
            cerr << "FAIL on accum_spec proc " << p << ", region " << indx << ", bin " << i << ", " << (size_t)check_data[ cur_bin ] << " != " << cur_bin << endl;
            exit(1);
          }
        }
      }

      ++indx;

    }
    
  }

  for ( size_t i = 0; i < nbin; ++i ) {
    size_t cur_spec = (size_t)( i / nlambda );
    size_t cur_lambda = i - (cur_spec * nlambda);
    if ( ( cur_spec >= overlap_spec ) && ( cur_spec < (nspec - overlap_spec) ) && ( cur_lambda >= overlap_lambda ) && ( cur_lambda < (nlambda - overlap_lambda) ) ) {
      if ( ( fabs ( ( check_data[i] - full_data[i] ) / full_data[i] ) ) > std::numeric_limits < double > :: epsilon() ) {
        cerr << "FAIL on spectral bin " << i << ", " << check_data[i] << " != " << full_data[i] << endl;
        exit(1);
      }
    }
  }

  full_data.resize(0);
  check_data.resize(0);

  cout << "  (PASSED)" << endl;

  cout << "Testing high-level, chunked extraction..." << endl;

  slice.reset ( new spec_slice ( 1, 0, overlap_lambda, nspec, nlambda-2*overlap_lambda, chunk_spec, chunk_lambda, overlap_spec, overlap_lambda ) );

  vector_double truth ( nbin );
  vector_double Rtruth ( nbin );
  vector_double f ( nbin );
  vector_double Rf ( nbin );
  vector_double err ( nbin );

  double first_lambda = 8000.0;
  double last_lambda = 8200.0;

  // create spec and read

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
  spec_props.put ( "skymod", nspec );
  spec_p testspec ( reg.create_spec ( "sim", spec_props ) );

  vector_double lambda;

  testspec->values ( truth );
  testspec->lambda ( lambda );
  vector_double fakeinvvar ( truth );
  fakeinvvar.clear();

  // instantiate the PSF

  boost::property_tree::ptree gauss_props;
  gauss_props.put ( "lambda_spec_type", "sim" );
  gauss_props.put_child ( "lambda_spec", spec_props );
  gauss_props.put ( "bundle_size", nspec );
  gauss_props.put ( "nbundle", 1 );

  psf_p gauss_psf ( reg.create_psf ( "gauss_sim", gauss_props ) );

  // instantiate image and read

  boost::property_tree::ptree img_props;
  img_props.put ( "spec_type", "sim" );
  img_props.put_child ( "spec", spec_props );
  img_props.put ( "psf_type", "gauss_sim" );
  img_props.put_child ( "psf", gauss_props );

  image_p img ( reg.create_image ( "sim", img_props ) );

  vector_double img_data;
  vector_double img_inv;

  img->values ( img_data );
  img->inv_variance ( img_inv );

  string outfile = datadir + "/extract_image_input.fits.out";
  image_fits::write ( outfile, img->n_rows(), img_data, img_inv );

  outfile = datadir + "/extract_spec_truth.fits.out";
  spec_fits::write ( outfile, truth, fakeinvvar, lambda );

  // do extraction

  vector < spec_slice_region > regions = slice->regions ( 0 );

  string prefix = "  extract:  ";

  cout << prefix << "Extracting " << regions.size() << " spectral chunks, each with " << chunk_spec << " spectra (overlap = " << overlap_spec << ") and " << chunk_lambda << " lambda points (overlap = " << overlap_lambda << ")" << endl;

  map < string, double > timing;

  bool lambda_mask = true;
  //bool lambda_mask = false;

  extract_slices ( slice, gauss_psf, img_data, img_inv, truth, Rf, f, err, Rtruth, timing, false, lambda_mask, prefix );
     
  cout << prefix << "Aggregate Timings:" << endl;
  cout << prefix << "  Build design matrix = " << timing["design"] << " seconds" << endl;
  cout << prefix << "  Build inverse covariance = " << timing["inverse"] << " seconds" << endl;
  cout << prefix << "  Eigendecompose inverse = " << timing["eigen"] << " seconds" << endl;
  cout << prefix << "  Compute column norm = " << timing["norm"] << " seconds" << endl;
  cout << prefix << "  Compute noise weighted spec = " << timing["nsespec"] << " seconds" << endl;
  cout << prefix << "  Extract spectra = " << timing["extract"] << " seconds" << endl;

  double chisq_reduced = 0.0;

  for ( size_t i = 0; i < nbin; ++i ) {
    if ( err[i] > std::numeric_limits < double > :: epsilon() ) {
      double val = ( Rf[i] - Rtruth[i] ) / err[i];
      chisq_reduced += val * val;
    }
  }

  chisq_reduced /= (double)( nbin );

  cout << prefix << "Reduced Chi square = " << chisq_reduced << endl;

  outfile = datadir + "/extract_spec_Rf.fits.out";
  spec_fits::write ( outfile, Rf, err, lambda );

  outfile = datadir + "/extract_spec_Rtruth.fits.out";
  spec_fits::write ( outfile, Rtruth, fakeinvvar, lambda );

  outfile = datadir + "/extract_spec_f.fits.out";
  spec_fits::write ( outfile, f, fakeinvvar, lambda );

  matrix_double_sparse AT;
  gauss_psf->project_transpose ( AT );

  vector_double f_projected;
  sparse_mv_trans ( AT, f, f_projected );

  outfile = datadir + "/extract_image_f-project.fits.out";
  image_fits::write ( outfile, img->n_rows(), f_projected, img_inv );



  // now do this for a fake flat

  spec_props.clear();
  spec_props.put ( "nspec", nspec );
  spec_props.put ( "lambda_n", nlambda );
  spec_props.put ( "lambda_start", first_lambda );
  spec_props.put ( "lambda_stop", last_lambda );
  spec_props.put ( "back", 500.0 );
  spec_props.put ( "atm", 0.0 );
  spec_props.put ( "obj", 0.0 );
  spec_props.put ( "atmspace", 12 );
  spec_props.put ( "skymod", nspec );
  testspec.reset ( reg.create_spec ( "sim", spec_props ) );

  testspec->values ( truth );
  testspec->lambda ( lambda );

  // instantiate image and read

  img_props.clear();
  img_props.put ( "spec_type", "sim" );
  img_props.put_child ( "spec", spec_props );
  img_props.put ( "psf_type", "gauss_sim" );
  img_props.put_child ( "psf", gauss_props );

  img.reset ( reg.create_image ( "sim", img_props ) );

  img->values ( img_data );
  img->inv_variance ( img_inv );

  outfile = datadir + "/flat_image_input.fits.out";
  image_fits::write ( outfile, img->n_rows(), img_data, img_inv );

  outfile = datadir + "/flat_spec_truth.fits.out";
  spec_fits::write ( outfile, truth, fakeinvvar, lambda );


  // do extraction

  regions = slice->regions ( 0 );

  prefix = "  flat:  ";

  cout << prefix << "Extracting " << regions.size() << " spectral chunks, each with " << chunk_spec << " spectra (overlap = " << overlap_spec << ") and " << chunk_lambda << " lambda points (overlap = " << overlap_lambda << ")" << endl;

  timing.clear();

  extract_slices ( slice, gauss_psf, img_data, img_inv, truth, Rf, f, err, Rtruth, timing, false, lambda_mask, prefix );
     
  cout << prefix << "Aggregate Timings:" << endl;
  cout << prefix << "  Build design matrix = " << timing["design"] << " seconds" << endl;
  cout << prefix << "  Build inverse covariance = " << timing["inverse"] << " seconds" << endl;
  cout << prefix << "  Eigendecompose inverse = " << timing["eigen"] << " seconds" << endl;
  cout << prefix << "  Compute column norm = " << timing["norm"] << " seconds" << endl;
  cout << prefix << "  Compute noise weighted spec = " << timing["nsespec"] << " seconds" << endl;
  cout << prefix << "  Extract spectra = " << timing["extract"] << " seconds" << endl;

  chisq_reduced = 0.0;

  for ( size_t i = 0; i < nbin; ++i ) {
    if ( err[i] > std::numeric_limits < double > :: epsilon() ) {
      double val = ( Rf[i] - Rtruth[i] ) / err[i];
      chisq_reduced += val * val;
    }
  }

  chisq_reduced /= (double)( nbin );

  cout << prefix << "Reduced Chi square = " << chisq_reduced << endl;

  outfile = datadir + "/flat_spec_Rf.fits.out";
  spec_fits::write ( outfile, Rf, err, lambda );

  outfile = datadir + "/flat_spec_Rtruth.fits.out";
  spec_fits::write ( outfile, Rtruth, fakeinvvar, lambda );

  outfile = datadir + "/flat_spec_f.fits.out";
  spec_fits::write ( outfile, f, fakeinvvar, lambda );


  gauss_psf->project_transpose ( AT );

  sparse_mv_trans ( AT, f, f_projected );

  outfile = datadir + "/flat_image_f-project.fits.out";
  image_fits::write ( outfile, img->n_rows(), f_projected, img_inv );


  cout << "  (PASSED)" << endl;

  return;
}



