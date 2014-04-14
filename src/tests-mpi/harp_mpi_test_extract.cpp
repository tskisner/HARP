#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include <harp_mpi_test.hpp>


using namespace std;
using namespace harp;



void mpi_test_extract_subaccum ( mpi_spec_slice_p slice, mpi_matrix & full_data, mpi_matrix & check_data ) {

  boost::mpi::communicator comm = slice->comm();

  int np = comm.size();
  int myp = comm.rank();

  spec_slice_region full_region = slice->full_region();

  std::vector < spec_slice_region > myslice = slice->regions();

  size_t indx = 0;

  mpi_matrix_zero ( check_data );

  for ( std::vector < spec_slice_region > :: const_iterator sit = myslice.begin(); sit != myslice.end(); ++sit ) {

    size_t sub_nbin = sit->n_spec * sit->n_lambda;

    mpi_matrix sub_data ( sub_nbin, 1 );
    mpi_matrix_zero ( sub_data );

    mpi_sub_spec ( full_region, (*sit), full_data, false, sub_data );

    for ( size_t i = 0; i < sub_nbin; ++i ) {
      size_t spec_offset = (size_t)( i / sit->n_lambda );
      size_t cur_spec = sit->first_spec + spec_offset;
      size_t cur_lambda = sit->first_lambda + i - ( spec_offset * sit->n_lambda );
      size_t cur_bin = cur_spec * full_region.n_lambda + cur_lambda;
      if ( (size_t)sub_data.Get(i,0) != cur_bin ) {
        cerr << "FAIL on sub_spec proc " << myp << ", region " << indx << ", bin " << i << ", " << (size_t)sub_data.Get(i,0) << " != " << cur_bin << endl;
        //exit(1);
      }
    }

    mpi_accum_spec ( (*sit), full_region, sub_data, true, check_data );

    for ( size_t i = 0; i < sub_nbin; ++i ) {
      size_t spec_offset = (size_t)( i / sit->n_lambda );
      size_t cur_spec = sit->first_spec + spec_offset;
      size_t cur_lambda = sit->first_lambda + i - ( spec_offset * sit->n_lambda );
      size_t cur_bin = cur_spec * full_region.n_lambda + cur_lambda;
      if ( ( cur_spec >= sit->first_good_spec ) && ( cur_spec < sit->first_good_spec + sit->n_good_spec ) && ( cur_lambda >= sit->first_good_lambda ) && ( cur_lambda < sit->first_good_lambda + sit->n_good_lambda ) ) {
        if ( (size_t)check_data.Get( cur_bin, 0 ) != cur_bin ) {
          cerr << "FAIL on accum_spec proc " << myp << ", region " << indx << ", bin " << i << ", " << (size_t)check_data.Get( cur_bin, 0 ) << " != " << cur_bin << endl;
          exit(1);
        }
      }
    }

    ++indx;

  }

  for ( size_t i = 0; i < full_data.Height(); ++i ) {
    if ( fabs ( full_data.Get(i, 0) > std::numeric_limits < double > :: epsilon() ) ) {
      if ( ( fabs ( ( check_data.Get(i, 0) - full_data.Get(i, 0) ) / full_data.Get(i, 0) ) ) > std::numeric_limits < double > :: epsilon() ) {
        cerr << "FAIL on spectral bin " << i << ", " << check_data.Get(i,0) << " != " << full_data.Get(i,0) << endl;
        exit(1);
      }
    }
  }

  return;
}




void harp::mpi_test_extract ( string const & datadir ) {

  boost::mpi::communicator comm;

  int np = comm.size();
  int myp = comm.rank();

  bool reg_mpi = true;
  plugin_registry & reg = plugin_registry::get( reg_mpi );

  if ( myp == 0 ) {
    cout << "Testing extraction spectral sub/accum functions..." << endl;
  }

  size_t nspec = 5;
  size_t nlambda = 200;
  size_t chunk_spec = 5;
  size_t overlap_spec = 0;
  size_t chunk_lambda = 20;
  size_t overlap_lambda = 10;

  size_t nbin = nspec * nlambda;

  mpi_matrix full_data ( nbin, 1 );
  mpi_matrix_zero ( full_data );

  mpi_matrix check_data ( nbin, 1 );

  for ( size_t i = 0; i < nbin; ++i ) {
    full_data.Set ( i, 0, (double)i );
  }
  
  mpi_spec_slice_p slice ( new mpi_spec_slice ( comm, nspec, nlambda, chunk_spec, chunk_lambda, overlap_spec, overlap_lambda ) );

  mpi_test_extract_subaccum ( slice, full_data, check_data );

  cout << "  (PASSED)" << endl;

  return;

  
  if ( myp == 0 ) {
    cout << "Testing gang-parallel slice and accum..." << endl;
  }

  // split communicator

  elem::Grid grid ( elem::mpi::COMM_WORLD );

  int gangsize = (int)( np / 2 );
  if ( gangsize < 1 ) {
    gangsize = 1;
  }

  int ngang = (int)( np / gangsize );
  int gangtot = ngang * gangsize;
  if ( myp == 0 ) {
    cout << "  Using " << ngang << " gangs of " << gangsize << " processes each" << endl;
  }
  if ( gangtot < np ) {
    if ( myp == 0 ) {
      cout << "  WARNING: " << (np-gangtot) << " processes are idle" << endl;
    }
  }
  int gang = (int)( myp / gangsize );
  int grank = myp % gangsize;
  if ( gang >= ngang ) {
    gang = MPI_UNDEFINED;
    grank = MPI_UNDEFINED;
  }

  boost::mpi::communicator gcomm = comm.split ( gang, grank );

  elem::Grid gang_grid ( gcomm );

  // setup data

  mpi_matrix_zero ( check_data );

  mpi_matrix gang_full_data ( nbin, 1, gang_grid );
  mpi_matrix gang_check_data ( nbin, 1, gang_grid );

  // world --> gang

  mpi_gang_distribute ( full_data, gang_full_data );

  // define slices

  mpi_spec_slice_p gang_slice ( new mpi_spec_slice ( gcomm, nspec, nlambda, chunk_spec, chunk_lambda, overlap_spec, overlap_lambda ) );

  // do sub / accum within each gang

  mpi_test_extract_subaccum ( gang_slice, gang_full_data, gang_check_data );

  // accum from all gangs

  mpi_gang_accum ( gang_check_data, check_data );

  for ( size_t i = 0; i < nbin; ++i ) {
    if ( ( fabs ( ( check_data.Get(i,0) - (double)ngang * full_data.Get(i,0) ) / (double)ngang * full_data.Get(i,0) ) ) > std::numeric_limits < double > :: epsilon() ) {
      cerr << "FAIL on spectral bin " << i << ", " << check_data.Get(i,0) << " != " << full_data.Get(i,0) << endl;
      exit(1);
    }
  }



  if ( myp == 0 ) {
    cout << "  (PASSED)" << endl;
  }


  /*

  cout << "Testing high-level, chunked extraction..." << endl;

  slice.reset ( new spec_slice ( 1, nspec, nlambda, chunk_spec, chunk_lambda, overlap_spec, overlap_lambda ) );

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
  vector < target > target_list;

  testspec->values ( truth );
  testspec->lambda ( lambda );
  testspec->targets ( target_list );

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

  // do extraction

  vector < spec_slice_region > regions = slice->regions ( 0 );

  string prefix = "  extract:  ";

  cout << prefix << "Extracting " << regions.size() << " spectral chunks, each with " << chunk_spec << " spectra (overlap = " << overlap_spec << ") and " << chunk_lambda << " lambda points (overlap = " << overlap_lambda << ")" << endl;

  map < string, double > timing;

  extract_slices ( slice, gauss_psf, img_data, img_inv, truth, Rf, f, err, Rtruth, timing, false, prefix );
     
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

  cout << "  (PASSED)" << endl;

  */

  return;
}



