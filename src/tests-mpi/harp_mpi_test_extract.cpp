#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include <harp_mpi_test.hpp>


using namespace std;
using namespace harp;


void harp::mpi_test_extract ( string const & datadir ) {

  bool reg_mpi = true;
  plugin_registry & reg = plugin_registry::get( reg_mpi );

  cout << "Testing extraction spectral sub/accum functions..." << endl;

  size_t nspec = 20;
  size_t nlambda = 50;
  size_t chunk_spec = 4;
  size_t overlap_spec = 1;
  size_t chunk_lambda = 10;
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

  spec_slice_p slice ( new spec_slice ( procs, nspec, nlambda, chunk_spec, chunk_lambda, overlap_spec, overlap_lambda ) );

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
    if ( ( fabs ( ( check_data[i] - full_data[i] ) / full_data[i] ) ) > std::numeric_limits < double > :: epsilon() ) {
      cerr << "FAIL on spectral bin " << i << ", " << check_data[i] << " != " << full_data[i] << endl;
      exit(1);
    }
  }

  full_data.resize(0);
  check_data.resize(0);

  cout << "  (PASSED)" << endl;


  
  if ( myp == 0 ) {
    cout << "Testing gang-parallel slice and accum..." << endl;
  }

  mpi_matrix fullspec ( NSPEC * SPECSIZE, 1, grid );
  mpi_matrix comp_fullspec ( NSPEC * SPECSIZE, 1, grid );

  for ( size_t i = 0; i < NSPEC; ++i ) {
    for ( size_t j = 0; j < SPECSIZE; ++j ) {
      fullspec.Set ( i * SPECSIZE + j, 0, 100.0 + (double)j );
    }
  }

  mpi_matrix gang_fullspec ( NSPEC * SPECSIZE, 1, gang_grid );

  mpi_gang_distribute ( fullspec, gang_fullspec );

  size_t nspec_chunk = (size_t)( NSPEC / 2 );
  size_t first_spec = (size_t)( NSPEC / 4 );
  size_t nlambda_chunk = (size_t)( SPECSIZE / 4 );
  size_t first_lambda = (size_t)( SPECSIZE / 8 );

  if ( myp == 0 ) {
    cout << "spec range = " << first_spec << " - " << (first_spec + nspec_chunk - 1) << endl;
    cout << "lambda range = " << first_lambda << " - " << (first_lambda + nlambda_chunk - 1) << endl;
  }

  mpi_matrix gang_subspec ( nspec_chunk * nlambda_chunk, 1, gang_grid );

  mpi_sub_spec ( gang_fullspec, NSPEC, first_spec, nspec_chunk, first_lambda, nlambda_chunk, gang_subspec );

  MPI_Barrier ( MPI_COMM_WORLD );

  mpi_accum_spec ( gang_fullspec, NSPEC, first_spec, nspec_chunk, first_lambda, nlambda_chunk, gang_subspec );

  mpi_gang_accum ( gang_fullspec, comp_fullspec );

  MPI_Barrier ( MPI_COMM_WORLD );

  for ( size_t i = 0; i < NSPEC; ++i ) {
    for ( size_t j = 0; j < SPECSIZE; ++j ) {
      inval = fullspec.Get ( i * SPECSIZE + j, 0 ) * (double)ngang;
      if ( ( ( i >= first_spec ) && ( i < first_spec + nspec_chunk ) ) && ( ( j >= first_lambda ) && ( j < first_lambda + nlambda_chunk ) ) ) {
        inval *= 2.0;
      }
      outval = comp_fullspec.Get ( i * SPECSIZE + j, 0 );
      relerr = fabs ( outval - inval ) / inval;
      if ( relerr > TOL ) {
        cerr << "FAIL on spectrum " << i << ", wavelength " << j << ": " << outval << " != " << inval << endl;
        exit(1);
      }

    }
  }


  if ( myp == 0 ) {
    cout << "  (PASSED)" << endl;
  }


  

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

  return;
}



