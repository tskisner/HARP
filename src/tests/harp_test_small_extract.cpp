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


void harp::test_small_extract ( string const & datadir ) {

  int np;
  int myp;
  int ret;

  MPI_Comm_size ( MPI_COMM_WORLD, &np );
  MPI_Comm_rank ( MPI_COMM_WORLD, &myp );

  if ( myp == 0 ) {
    cerr << "Testing small spectral extraction..." << endl;
  }


  string psffile = datadir + "/test_small_psf.fits";
  string specfile = datadir + "/test_small_spec.fits";
  string imgfile = datadir + "/small_extract_image.fits.out";

  // create spec

  boost::property_tree::ptree spec_props;
  spec_props.clear();
  spec_props.put ( "format", "specter" );
  spec_props.put ( "path", specfile );
  spec_p testspec ( spec::create ( spec_props ) );

  // Create PSF

  boost::property_tree::ptree psf_props;
  psf_props.put ( "format", "gauss" );
  psf_props.put ( "path", psffile );
  psf_props.put ( "corr", 10 );
  psf_props.put ( "imgrows", 100 );
  psf_props.put ( "imgcols", 200 );
  psf_p testpsf ( psf::create ( psf_props ) );

  size_t nlambda = testpsf->nlambda();
  size_t nspec = testpsf->nspec();
  size_t nbins = nspec * nlambda;

  // Create image

  boost::property_tree::ptree img_props;
  img_props.put ( "format", "sim" );
  img_props.put_child ( "psf", psf_props );
  img_props.put_child ( "spec", spec_props );
  img_props.put ( "debug", imgfile );
  image_p testimg ( image::create ( img_props ) );

  // Read truth spec

  matrix_dist truth ( nbins, 1 );
  vector < double > lambda;
  vector < bool > sky;

  testspec->read( truth, lambda, sky );

  string outfile = datadir + "/small_extract_truth.out";
  elem::Write ( truth, outfile );

  // Read image and noise covariance

  size_t rows = testimg->rows();
  size_t cols = testimg->cols();
  size_t npix = rows * cols;

  matrix_local measured ( npix, 1 );  
  local_matrix_zero ( measured );

  matrix_local invnoise ( npix, 1 );
  local_matrix_zero ( invnoise );

  testimg->read ( measured );
  testimg->read_noise ( invnoise );

  // generate design matrix

  matrix_sparse design;

  double tstart;
  double tstop;

  tstart = MPI_Wtime();
  testpsf->projection ( 0, nspec - 1, 0, nlambda - 1, design );
  tstop = MPI_Wtime();
  if ( myp == 0 ) {
    cerr << "  Time for PSF creation = " << tstop-tstart << " seconds" << endl;
  }

  elem::Grid grid ( elem::mpi::COMM_WORLD );
  
  matrix_dist inv ( nbins, nbins, grid );

  tstart = MPI_Wtime();
  inverse_covariance ( design, invnoise, inv );
  tstop = MPI_Wtime();
  if ( myp == 0 ) {
    cerr << "  Time for inverse covariance building = " << tstop-tstart << " seconds" << endl;
  }

  matrix_dist W ( nbins, nbins, grid );
  matrix_dist D ( nbins, 1, grid );

  tstart = MPI_Wtime();
  eigen_decompose ( inv, D, W );
  tstop = MPI_Wtime();
  if ( myp == 0 ) {
    cerr << "  Time for eigendecomposition = " << tstop-tstart << " seconds" << endl;
  }

  outfile = datadir + "/small_extract_eigen.out";
  elem::Write ( D, outfile );

  matrix_dist S ( nbins, 1, grid );

  tstart = MPI_Wtime();
  norm ( D, W, S );
  tstop = MPI_Wtime();
  if ( myp == 0 ) {
    cerr << "  Time to compute S norm = " << tstop-tstart << " seconds" << endl;
  }

  outfile = datadir + "/small_extract_colnorm.out";
  elem::Write ( S, outfile );

  matrix_dist z ( nbins, 1 );
  matrix_dist Rf ( nbins, 1 );
  matrix_dist f ( nbins, 1 );

  tstart = MPI_Wtime();
  noise_weighted_spec ( design, invnoise, measured, z );
  tstop = MPI_Wtime();
  if ( myp == 0 ) {
    cerr << "  Time to compute RHS = " << tstop-tstart << " seconds" << endl;
  }

  outfile = datadir + "/small_extract_rhs.out";
  elem::Write ( z, outfile );

  tstart = MPI_Wtime();
  extract ( D, W, S, z, Rf, f );
  tstop = MPI_Wtime();
  if ( myp == 0 ) {
    cerr << "  Time for extraction = " << tstop-tstart << " seconds" << endl;
  }

  outfile = datadir + "/small_extract_Rf.out";
  elem::Write ( Rf, outfile );

  matrix_dist Rtruth ( nbins, 1 );
  matrix_dist R ( nbins, nbins, grid );

  tstart = MPI_Wtime();
  resolution ( D, W, S, R );
  tstop = MPI_Wtime();

  outfile = datadir + "/small_extract_colnorm-res.out";
  elem::Write ( S, outfile );

  if ( myp == 0 ) {
    cerr << "  Time for R matrix construction = " << tstop-tstart << " seconds" << endl;
  }

  // resolution-convolved truth

  elem::Gemv ( elem::NORMAL, 1.0, R, truth, 0.0, Rtruth );

  outfile = datadir + "/small_extract_Rtruth.out";
  elem::Write ( Rtruth, outfile );

  // do some sub spec and accum operations to test those...

  matrix_dist test_Rtruth ( nbins, 1 );
  dist_matrix_zero ( test_Rtruth );

  matrix_dist test_Rf ( nbins, 1 );
  dist_matrix_zero ( test_Rf );

  matrix_dist test_truth ( nbins, 1 );
  dist_matrix_zero ( test_truth );

  matrix_dist out_spec ( nspec * 20, 1 );

  for ( size_t b = 0; b < 3; ++b ) {

    sub_spec ( Rtruth, nspec, 0, nspec, 20 * b, 20, out_spec );
    accum_spec ( test_Rtruth, nspec, 0, nspec, 20 * b, 20, out_spec );

    sub_spec ( truth, nspec, 0, nspec, 20 * b, 20, out_spec );
    accum_spec ( test_truth, nspec, 0, nspec, 20 * b, 20, out_spec );

    sub_spec ( Rf, nspec, 0, nspec, 20 * b, 20, out_spec );
    accum_spec ( test_Rf, nspec, 0, nspec, 20 * b, 20, out_spec );

  }

  Rtruth = test_Rtruth;
  truth = test_truth;
  Rf = test_Rf;

  // write outputs

  fstream fout;
  fout.precision(16);

  matrix_local full;

  fitsfile * fp;

  if ( myp == 0 ) {
    outfile = datadir + "/small_extract.out";
    fout.open ( outfile.c_str(), ios::out );

    full.ResizeTo ( nbins, 1 );
    string outspec = datadir + "/small_extract.fits.out";
    fits::create ( fp, outspec );
  }

  for ( size_t i = 0; i < nbins; ++i ) {
    double dtemp = truth.Get(i,0);
    if ( myp == 0 ) {
      full.Set ( i, 0, dtemp );
    }
  }
  if ( myp == 0 ) {
    fits::img_append ( fp, nspec, nlambda );
    fits::img_write ( fp, full );
  }

  for ( size_t i = 0; i < nbins; ++i ) {
    double dtemp = Rtruth.Get(i,0);
    if ( myp == 0 ) {
      full.Set ( i, 0, dtemp );
    }
  }
  if ( myp == 0 ) {
    fits::img_append ( fp, nspec, nlambda );
    fits::img_write ( fp, full );
  }

  for ( size_t i = 0; i < nbins; ++i ) {
    double dtemp = Rf.Get(i,0);
    if ( myp == 0 ) {
      full.Set ( i, 0, dtemp );
    }
  }
  if ( myp == 0 ) {
    fits::img_append ( fp, nspec, nlambda );
    fits::img_write ( fp, full );
  }

  for ( size_t i = 0; i < nbins; ++i ) {
    double dtemp = inv.Get(i,i);
    if ( myp == 0 ) {
      full.Set ( i, 0, sqrt(1.0/dtemp) );
    }
  }
  if ( myp == 0 ) {
    fits::img_append ( fp, nspec, nlambda );
    fits::img_write ( fp, full );
  }

  for ( size_t i = 0; i < nbins; ++i ) {
    double out_rf = Rf.Get(i,0);
    double out_rt = Rtruth.Get(i,0);
    double out_tr = truth.Get(i,0);
    double errval = S.Get ( i, 0 );
    if ( myp == 0 ) {
      fout << i << " " << out_tr << " " << out_rt << " " << out_rf << " " << errval << endl;
    }
  }

  if ( myp == 0 ) {
    fout.close();
    fits::close ( fp );
  }

  if ( myp == 0 ) {
    cerr << "  (PASSED)" << endl;
  }
     
  return;
}



