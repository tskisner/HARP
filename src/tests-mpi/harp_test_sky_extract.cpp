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


void harp::test_sky_extract ( string const & datadir ) {

  int np;
  int myp;
  int ret;

  MPI_Comm_size ( MPI_COMM_WORLD, &np );
  MPI_Comm_rank ( MPI_COMM_WORLD, &myp );

  if ( myp == 0 ) {
    cerr << "Testing sky subtraction..." << endl;
  }

  size_t nlambda = 60;
  size_t nspec = 20;
  size_t first_lambda = 8000.0;
  size_t last_lambda = 8004.9;

  size_t nbins = nspec * nlambda;

  // create simulated spec

  boost::property_tree::ptree spec_props;
  spec_props.clear();
  spec_props.put ( "format", "sim" );
  spec_props.put ( "nspec", nspec );
  spec_props.put ( "nlambda", nlambda );
  spec_props.put ( "first_lambda", first_lambda );
  spec_props.put ( "last_lambda", last_lambda );
  spec_props.put ( "back", 10.0 );
  spec_props.put ( "atm", 5000.0 );
  spec_props.put ( "obj", 1000.0 );
  spec_props.put ( "atmspace", 12 );
  spec_props.put ( "skymod", 2 );
  spec_p testspec ( spec::create ( spec_props ) );

  // Create PSF

  boost::property_tree::ptree psf_props;
  psf_props.put ( "format", "gauss" );
  psf_props.put ( "corr", 10 );
  psf_props.put ( "fake", "TRUE" );
  psf_props.put_child ( "spec", spec_props );
  psf_props.put ( "bundle_size", 20 );
  psf_props.put ( "nbundle", 1 );
  psf_props.put ( "fwhm", 2.4 );
  psf_props.put ( "margin", 10 );
  psf_props.put ( "gap", 7 );
  psf_p testpsf ( psf::create ( psf_props ) );

  // Create image

  string imgfile = datadir + "/sky_extract_image.fits.out";

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

  string outfile = datadir + "/sky_extract_truth.out";
  El::Write ( truth, "truth", outfile );

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

  fstream fout;
  fout.precision(16);

  El::Grid grid ( El::mpi::COMM_WORLD );
  
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

  outfile = datadir + "/sky_extract_orig_eigen.out";
  El::Write ( D, "eigen", outfile );

  matrix_dist S ( nbins, 1, grid );

  tstart = MPI_Wtime();
  norm ( D, W, S );
  tstop = MPI_Wtime();
  if ( myp == 0 ) {
    cerr << "  Time to compute S norm = " << tstop-tstart << " seconds" << endl;
  }

  outfile = datadir + "/sky_extract_orig_colnorm.out";
  El::Write ( S, "colnorm", outfile );

  matrix_dist z ( nbins, 1 );
  matrix_dist Rf ( nbins, 1 );
  matrix_dist f ( nbins, 1 );

  tstart = MPI_Wtime();
  noise_weighted_spec ( design, invnoise, measured, z );
  tstop = MPI_Wtime();
  if ( myp == 0 ) {
    cerr << "  Time to compute RHS = " << tstop-tstart << " seconds" << endl;
  }

  outfile = datadir + "/sky_extract_orig_rhs.out";
  El::Write ( z, "rhs", outfile );

  tstart = MPI_Wtime();
  extract ( D, W, S, z, Rf, f );
  tstop = MPI_Wtime();
  if ( myp == 0 ) {
    cerr << "  Time for extraction = " << tstop-tstart << " seconds" << endl;
  }

  outfile = datadir + "/sky_extract_orig_Rf.out";
  El::Write ( Rf, "Rf", outfile );

  matrix_dist Rtruth ( nbins, 1 );
  matrix_dist R ( nbins, nbins, grid );

  tstart = MPI_Wtime();
  resolution ( D, W, S, R );
  tstop = MPI_Wtime();
  if ( myp == 0 ) {
    cerr << "  Time for R matrix construction = " << tstop-tstart << " seconds" << endl;
  }

  // resolution-convolved truth

  El::Gemv ( El::NORMAL, 1.0, R, truth, 0.0, Rtruth );

  outfile = datadir + "/sky_extract_orig_Rtruth.out";
  El::Write ( Rtruth, "Rtruth", outfile );


  // write outputs

  matrix_local full;

  fitsfile * fp;

  if ( myp == 0 ) {
    outfile = datadir + "/sky_extract_orig.out";
    fout.open ( outfile.c_str(), ios::out );

    full.ResizeTo ( nbins, 1 );
    string outspec = datadir + "/sky_extract_orig.fits.out";
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


  // get sky truth

  spec_sim * testsimspec = dynamic_cast < spec_sim * > ( testspec.get() );

  matrix_dist truth_sky;

  testsimspec->sky_truth ( truth_sky );


  // generate sky subtraction design matrix

  size_t nspec_obj = 0;

  vector < bool > :: const_iterator itsky;
  for ( itsky = sky.begin(); itsky != sky.end(); ++itsky ) {
    if ( ! (*itsky) ) {
      ++nspec_obj;
    }
  }

  size_t nspec_sky = nspec_obj + 1;
  size_t nbins_sky = nspec_sky * nlambda;

  cerr << "  " << nspec_obj << " object spectra plus one sky" << endl;
  cerr << "   (" << nbins_sky << ") bins" << endl; 

  matrix_sparse sdesign;
  sky_design ( design, sky, sdesign );


  outfile = datadir + "/sky_design_orig.out";
  fout.open ( outfile.c_str(), ios::out );
  cliq::Print ( design, "original", fout );
  fout.close();

  outfile = datadir + "/sky_design_new.out";
  fout.open ( outfile.c_str(), ios::out );
  cliq::Print ( sdesign, "sky", fout );
  fout.close();
  
  matrix_dist inv_sky ( nbins_sky, nbins_sky, grid );

  tstart = MPI_Wtime();
  inverse_covariance ( sdesign, invnoise, inv_sky );
  tstop = MPI_Wtime();
  if ( myp == 0 ) {
    cerr << "  Time for inverse covariance building = " << tstop-tstart << " seconds" << endl;
  }

  matrix_dist W_sky ( nbins_sky, nbins_sky, grid );
  matrix_dist D_sky ( nbins_sky, 1, grid );

  tstart = MPI_Wtime();
  eigen_decompose ( inv_sky, D_sky, W_sky );
  tstop = MPI_Wtime();
  if ( myp == 0 ) {
    cerr << "  Time for eigendecomposition = " << tstop-tstart << " seconds" << endl;
  }

  outfile = datadir + "/sky_extract_eigen.out";
  El::Write ( D_sky, "eigen", outfile );

  matrix_dist S_sky ( nbins_sky, 1, grid );

  tstart = MPI_Wtime();
  norm ( D_sky, W_sky, S_sky );
  tstop = MPI_Wtime();
  if ( myp == 0 ) {
    cerr << "  Time to compute S norm = " << tstop-tstart << " seconds" << endl;
  }

  outfile = datadir + "/sky_extract_colnorm.out";
  El::Write ( S_sky, "colnorm", outfile );

  matrix_dist z_sky ( nbins_sky, 1 );
  matrix_dist Rf_sky ( nbins_sky, 1 );
  matrix_dist f_sky ( nbins_sky, 1 );

  tstart = MPI_Wtime();
  noise_weighted_spec ( sdesign, invnoise, measured, z_sky );
  tstop = MPI_Wtime();
  if ( myp == 0 ) {
    cerr << "  Time to compute RHS = " << tstop-tstart << " seconds" << endl;
  }

  outfile = datadir + "/sky_extract_rhs.out";
  El::Write ( z_sky, "rhs", outfile );

  tstart = MPI_Wtime();
  extract ( D_sky, W_sky, S_sky, z_sky, Rf_sky, f_sky );
  tstop = MPI_Wtime();
  if ( myp == 0 ) {
    cerr << "  Time for extraction = " << tstop-tstart << " seconds" << endl;
  }

  outfile = datadir + "/sky_extract_Rf.out";
  El::Write ( Rf_sky, "Rf", outfile );

  matrix_dist R_sky ( nbins_sky, nbins_sky, grid );

  tstart = MPI_Wtime();
  resolution ( D_sky, W_sky, S_sky, R_sky );
  tstop = MPI_Wtime();
  if ( myp == 0 ) {
    cerr << "  Time for R matrix construction = " << tstop-tstart << " seconds" << endl;
  }

  // resolution-convolved truth

  matrix_dist Rtruth_sky ( nbins_sky, 1, grid );

  El::Gemv ( El::NORMAL, 1.0, R_sky, truth_sky, 0.0, Rtruth_sky );

  // projected truth

  matrix_local projected ( npix, 1 );
  spec_project ( sdesign, truth_sky, projected );

  if ( myp == 0 ) {
    string outimg = datadir + "/sky_extract_truth-project.fits.out";
    fits::create ( fp, outimg );
    fits::img_append ( fp, rows, cols );
    fits::write_key ( fp, "EXTNAME", "truth", "sky subtraction truth" );
    fits::img_write ( fp, projected );
    fits::close ( fp );
  }


  // write outputs

  if ( myp == 0 ) {
    outfile = datadir + "/sky_extract.out";
    fout.open ( outfile.c_str(), ios::out );

    full.ResizeTo ( nbins_sky, 1 );
    string outspec = datadir + "/sky_extract.fits.out";
    fits::create ( fp, outspec );
  }

  for ( size_t i = 0; i < nbins_sky; ++i ) {
    double dtemp = Rf_sky.Get(i,0);
    if ( myp == 0 ) {
      full.Set ( i, 0, dtemp );
    }
  }
  if ( myp == 0 ) {
    fits::img_append ( fp, nspec_sky, nlambda );
    fits::img_write ( fp, full );
  }

  for ( size_t i = 0; i < nbins_sky; ++i ) {
    double dtemp = inv_sky.Get(i,i);
    if ( myp == 0 ) {
      full.Set ( i, 0, sqrt(1.0/dtemp) );
    }
  }
  if ( myp == 0 ) {
    fits::img_append ( fp, nspec_sky, nlambda );
    fits::img_write ( fp, full );
  }

  for ( size_t i = 0; i < nbins_sky; ++i ) {
    double out_tr = truth_sky.Get(i,0);
    double out_rt = Rtruth_sky.Get(i,0);
    double out_rf = Rf_sky.Get(i,0);
    double errval = S_sky.Get ( i, 0 );
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



