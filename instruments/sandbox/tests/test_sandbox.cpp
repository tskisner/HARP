#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include <harp_test.hpp>

#include <boost/random.hpp>


extern "C" {
#include <unistd.h>
}

using namespace std;
using namespace harp;


void sandbox_profile ( string const & name, string const & desc, double & totaltime, double & opencltime, map < string, long long int > & papi ) {
  
  cerr << "Profiling:   " << desc << ":  " << totaltime << " seconds" << endl;
  
  return;
}


void harp::test_sandbox ( string const & datadir ) {

  int np;
  int myp;
  int ret;

  MPI_Comm_size ( MPI_COMM_WORLD, &np );
  MPI_Comm_rank ( MPI_COMM_WORLD, &myp );

  string filepath = datadir + "/psf_g2d_2012-02-11.fits";

  int statret;
  if ( myp == 0 ) {
    struct stat statbuf;
    statret = stat ( filepath.c_str(), &statbuf );
  }

  ret = MPI_Bcast ( (void*)(&statret), 1, MPI_INT, 0, MPI_COMM_WORLD );
  mpi_check ( MPI_COMM_WORLD, ret );
    
  if ( statret == 0 ) {
    if ( myp == 0 ) {
      cerr << "Testing sandbox elliptical gaussian PSF..." << endl;
    }

    boost::property_tree::ptree gauss_props;
    gauss_props.put ( "format", "sandbox" );
    gauss_props.put ( "path", filepath );
    gauss_props.put ( "corr", 7 );
    gauss_props.put ( "imgrows", 4697 );
    gauss_props.put ( "imgcols", 4110 );

    psf_p gauss_psf ( psf::create ( gauss_props ) );

    size_t gauss_npix = gauss_psf->pixrows() * gauss_psf->pixcols();

    if ( gauss_npix != 4697 * 4110 ) {
      cerr << "gauss PSF image size (" << gauss_psf->pixrows() << " x " << gauss_psf->pixcols() << ") does not agree with parameters (4697 x 4110)" << endl;
      MPI_Abort( MPI_COMM_WORLD, 1 );
    }

    size_t gauss_nspec = gauss_psf->nspec();

    if ( gauss_nspec != 500 ) {
      cerr << "gauss PSF nspec (" << gauss_nspec << ") does not agree with file (500)" << endl;
      MPI_Abort( MPI_COMM_WORLD, 1 );
    }

    size_t gauss_nlambda = gauss_psf->nlambda();

    if ( gauss_nlambda != 4697 ) {
      cerr << "gauss PSF nlambda (" << gauss_nlambda << ") does not agree with file (4697)" << endl;
      MPI_Abort( MPI_COMM_WORLD, 1 );
    }

    vector < double > gauss_lambda = gauss_psf->lambda();

    if ( ( fabs( gauss_lambda[0] - 7460.0) > 1.0e-6 ) || ( fabs( gauss_lambda[gauss_nlambda-1] - 9808.0) > 1.0e-6 ) ) {
      cerr << "gauss PSF lambda (" << gauss_lambda[0] << " -- " << gauss_lambda[gauss_nlambda-1] << ") does not agree with file (7460.0 -- 9808.0)" << endl;
      MPI_Abort( MPI_COMM_WORLD, 1 );
    }

    // make fake spectra

    size_t fibers = 25;

    size_t gauss_global = gauss_nspec * gauss_nlambda;

    matrix_dist gauss_spec ( gauss_global, 1 );

    size_t delta = 10;

    for ( size_t i = 0; i < gauss_global; ++i ) {
      double val = 0.0;
      size_t spec = (size_t)(i / gauss_nlambda);
      size_t l = i - ( spec * gauss_nlambda );
      if ( l % delta == 0 ) {
        val = 100.0;
      }
      gauss_spec.Set ( i, 0, val );
    }

    // get design matrix from PSF and project

    matrix_sparse gauss_design;

    matrix_local gauss_image ( gauss_npix, 1 );
    local_matrix_zero ( gauss_image );

    gauss_psf->projection ( 0, fibers-1, 0, gauss_nlambda - 1, gauss_design );

    matrix_dist gauss_spec_block ( fibers * gauss_nlambda, 1 );

    sub_spec ( gauss_spec, gauss_nspec, 0, fibers, 0, gauss_nlambda, gauss_spec_block );

    spec_project ( gauss_design, gauss_spec_block, gauss_image );

    fitsfile * fp;

    if ( myp == 0 ) {
      string outimg = datadir + "/sandbox_gauss_project.fits.out";
      fits::create ( fp, outimg );
      fits::img_append ( fp, gauss_psf->pixrows(), gauss_psf->pixcols() );
      fits::img_write ( fp, gauss_image );
      fits::close ( fp );
    }

    if ( myp == 0 ) {
      cerr << "  (PASSED)" << endl;
    }
  } else {
    if ( myp == 0 ) {
      cerr << "Skipping sandbox PSF (file " << filepath << " not found)" << endl;
    }
  }


  if ( myp == 0 ) {
    cerr << "Testing sandbox fake PSF generation..." << endl;
  }

  // Create PSF

  boost::property_tree::ptree psf_props;
  psf_props.put ( "format", "sandbox" );
  psf_props.put ( "FAKE", "TRUE" );

  psf_p testpsf ( psf::create ( psf_props ) );

  size_t npix = testpsf->pixrows() * testpsf->pixcols();
  size_t nspec = testpsf->nspec();
  size_t nlambda = testpsf->nlambda();

  vector < double > psf_lambda = testpsf->lambda();

  size_t nbins = nspec * nlambda;

  boost::property_tree::ptree spec_props;
  spec_props.clear();
  spec_props.put ( "format", "sandbox" );
  spec_props.put ( "nspec", nspec );
  spec_props.put ( "nlambda", nlambda );
  spec_props.put ( "first_lambda", psf_lambda[0] );
  spec_props.put ( "last_lambda", psf_lambda[ nlambda - 1 ] );
  spec_props.put ( "back", 10.0 );
  spec_props.put ( "atm", 500.0 );
  spec_props.put ( "obj", 80.0 );
  spec_props.put ( "atmspace", 12 );
  spec_props.put ( "skymod", 25 );
  spec_p testspec ( spec::create ( spec_props ) );

  matrix_dist truth ( nbins, 1 );
  vector < double > lambda;
  vector < bool > sky;

  testspec->read( truth, lambda, sky );

  truth.Write("sndbx_truth");

  boost::property_tree::ptree img_props;
  img_props.put ( "format", "sandbox_sim" );
  img_props.put_child ( "psf", psf_props );
  img_props.put_child ( "spec", spec_props );

  image_p testimg ( image::create ( img_props ) );

  matrix_local measured ( npix, 1 );  
  local_matrix_zero ( measured );

  matrix_local invnoise ( npix, 1 );
  local_matrix_zero ( invnoise );

  testimg->read ( measured );
  testimg->read_noise ( invnoise );


  // generate design matrix and spectral subset for 
  // one slice of wavelength points.

  size_t nspec_slice = 25;
  size_t nlambda_slice = 100;
  size_t nbins_slice = nspec_slice * nlambda_slice;

  matrix_sparse design;

  double tstart;
  double tstop;

  tstart = MPI_Wtime();
  testpsf->projection ( 0, nspec_slice - 1, 0, nlambda_slice - 1, design );
  tstop = MPI_Wtime();
  if ( myp == 0 ) {
    cerr << "  Time for PSF creation = " << tstop-tstart << " seconds" << endl;
  }

  if ( myp == 0 ) {
    cerr << "Testing sandbox inverse covariance calculation..." << endl;
  }

  elem::Grid grid ( elem::mpi::COMM_WORLD );
  
  matrix_dist inv ( nbins_slice, nbins_slice, grid );

  tstart = MPI_Wtime();
  inverse_covariance ( design, invnoise, inv );
  tstop = MPI_Wtime();
  if ( myp == 0 ) {
    cerr << "  Time for inverse covariance building = " << tstop-tstart << " seconds" << endl;
  }

  if ( myp == 0 ) {
    cerr << "Testing sandbox inverse covariance eigendecomposition..." << endl;
  }

  matrix_dist W ( nbins_slice, nbins_slice, grid );
  matrix_dist D ( nbins_slice, 1, grid );

  tstart = MPI_Wtime();
  eigen_decompose ( inv, D, W );
  tstop = MPI_Wtime();
  if ( myp == 0 ) {
    cerr << "  Time for eigendecomposition = " << tstop-tstart << " seconds" << endl;
  }

  //W.Print("eigenvectors");

  ostringstream os;

  os.str("");
  os << "sndbx_D_" << np;
  D.Write( os.str() );

  if ( myp == 0 ) {
    cerr << "Testing sandbox resolution matrix..." << endl;
  }

  matrix_dist S ( nbins_slice, 1, grid );

  tstart = MPI_Wtime();
  norm ( D, W, S );
  tstop = MPI_Wtime();
  if ( myp == 0 ) {
    cerr << "  Time to compute S norm = " << tstop-tstart << " seconds" << endl;
  }

  os.str("");
  os << "sndbx_S_" << np;
  S.Write( os.str() );

  if ( myp == 0 ) {
    cerr << "Testing sandbox extraction..." << endl;
  }

  matrix_dist z ( nbins_slice, 1 );
  matrix_dist Rf ( nbins_slice, 1 );
  matrix_dist err ( nbins_slice, 1 );

  tstart = MPI_Wtime();
  noise_weighted_spec ( design, invnoise, measured, z );
  tstop = MPI_Wtime();
  if ( myp == 0 ) {
    cerr << "  Time to compute RHS = " << tstop-tstart << " seconds" << endl;
  }

  os.str("");
  os << "sndbx_z_" << np;
  z.Write( os.str() );

  tstart = MPI_Wtime();
  extract ( D, W, S, z, Rf, err );
  tstop = MPI_Wtime();
  if ( myp == 0 ) {
    cerr << "  Time for extraction = " << tstop-tstart << " seconds" << endl;
  }

  os.str("");
  os << "sndbx_Rf_" << np;
  Rf.Write( os.str() );

  matrix_dist Rtruth ( nbins_slice, 1 );
  matrix_dist R ( nbins_slice, nbins_slice, grid );

  tstart = MPI_Wtime();
  resolution ( D, W, S, R );
  tstop = MPI_Wtime();
  if ( myp == 0 ) {
    cerr << "  Time for R matrix construction = " << tstop-tstart << " seconds" << endl;
  }

  // (select truncated truth spectrum)

  matrix_dist truth_slice ( nbins_slice, 1 );
  dist_matrix_zero ( truth_slice );

  sub_spec ( truth, nbins, 0, nspec_slice, 0, nlambda_slice, truth_slice );

  elem::Gemv ( elem::NORMAL, 1.0, R, truth_slice, 0.0, Rtruth ); 

  os.str("");
  os << "sndbx_Rtruth_" << np;
  Rtruth.Write( os.str() );

  fstream fout;
  fout.precision(16);

  matrix_local full;

  fitsfile * fp;

  if ( myp == 0 ) {
    os.str("");
    os << "sndbx_results_" << np;
    fout.open ( os.str().c_str(), ios::out );

    full.ResizeTo ( nbins_slice, 1 );
    string outspec = datadir + "/sandbox_results.fits.out";
    fits::create ( fp, outspec );
  }

  for ( size_t i = 0; i < nbins_slice; ++i ) {
    double dtemp = truth.Get(i,0);
    if ( myp == 0 ) {
      full.Set(i,0, dtemp);
    }
  }
  if ( myp == 0 ) {
    fits::img_append ( fp, nspec, nlambda_slice );
    fits::img_write ( fp, full );
  }

  for ( size_t i = 0; i < nbins_slice; ++i ) {
    double dtemp = Rtruth.Get(i,0);
    if ( myp == 0 ) {
      full.Set(i,0, dtemp);
    }
  }
  if ( myp == 0 ) {
    fits::img_append ( fp, nspec, nlambda_slice );
    fits::img_write ( fp, full );
  }

  for ( size_t i = 0; i < nbins_slice; ++i ) {
    double dtemp = Rf.Get(i,0);
    if ( myp == 0 ) {
      full.Set(i,0, dtemp);
    }
  }
  if ( myp == 0 ) {
    fits::img_append ( fp, nspec, nlambda_slice );
    fits::img_write ( fp, full );
  }

  for ( size_t i = 0; i < nbins_slice; ++i ) {
    double dtemp = inv.Get(i,i);
    if ( myp == 0 ) {
      full.Set(i,0, sqrt(1.0/dtemp));
    }
  }
  if ( myp == 0 ) {
    fits::img_append ( fp, nspec, nlambda_slice );
    fits::img_write ( fp, full );
  }

  for ( size_t i = 0; i < nbins_slice; ++i ) {
    double out_rf = Rf.Get(i,0);
    double out_rt = Rtruth.Get(i,0);
    double out_tr = truth.Get(i,0);
    double errval = err.Get ( i, 0 );
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



