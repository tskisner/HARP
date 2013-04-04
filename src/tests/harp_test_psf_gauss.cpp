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
      cerr << "Testing elliptical gaussian PSF..." << endl;
    }

    boost::property_tree::ptree gauss_props;
    gauss_props.put ( "format", "gauss" );
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

    if ( myp == 0 ) {
      cerr << "  (PASSED)" << endl;
    }
  } else {
    if ( myp == 0 ) {
      cerr << "Skipping gauss PSF (file " << filepath << " not found)" << endl;
    }
  }

  return;
}


