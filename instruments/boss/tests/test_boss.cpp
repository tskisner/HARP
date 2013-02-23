#include <iostream>
#include <fstream>
#include <sstream>

#include <harp_test.hpp>

#include <boost/random.hpp>

extern "C" {
  #include <unistd.h>
  #include <sys/stat.h>
}

using namespace std;
using namespace harp;


void harp::test_boss ( string const & datadir ) {

  int np;
  int myp;
  int ret;

  MPI_Comm_size ( MPI_COMM_WORLD, &np );
  MPI_Comm_rank ( MPI_COMM_WORLD, &myp );

  string filepath = datadir + "/ted_spec_2012-02-12.fits";

  int statret;
  if ( myp == 0 ) {
    struct stat statbuf;
    statret = stat ( filepath.c_str(), &statbuf );
  }

  ret = MPI_Bcast ( (void*)(&statret), 1, MPI_INT, 0, MPI_COMM_WORLD );
  mpi_check ( MPI_COMM_WORLD, ret );
    
  if ( statret == 0 ) {
    if ( myp == 0 ) {
      cerr << "Testing boss_specter operations..." << endl;
    }

    boost::property_tree::ptree props;
    props.put ( "format", "boss_specter" );
    props.put ( "path", filepath );

    spec_p testspec ( spec::create ( props ) );

    size_t nspec = testspec->nspectrum();
    if ( nspec != 500 ) {
      cerr << "FAIL:  number of spectra (" << nspec << ") is not 500" << endl;
      exit(1);
    }

    size_t nlambda = testspec->nlambda();
    if ( nlambda != 4697 ) {
      cerr << "FAIL:  number of wavelength points (" << nlambda << ") is not 4697" << endl;
      exit(1);
    }

    size_t nglobal = nspec * nlambda;

    matrix_dist specdata ( nglobal, 1 );
    vector < double > lambda;
    vector < bool > sky;

    testspec->read ( specdata, lambda, sky );

    string outfile = datadir + "/boss_specter_data.out";

    specdata.Write( outfile );

    if ( myp == 0 ) {
      fstream fout;
      fout.precision(16);

      outfile = datadir + "/boss_specter_lambda.out";

      fout.open ( outfile.c_str(), ios::out );
      for ( size_t i = 0; i < nlambda; ++i ) {
        fout << lambda[i] << endl;
      }
      fout.close();

      outfile = datadir + "/boss_specter_sky.out";

      fout.open ( outfile.c_str(), ios::out );
      for ( size_t i = 0; i < nspec; ++i ) {
        fout << sky[i] << endl;
      }
      fout.close();

    }
    
    if ( myp == 0 ) {
      cerr << "  (PASSED)" << endl;
    }
  } else {
    if ( myp == 0 ) {
      cerr << "Skipping boss_specter tests (file " << filepath << " not found)" << endl;
    }
  }
     
  return;
}



