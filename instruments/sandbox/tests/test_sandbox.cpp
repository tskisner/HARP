#include <iostream>
#include <fstream>

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

  MPI_Comm_size ( MPI_COMM_WORLD, &np );
  MPI_Comm_rank ( MPI_COMM_WORLD, &myp );
  
  if ( myp == 0 ) {
    cerr << "Testing sandbox fake PSF generation..." << endl;
  }

  boost::property_tree::ptree props;
  props.put ( "FAKE", "TRUE" );

  psf_p testpsf ( psf::create ( "sandbox", props ) );

  size_t npix = testpsf->pixrows() * testpsf->pixcols();
  size_t nspec = testpsf->nspec();
  size_t specsize = testpsf->specsize();
  size_t nbins = nspec * specsize;

  matrix_sparse design;

  cerr << "constructing design matrix" << endl;
  testpsf->projection ( 0, specsize - 1, design );

  cerr << "allocate truth" << endl;
  matrix_dist truth ( nbins, 1 );

  cerr << "pointer cast" << endl;
  boost::shared_ptr < psf_sandbox > sndpsf = testpsf->shared_ref < psf_sandbox > ();
  cerr << "getting fake spec" << endl;
  sndpsf->fake_spec ( truth );

  truth.Write("sndbx_truth");

  matrix_local signal ( npix, 1 );
  cerr << "projecting to pixels" << endl;
  spec_project ( design, truth, signal );

  matrix_local noise ( npix, 1 );  
  local_matrix_zero ( noise );

  matrix_local measured ( npix, 1 );  
  local_matrix_zero ( measured );

  matrix_local invnoise ( npix, 1 );
  local_matrix_zero ( invnoise );

  typedef boost::ecuyer1988 base_generator_type;
  base_generator_type generator(42u);

  double rms;

  for ( size_t i = 0; i < npix; ++i ) {
    rms = sqrt( 16.0 + signal.Get ( i, 0 ) );
    invnoise.Set ( i, 0, 1.0 / (rms*rms) );

    boost::normal_distribution < double > dist ( 0.0, rms );
    boost::variate_generator < base_generator_type&, boost::normal_distribution < double > > gauss ( generator, dist );

    noise.Set ( i, 0, gauss() );
    measured.Set ( i, 0, signal.Get(i,0) + noise.Get(i,0) );
  }

  if ( myp == 0 ) {
    cerr << "Testing sandbox inverse covariance calculation..." << endl;
  }

  elem::Grid grid ( elem::mpi::COMM_WORLD );
  
  matrix_dist inv ( nbins, nbins, grid );

  inverse_covariance ( design, invnoise, inv );

  if ( myp == 0 ) {
    cerr << "Testing sandbox inverse covariance eigendecomposition..." << endl;
  }

  cerr << "matrix dims = " << nbins << " x " << nbins << endl;

  matrix_dist W ( nbins, nbins, grid );
  matrix_dist D ( nbins, 1, grid );
  eigen_decompose ( inv, D, W );

  D.Write("sndbx_D");

  if ( myp == 0 ) {
    cerr << "Testing sandbox resolution matrix..." << endl;
  }

  matrix_dist R ( nbins, nbins, grid );
  matrix_dist S ( nbins, 1, grid );

  norm ( D, W, S );
  //resolution ( D, W, S, R );

  S.Write("sndbx_S");

  if ( myp == 0 ) {
    cerr << "Testing sandbox extraction..." << endl;
  }

  matrix_dist z ( nbins, 1 );
  matrix_dist Rf ( nbins, 1 );

  noise_weighted_spec ( design, invnoise, measured, z );

  z.Write("sndbx_z");

  extract ( D, W, S, z, Rf );

  Rf.Write("sndbx_Rf");

  if ( myp == 0 ) {
    cerr << "  (PASSED)" << endl;
  }
     
  return;
}



