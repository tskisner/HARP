#include <iostream>
#include <fstream>
#include <sstream>

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

  double tstart;
  double tstop;

  tstart = MPI_Wtime();
  testpsf->projection ( 0, specsize - 1, design );
  tstop = MPI_Wtime();
  if ( myp == 0 ) {
    cerr << "  Time for PSF creation = " << tstop-tstart << " seconds" << endl;
  }

  matrix_dist truth ( nbins, 1 );

  boost::shared_ptr < psf_sandbox > sndpsf = testpsf->shared_ref < psf_sandbox > ();

  sndpsf->fake_spec ( truth );

  truth.Write("sndbx_truth");

  matrix_local signal ( npix, 1 );

  tstart = MPI_Wtime();
  spec_project ( design, truth, signal );
  tstop = MPI_Wtime();
  if ( myp == 0 ) {
    cerr << "  Time for spec to image projection = " << tstop-tstart << " seconds" << endl;
  }

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
    string outimg = datadir + "/sandbox_measured.fits.out";
    fitsfile * fp;
    fits::create ( fp, outimg );
    fits::img_append ( fp, testpsf->pixrows(), testpsf->pixcols() );
    fits::img_write ( fp, measured );
    fits::close ( fp );

    outimg = datadir + "/sandbox_signal.fits.out";
    fits::create ( fp, outimg );
    fits::img_append ( fp, testpsf->pixrows(), testpsf->pixcols() );
    fits::img_write ( fp, signal );
    fits::close ( fp );

    outimg = datadir + "/sandbox_noise.fits.out";
    fits::create ( fp, outimg );
    fits::img_append ( fp, testpsf->pixrows(), testpsf->pixcols() );
    fits::img_write ( fp, noise );
    fits::close ( fp );
  }

  if ( myp == 0 ) {
    cerr << "Testing sandbox inverse covariance calculation..." << endl;
  }

  elem::Grid grid ( elem::mpi::COMM_WORLD );
  
  matrix_dist inv ( nbins, nbins, grid );

  tstart = MPI_Wtime();
  inverse_covariance ( design, invnoise, inv );
  tstop = MPI_Wtime();
  if ( myp == 0 ) {
    cerr << "  Time for inverse covariance building = " << tstop-tstart << " seconds" << endl;
  }

  if ( myp == 0 ) {
    cerr << "Testing sandbox inverse covariance eigendecomposition..." << endl;
  }

  matrix_dist W ( nbins, nbins, grid );
  matrix_dist D ( nbins, 1, grid );

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

  matrix_dist S ( nbins, 1, grid );

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

  matrix_dist z ( nbins, 1 );
  matrix_dist Rf ( nbins, 1 );

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
  extract ( D, W, S, z, Rf );
  tstop = MPI_Wtime();
  if ( myp == 0 ) {
    cerr << "  Time for extraction = " << tstop-tstart << " seconds" << endl;
  }

  os.str("");
  os << "sndbx_Rf_" << np;
  Rf.Write( os.str() );

  matrix_dist Rtruth ( nbins, 1 );
  matrix_dist R ( nbins, nbins, grid );

  tstart = MPI_Wtime();
  resolution ( D, W, S, R );
  tstop = MPI_Wtime();
  if ( myp == 0 ) {
    cerr << "  Time for R matrix construction = " << tstop-tstart << " seconds" << endl;
  }

  elem::Gemv ( elem::NORMAL, 1.0, R, truth, 0.0, Rtruth ); 

  os.str("");
  os << "sndbx_Rtruth_" << np;
  Rtruth.Write( os.str() );

  fstream fout;
  fout.precision(16);

  if ( myp == 0 ) {
    os.str("");
    os << "sndbx_results_" << np;
    fout.open ( os.str().c_str(), ios::out );
  }

  for ( size_t i = 0; i < nbins; ++i ) {
    double out_rf = Rf.Get(i,0);
    double out_rt = Rtruth.Get(i,0);
    double out_tr = truth.Get(i,0);
    double err = inv.Get ( i, i );
    if ( myp == 0 ) {
      err = sqrt( 1.0 / err );
      fout << i << " " << out_tr << " " << out_rt << " " << out_rf << " " << err << endl;
    }
  }

  if ( myp == 0 ) {
    fout.close();
  }


  if ( myp == 0 ) {
    cerr << "  (PASSED)" << endl;
  }
     
  return;
}



