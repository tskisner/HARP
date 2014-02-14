
#include <iostream>

#include <harp_mpi_test.hpp>

#include <boost/random.hpp>

using namespace std;
using namespace harp;

#define SIZE 15
#define SPECSIZE 100
#define NSPEC 10
#define MAX 1000.0
#define TOL 1.0e-6


void harp::mpi_test_linalg ( string const & datadir ) {

  boost::mpi::communicator comm;

  int np = comm.size();
  int myp = comm.rank();

  if ( myp == 0 ) {
    cout << "Testing elemental eigendecomposition..." << endl;
  }

  cerr.precision(16);
  
  // construct random matrices

  elem::Grid grid ( comm );

  mpi_matrix a1 ( SIZE, SIZE, grid );
  mpi_matrix a2 ( SIZE, SIZE, grid );

  typedef boost::ecuyer1988 base_generator_type;
  typedef boost::uniform_01<> distribution_type;
  typedef boost::variate_generator < base_generator_type&, distribution_type > gen_type;

  base_generator_type generator(42u);
  gen_type gen ( generator, distribution_type() );

  for ( size_t i = 0; i < SIZE; ++i ) {
    for ( size_t j = 0; j < SIZE; ++j ) {
      a1.Set( i, j, MAX * gen() );
      a2.Set( i, j, MAX * gen() );
    }
  }

  // construct symmetric test matrix

  mpi_matrix sym ( SIZE, SIZE, grid );

  elem::Gemm ( elem::TRANSPOSE, elem::NORMAL, 1.0, a1, a1, 0.0, sym );
  elem::Gemm ( elem::TRANSPOSE, elem::NORMAL, 1.0, a2, a2, 1.0, sym );

  // get eigenvectors and eigenvalues

  mpi_matrix w;
  mpi_matrix Z;

  mpi_eigen_decompose ( sym, w, Z );

  mpi_matrix symprod ( SIZE, SIZE, grid );
  mpi_matrix eprod ( SIZE, SIZE, grid );
  mpi_matrix wdiag ( SIZE, SIZE, grid );

  mpi_matrix_zero ( wdiag );
  for ( size_t i = 0; i < SIZE; ++i ) {
    double wval = w.Get(i,0);
    wdiag.Set ( i, i, wval );
  }

  elem::Gemm ( elem::NORMAL, elem::NORMAL, 1.0, sym, Z, 0.0, symprod );
  elem::Gemm ( elem::NORMAL, elem::NORMAL, 1.0, Z, wdiag, 0.0, eprod );

  double relerr;
  double eval, sval;

  for ( size_t i = 0; i < SIZE; ++i ) {
    for ( size_t j = 0; j < SIZE; ++j ) {
      eval = eprod.Get ( j, i );
      sval = symprod.Get ( j, i );
      relerr = fabs ( eval - sval ) / sval;
      if ( relerr > TOL ) {
        cerr << "FAIL on matrix element (" << j << ", " << i << ") Av = " << sval << ", ev = " << eval << " rel err = " << relerr << endl;
        exit(1);
      }
    }
  }

  if ( myp == 0 ) {
    cout << "  (PASSED)" << endl;
  }

  if ( myp == 0 ) {
    cout << "Testing re-composition..." << endl;
  }

  mpi_matrix outcomp;

  mpi_eigen_compose ( EIG_NONE, w, Z, outcomp );

  double inval;
  double outval;

  for ( size_t i = 0; i < SIZE; ++i ) {
    for ( size_t j = 0; j < SIZE; ++j ) {
      inval = sym.Get ( j, i );
      outval = outcomp.Get ( j, i );
      relerr = fabs ( outval - inval ) / inval;
      if ( relerr > TOL ) {
        cerr << "FAIL on matrix element (" << j << ", " << i << ") input = " << inval << ", output = " << outval << " rel err = " << relerr << endl;
        exit(1);
      }
    }
  }

  mpi_matrix mat_rt;
  mpi_eigen_compose ( EIG_SQRT, w, Z, mat_rt );

  mpi_matrix mat_invrt;
  mpi_eigen_compose ( EIG_INVSQRT, w, Z, mat_invrt );
  
  mpi_matrix w_inv;
  mpi_matrix Z_inv;

  mpi_eigen_decompose ( mat_invrt, w_inv, Z_inv );

  for ( size_t i = 0; i < SIZE; ++i ) {
    double val = w_inv.Get(i,0);
    val *= val;
    w_inv.Set ( i, 0, val );
  }

  mpi_matrix comp_rt;

  mpi_eigen_compose ( EIG_INVSQRT, w_inv, Z_inv, comp_rt );

  for ( size_t i = 0; i < SIZE; ++i ) {
    for ( size_t j = 0; j < SIZE; ++j ) {
      inval = mat_rt.Get ( j, i );
      outval = comp_rt.Get ( j, i );
      relerr = fabs ( outval - inval ) / inval;
      if ( relerr > TOL ) {
        cerr << "FAIL on doubly inverted sqrt matrix element (" << j << ", " << i << ") original = " << inval << ", output = " << outval << " rel err = " << relerr << endl;
        exit(1);
      }
    }
  }

  if ( myp == 0 ) {
    cout << "  (PASSED)" << endl;
  }

  if ( myp == 0 ) {
    cout << "Testing gang-parallel (re)distribution..." << endl;
  }

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

  // create gang process grids

  elem::Grid gang_grid ( gcomm );

  mpi_matrix redist_comp ( outcomp );
  mpi_matrix gout_comp ( SIZE, SIZE, gang_grid );

  mpi_gang_distribute ( outcomp, gout_comp );

  mpi_gang_accum ( gout_comp, redist_comp );

  for ( size_t i = 0; i < SIZE; ++i ) {
    for ( size_t j = 0; j < SIZE; ++j ) {
      inval = (double)ngang * outcomp.Get ( j, i );
      outval = redist_comp.Get ( j, i );
      relerr = fabs ( outval - inval ) / inval;
      if ( relerr > TOL ) {
        cerr << "FAIL on matrix element (" << j << ", " << i << ") input = " << inval << ", output = " << outval << " rel err = " << relerr << endl;
        exit(1);
      }
    }
  }

  if ( myp == 0 ) {
    cout << "  (PASSED)" << endl;
  }


  return;
}
