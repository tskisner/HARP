
#include <iostream>

#include <harp_test.hpp>

#include <boost/random.hpp>

using namespace std;
using namespace harp;

#define SIZE 10
#define MAX 1000.0
#define TOL 1.0e-6

void harp::test_elemental ( string const & datadir ) {

  int np;
  int myp;

  MPI_Comm_size ( MPI_COMM_WORLD, &np );
  MPI_Comm_rank ( MPI_COMM_WORLD, &myp );

  if ( myp == 0 ) {
    cerr << "Testing elemental eigendecomposition..." << endl;
  }

  cerr.precision(16);
  
  // construct random matrices

  elem::Grid grid;

  matrix_dist a1 ( SIZE, SIZE, grid );
  matrix_dist a2 ( SIZE, SIZE, grid );

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

  matrix_dist sym ( SIZE, SIZE, grid );

  elem::Gemm ( elem::TRANSPOSE, elem::NORMAL, 1.0, a1, a1, 0.0, sym );
  elem::Gemm ( elem::TRANSPOSE, elem::NORMAL, 1.0, a2, a2, 1.0, sym );

  // get eigenvectors and eigenvalues

  matrix_dist w;
  matrix_dist Z;

  eigen_decompose ( sym, w, Z );

  matrix_dist symprod ( SIZE, SIZE, grid );
  matrix_dist eprod ( SIZE, SIZE, grid );
  matrix_dist wdiag ( SIZE, SIZE, grid );

  dist_matrix_zero ( wdiag );
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
    cerr << "Testing re-composition..." << endl;
  }

  matrix_dist outcomp;

  eigen_compose ( EIG_NONE, w, Z, outcomp );

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


  if ( myp == 0 ) {
    cerr << "  (PASSED)" << endl;
  }


  return;
}