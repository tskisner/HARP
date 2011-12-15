
#include <iostream>

#include <harp_test.hpp>

#include <boost/random.hpp>

using namespace std;
using namespace harp;

#define SIZE 100
#define MAX 1000.0
#define TOL 1.0e-6

void harp::test_eigen ( string const & datadir ) {

  cerr << "Testing dense eigendecomposition..." << endl;

  cerr.precision(16);
  
  // construct random matrices

  mat_denserow a1 ( SIZE, SIZE );
  mat_denserow a2 ( SIZE, SIZE );

  typedef boost::ecuyer1988 base_generator_type;
  typedef boost::uniform_01<> distribution_type;
  typedef boost::variate_generator < base_generator_type&, distribution_type > gen_type;

  base_generator_type generator(42u);
  gen_type gen ( generator, distribution_type() );

  mat_denserow :: iterator1 rowit;
  mat_denserow :: iterator2 colit;

  for ( size_t i = 0; i < SIZE; ++i ) {
    for ( size_t j = 0; j < SIZE; ++j ) {
      a1( i, j ) = MAX * gen();
      a2( i, j ) = MAX * gen();
    }
  }

  // construct symmetric test matrix

  mat_denserow sym ( SIZE, SIZE );

  boost::numeric::ublas::axpy_prod ( boost::numeric::ublas::trans ( a1 ), a1, sym, true );
  boost::numeric::ublas::axpy_prod ( boost::numeric::ublas::trans ( a2 ), a2, sym, false );

  // get eigenvectors and eigenvalues

  boost::numeric::ublas::EigenvalueDecomposition eig ( sym );

  boost::numeric::ublas::matrix < double > eigval = eig.getD();
  boost::numeric::ublas::matrix < double > eigvec = eig.getV();

  // Test for solution to eigenvalue equation

  mat_denserow symprod ( SIZE, SIZE );
  mat_denserow eprod ( SIZE, SIZE );

  boost::numeric::ublas::axpy_prod ( sym, eigvec, symprod, true );
  boost::numeric::ublas::axpy_prod ( eigvec, eigval, eprod, true );

  double relerr;

  for ( size_t i = 0; i < SIZE; ++i ) {
    for ( size_t j = 0; j < SIZE; ++j ) {
      relerr = fabs ( eprod(i,j) - symprod(i,j) ) / symprod(i,j);
      if ( relerr > TOL ) {
         cerr << "FAIL on matrix element (" << i << ", " << j << ") Av = " << symprod(i,j) << ", ev = " << eprod(i,j) << " rel err = " << relerr << endl;
         exit(1);
      }
    }
  }
  
  cerr << "  (PASSED)" << endl;

  return;
}
