#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <limits>

#include <harp_test.hpp>

#include <boost/random.hpp>

using namespace std;
using namespace harp;


void harp::test_linalg ( string const & datadir ) {

  cout << "Testing ublas LAPACK bindings..." << endl;

  // matrix-vector multiply

  size_t dim = 100;

  matrix_double mat ( dim, dim );
  vector_double input ( dim );
  vector_double output ( dim );
  vector_double check ( dim );

  typedef boost::ecuyer1988 base_generator_type;
  base_generator_type generator(42u);

  boost::uniform_real < double > dist ( -100.0, 100.0 );
  boost::variate_generator < base_generator_type&, boost::uniform_real < double > > uni ( generator, dist );

  for ( size_t i = 0; i < dim; ++i ) {
    input[i] = 10.0;
    for ( size_t j = 0; j < dim; ++j ) {
      mat( j, i ) = uni();
    }
  }

  check.clear();

  for ( size_t i = 0; i < dim; ++i ) {
    for ( size_t j = 0; j < dim; ++j ) {
      check[i] += mat( i, j ) * input[j];
    }
  }

  boost::numeric::bindings::blas::gemv ( 1.0, mat, input, 0.0, output );

  for ( size_t i = 0; i < dim; ++i ) {
    if ( fabs ( ( output[i] - check[i] ) / check[i] ) > std::numeric_limits < double > :: epsilon() ) {
      cerr << "FAIL:  blas::gemv output element " << i << " is wrong (" << output[i]
 << " != " << check[i] << ")" << endl;
      exit(1);
    }
  }
  
  cout << "  (PASSED)" << endl;

  cout << "Testing eigen-decomposition..." << endl;

  // construct random matrices

  double rngmax = 1000.0;

  matrix_double a1 ( dim, dim );
  matrix_double a2 ( dim, dim );

  typedef boost::ecuyer1988 base_generator_type;
  typedef boost::uniform_01<> distribution_type;
  typedef boost::variate_generator < base_generator_type&, distribution_type > gen_type;

  gen_type gen ( generator, distribution_type() );

  for ( size_t i = 0; i < dim; ++i ) {
    for ( size_t j = 0; j < dim; ++j ) {
      a1 ( i, j ) = rngmax * gen();
      a2 ( i, j ) = rngmax * gen();
    }
  }

  // construct symmetric test matrix

  matrix_double sym ( dim, dim );

  boost::numeric::bindings::blas::gemm ( 1.0, boost::numeric::bindings::trans ( a1 ), a1, 0.0, sym );
  boost::numeric::bindings::blas::gemm ( 1.0, boost::numeric::bindings::trans ( a2 ), a2, 1.0, sym );

  // get eigenvectors and eigenvalues

  vector_double w;
  matrix_double Z;

  eigen_decompose ( sym, w, Z );

  matrix_double symprod ( dim, dim );
  matrix_double eprod ( dim, dim );
  matrix_double wdiag ( dim, dim );

  wdiag.clear();

  for ( size_t i = 0; i < dim; ++i ) {
    wdiag ( i, i ) = w[i];
  }

  boost::numeric::bindings::blas::gemm ( 1.0, sym, Z, 0.0, symprod );
  boost::numeric::bindings::blas::gemm ( 1.0, Z, wdiag, 0.0, eprod );

  double relerr;
  double eval, sval;

  for ( size_t i = 0; i < dim; ++i ) {
    for ( size_t j = 0; j < dim; ++j ) {
      eval = eprod ( j, i );
      sval = symprod ( j, i );
      relerr = fabs ( ( eval - sval ) / sval );
      if ( relerr > std::numeric_limits < float > :: epsilon() ) {
        cerr << "FAIL on matrix element (" << j << ", " << i << ") Av = " << sval << ", ev = " << eval << " rel err = " << relerr << endl;
        exit(1);
      }
    }
  }

  cout << "  (PASSED)" << endl;

  cout << "Testing re-composition..." << endl;

  matrix_double outcomp;

  eigen_compose ( EIG_NONE, w, Z, outcomp );

  double inval;
  double outval;

  for ( size_t i = 0; i < dim; ++i ) {
    for ( size_t j = 0; j < dim; ++j ) {
      inval = sym ( j, i );
      outval = outcomp ( j, i );
      relerr = fabs ( ( outval - inval ) / inval );
      if ( relerr > std::numeric_limits < float > :: epsilon() ) {
        cerr << "FAIL on matrix element (" << j << ", " << i << ") input = " << inval << ", output = " << outval << " rel err = " << relerr << endl;
        exit(1);
      }
    }
  }

  matrix_double mat_rt;
  eigen_compose ( EIG_SQRT, w, Z, mat_rt );

  matrix_double mat_invrt;
  eigen_compose ( EIG_INVSQRT, w, Z, mat_invrt );
  
  vector_double w_inv;
  matrix_double Z_inv;

  eigen_decompose ( mat_invrt, w_inv, Z_inv );

  for ( size_t i = 0; i < dim; ++i ) {
    w_inv[i] *= w_inv[i];
  }

  matrix_double comp_rt;

  eigen_compose ( EIG_INVSQRT, w_inv, Z_inv, comp_rt );

  for ( size_t i = 0; i < dim; ++i ) {
    for ( size_t j = 0; j < dim; ++j ) {
      inval = mat_rt ( j, i );
      outval = comp_rt ( j, i );
      relerr = fabs ( ( outval - inval ) / inval );
      if ( relerr > std::numeric_limits < float > :: epsilon() ) {
        cerr << "FAIL on doubly inverted sqrt matrix element (" << j << ", " << i << ") original = " << inval << ", output = " << outval << " rel err = " << relerr << endl;
        exit(1);
      }
    }
  }

  cout << "  (PASSED)" << endl;

  return;
}



