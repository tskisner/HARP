#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <limits>

#include <harp_test.hpp>

#include <boost/random.hpp>

using namespace std;
using namespace harp;


void harp::test_lapack ( string const & datadir ) {

  cerr << "Testing ublas LAPACK bindings..." << endl;

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
  

  cerr << "  (PASSED)" << endl;

  return;
}



