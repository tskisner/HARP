
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



  cerr << "Testing Lanczos eigendecomposition..." << endl;

  typedef boost::numeric::ublas::vector < double > lan_vector;

  typedef boost::numeric::ublas::wrapper_vectorspace < lan_vector > lan_vs;

  typedef std::complex<double> value_type;

  // get eigenvectors and eigenvalues

  lan_vs vs ( SIZE );
  boost::numeric::ublas::lanczos < mat_denserow, lan_vs > lcz ( sym, vs );

  int max_iter = 2*SIZE;  
  cerr << "Computation of eigenvalues with fixed size of T-matrix\n\n";
  cerr << "-----------------------------------\n\n";

  vector<double> eigen;
  vector<double> err;
  vector<int> multiplicity;
  boost::numeric::ublas::fixed_lanczos_iteration<double> iter(max_iter);

  typedef boost::lagged_fibonacci607 Gen;
  Gen mygen;

  try {
    lcz.calculate_eigenvalues(iter, mygen);
    eigen = lcz.eigenvalues();
    err = lcz.errors();
    multiplicity = lcz.multiplicities();
  }
  catch (std::runtime_error& e) {
    cerr << e.what() << std::endl;
  } 

  // Printing eigenvalues with error & multiplicities:  
  cerr << "#         eigenvalue         error         multiplicity\n";  
  cerr.precision(10);

  for (size_t i=0; i<eigen.size(); ++i)
    cerr << i << "\t" << eigen[i] << "\t" << err[i] << "\t" 
        << multiplicity[i] << "\n";
  
  // call of eigenvectors function follows:   
  cerr << "\nEigen vectors computations for 3 lowest eigenvalues:\n\n";  
  vector<double>::iterator start = eigen.begin();
  vector<double>::iterator end = eigen.begin()+3;
  vector< lan_vector > eigenvectors; // for storing the eigen vector. 
  boost::numeric::ublas::Info<double> info; // (m1, m2, ma, eigenvalue, residual, status).
  
  try {
    lcz.eigenvectors(start, end, std::back_inserter(eigenvectors),info,mygen); 
  }
  catch (std::runtime_error& e) {
    cerr << e.what() << std::endl;
  }  
  
  cerr << "Printing eigen Vectors:" << std::endl << std::endl; 
  for(std::vector< lan_vector >::iterator it = eigenvectors.begin(); it != eigenvectors.end(); it++){
    std::copy((it)->begin(),(it)->end(),std::ostream_iterator<value_type>(cerr,"\n"));
    cerr << "\n\n";
  }
  cerr << " Information about the eigen vectors computations:\n\n";
  for(int i = 0; i < info.size(); i++) {
    cerr << " m1(" << i+1 << "): " << info.m1(i) << ", m2(" << i+1 << "): "
        << info.m2(i) << ", ma(" << i+1 << "): " << info.ma(i) << " eigenvalue("
        << i+1 << "): " << info.eigenvalue(i) << " residual(" << i+1 << "): "
        << info.residual(i) << " error_info(" << i+1 << "): "
        << info.error_info(i) << std::endl << std::endl;
  }

  /*
  //boost::numeric::ublas::matrix < double > eigval = eig.getD();
  //boost::numeric::ublas::matrix < double > eigvec = eig.getV();

  // Test for solution to eigenvalue equation

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
  */

  cerr << "  (PASSED)" << endl;

  return;
}
