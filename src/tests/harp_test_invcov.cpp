
#include <iostream>

#include <harp_test.hpp>

#include <boost/random.hpp>

using namespace std;
using namespace harp;

#define DATASIZE 100
#define SIGSIZE 10
#define MAX 1000.0
#define TOL 1.0e-6

void harp::test_invcov ( string const & datadir ) {

  int np;
  int myp;

  MPI_Comm_size ( MPI_COMM_WORLD, &np );
  MPI_Comm_rank ( MPI_COMM_WORLD, &myp );

  if ( myp == 0 ) {
    cerr << "Testing inverse covariance construction..." << endl;
  }

  cerr.precision(16);
  
  // construct random sparse matrix

  matrix_sparse AT ( SIGSIZE, DATASIZE, elem::mpi::COMM_WORLD );

  matrix_local compAT ( SIGSIZE, DATASIZE );
  for ( size_t i = 0; i < SIGSIZE; ++i ) {
    for ( size_t j = 0; j < DATASIZE; ++j ) {
      compAT.Set ( i, j, 0.0 );
    }
  }

  size_t local_firstrow = AT.FirstLocalRow();
  size_t local_rows = AT.LocalHeight();

  typedef boost::ecuyer1988 base_generator_type;
  typedef boost::uniform_01<> distribution_type;
  typedef boost::variate_generator < base_generator_type&, distribution_type > gen_type;

  base_generator_type generator(42u);
  gen_type gen ( generator, distribution_type() );

  AT.StartAssembly();

  size_t nnz = 10;

  AT.Reserve ( nnz * local_rows );

  size_t col;

  for ( size_t i = 0; i < SIGSIZE; ++i ) {

    for ( size_t j = 0; j < nnz; ++j ) {

      col = (size_t)( 2 * i + j );
      //col = (size_t)( (double)DATASIZE * gen() );

      compAT.Set ( i, col, 1.0 );

      if ( ( i >= local_firstrow ) && ( i < local_firstrow + local_rows ) ) {
        //cerr << "update AT elem " << i << ", " << col << endl;
        AT.Update ( i, col, 1.0 );
      } else {
        //cerr << "ignore AT elem " << i << ", " << col << endl;
      }
    }

  }

  AT.StopAssembly();

  // fake noise covariance

  matrix_local invpix ( DATASIZE, 1 );
  local_matrix_zero ( invpix );

  for ( size_t i = 0; i < DATASIZE; ++i ) {
    invpix.Set ( i, 0, 0.1 );
  }

  // construct test output

  matrix_local compinv ( SIGSIZE, SIGSIZE );
  local_matrix_zero ( compinv );

  elem::Gemm ( elem::NORMAL, elem::TRANSPOSE, 1.0, compAT, compAT, 0.0, compinv );

  //compinv.Print ( "Serial inverse covariance" );

  // parallel implementation

  elem::Grid grid ( elem::mpi::COMM_WORLD );

  matrix_dist inv ( SIGSIZE, SIGSIZE, grid );

  inverse_covariance ( AT, invpix, inv );

  //inv.Print ( "MPI inverse covariance" );

  // compare results in lower triangle

  matrix_local local_inv ( SIGSIZE, SIGSIZE );
  local_matrix_zero ( local_inv );

  elem::AxpyInterface < double > globloc;
  globloc.Attach( elem::GLOBAL_TO_LOCAL, inv );

  globloc.Axpy ( 1.0, local_inv, 0, 0 );

  globloc.Detach();

  for ( size_t i = 0; i < SIGSIZE; ++i ) {
    for ( size_t j = i; j < SIGSIZE; ++j ) {
      double locval = local_inv.Get ( j, i );
      double compval = compinv.Get ( j, i );
      if ( fabs ( locval - compval ) / compval > TOL ) {
        cerr << "proc " << myp << " FAIL on element (" << j << ", " << i << "); " << locval << " != " << compval << endl;
        exit(1);
      }

    }
  }

  if ( myp == 0 ) {
    cerr << "  (PASSED)" << endl;
  }

  return;
}
