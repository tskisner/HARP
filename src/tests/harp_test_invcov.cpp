
#include <iostream>

#include <harp_test.hpp>

#include <boost/random.hpp>

using namespace std;
using namespace harp;

#define DATASIZE 100
#define SIGSIZE 6
#define MAX 1000.0
#define TOL 1.0e-6

void harp::test_invcov ( string const & datadir ) {

  int np;
  int myp;

  MPI_Comm_size ( MPI_COMM_WORLD, &np );
  MPI_Comm_rank ( MPI_COMM_WORLD, &myp );

  elem::Grid grid ( elem::mpi::COMM_WORLD );

  if ( myp == 0 ) {
    cerr << "Testing inverse covariance construction..." << endl;
  }

  cerr.precision(16);
  
  // construct random sparse matrix

  matrix_sparse AT ( SIGSIZE, DATASIZE, elem::mpi::COMM_WORLD );

  matrix_local compAT ( SIGSIZE, DATASIZE );
  local_matrix_zero ( compAT );

  size_t local_firstrow = AT.FirstLocalRow();
  size_t local_rows = AT.LocalHeight();

  typedef boost::ecuyer1988 base_generator_type;
  base_generator_type generator(42u);

  double rms = 0.5;
  boost::normal_distribution < double > dist ( 0.0, rms );
  boost::variate_generator < base_generator_type&, boost::normal_distribution < double > > gauss ( generator, dist );


  AT.StartAssembly();

  size_t nnz = 10;

  AT.Reserve ( nnz * local_rows );

  size_t col;

  for ( size_t i = 0; i < SIGSIZE; ++i ) {

    for ( size_t j = 0; j < nnz; ++j ) {

      col = (size_t)( 2 * i + 8*j );

      compAT.Set ( i, col, 1.0 );

      if ( ( i >= local_firstrow ) && ( i < local_firstrow + local_rows ) ) {
        AT.Update ( i, col, 1.0 );
      }

    }

  }

  AT.StopAssembly();

  // truth

  matrix_dist truth ( SIGSIZE, 1, grid );

  matrix_local signal ( DATASIZE, 1 );

  for ( int i = 0; i < SIGSIZE; ++i ) {
    truth.Set ( i, 0, 42.0 );
  }

  spec_project ( AT, truth, signal );
  
  // fake noise covariance and measured data

  matrix_local noise ( DATASIZE, 1 );

  matrix_local measured ( DATASIZE, 1 );

  matrix_local invpix ( DATASIZE, 1 );

  for ( size_t i = 0; i < DATASIZE; ++i ) {
    invpix.Set ( i, 0, 1.0/(rms * rms) );
    noise.Set ( i, 0, gauss() );
    measured.Set ( i, 0, signal.Get(i,0) + noise.Get(i,0) );
  }

  // RHS

  matrix_dist z ( SIGSIZE, 1, grid );

  noise_weighted_spec ( AT, invpix, measured, z );

  
  // construct test output

  matrix_local compinv ( SIGSIZE, SIGSIZE );
  local_matrix_zero ( compinv );

  matrix_local compATN ( compAT );

  for ( size_t i = 0; i < SIGSIZE; ++i ) {

    for ( size_t j = 0; j < nnz; ++j ) {

      col = (size_t)( 2 * i + 8*j );

      compATN.Set ( i, col, 1.0/(rms * rms) );

    }

  }

  elem::Gemm ( elem::NORMAL, elem::TRANSPOSE, 1.0, compATN, compAT, 0.0, compinv );

  if ( myp == 0 ) {
    compinv.Print ( "Serial inverse covariance" );
  }

  // parallel implementation

  matrix_dist inv ( SIGSIZE, SIGSIZE, grid );

  inverse_covariance ( AT, invpix, inv );

  inv.Print ( "MPI inverse covariance" );

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
    cerr << "Testing extraction..." << endl;
  }

  matrix_dist W ( SIGSIZE, SIGSIZE, grid );
  matrix_dist D ( SIGSIZE, 1, grid );
  eigen_decompose ( inv, D, W );

  matrix_dist sq ( SIGSIZE, SIGSIZE, grid );
  eigen_compose ( EIG_SQRT, D, W, sq );

  sq.Print ( "sqrt of invcov" );

  matrix_dist Rdirect ( sq );
  matrix_dist R ( SIGSIZE, SIGSIZE, grid );
  matrix_dist S ( SIGSIZE, 1, grid );

  norm ( D, W, S );

  S.Print ( "column norm from eigen decomposition" );

  apply_norm ( S, Rdirect );

  Rdirect.Print ( "direct resolution matrix" );

  resolution ( D, W, S, R );

  R.Print ( "resolution matrix" );

  matrix_dist f ( SIGSIZE, 1, grid );

  extract ( D, W, S, z, f );

  for ( size_t i = 0; i < SIGSIZE; ++i ) {
    double rfdir = 0.0;
    double tr = truth.Get ( i, 0 );
    double rf = f.Get ( i, 0 );
    double err = inv.Get ( i, i );
    double ze = z.Get ( i, 0 );
    for ( size_t j = 0; j < SIGSIZE; ++j ) {
      rfdir += R.Get(i,j) * f.Get(j,0);
    }
    if ( myp == 0 ) {
      err = sqrt( 1.0 / err );
      cout << "  truth = " << tr << ", z = " << ze << ", convolved = " << rfdir << ", Rf = " << rf << ", err = " << err << endl;
    }
  }
  

  if ( myp == 0 ) {
    cerr << "  (PASSED)" << endl;
  }


  return;
}