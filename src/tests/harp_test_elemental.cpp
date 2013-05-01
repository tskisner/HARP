
#include <iostream>

#include <harp_test.hpp>

#include <boost/random.hpp>

using namespace std;
using namespace harp;

#define SIZE 15
#define SPECSIZE 100
#define NSPEC 10
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

  elem::Grid grid ( MPI_COMM_WORLD );

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

  if ( myp == 0 ) {
    cerr << "Testing gang-parallel (re)distribution..." << endl;
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

  MPI_Comm gcomm;
  int ret = MPI_Comm_split ( MPI_COMM_WORLD, gang, grank, &gcomm );
  mpi_check ( MPI_COMM_WORLD, ret );

  // create gang process grids

  elem::Grid gang_grid ( gcomm );

  matrix_dist redist_comp ( outcomp );
  matrix_dist gout_comp ( SIZE, SIZE, gang_grid );

  gang_distribute ( outcomp, gout_comp );

  gang_accum ( gout_comp, redist_comp );

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
    cerr << "  (PASSED)" << endl;
  }

  if ( myp == 0 ) {
    cerr << "Testing gang-parallel slice and accum..." << endl;
  }

  matrix_dist fullspec ( NSPEC * SPECSIZE, 1, grid );
  matrix_dist comp_fullspec ( NSPEC * SPECSIZE, 1, grid );

  for ( size_t i = 0; i < NSPEC; ++i ) {
    for ( size_t j = 0; j < SPECSIZE; ++j ) {
      fullspec.Set ( i * SPECSIZE + j, 0, 100.0 + (double)j );
    }
  }

  matrix_dist gang_fullspec ( NSPEC * SPECSIZE, 1, gang_grid );

  gang_distribute ( fullspec, gang_fullspec );

  size_t nspec_chunk = (size_t)( NSPEC / 2 );
  size_t first_spec = (size_t)( NSPEC / 4 );
  size_t nlambda_chunk = (size_t)( SPECSIZE / 4 );
  size_t first_lambda = (size_t)( SPECSIZE / 8 );

  if ( myp == 0 ) {
    cout << "spec range = " << first_spec << " - " << (first_spec + nspec_chunk - 1) << endl;
    cout << "lambda range = " << first_lambda << " - " << (first_lambda + nlambda_chunk - 1) << endl;
  }

  matrix_dist gang_subspec ( nspec_chunk * nlambda_chunk, 1, gang_grid );

  sub_spec ( gang_fullspec, NSPEC, first_spec, nspec_chunk, first_lambda, nlambda_chunk, gang_subspec );

  MPI_Barrier ( MPI_COMM_WORLD );

  accum_spec ( gang_fullspec, NSPEC, first_spec, nspec_chunk, first_lambda, nlambda_chunk, gang_subspec );

  gang_accum ( gang_fullspec, comp_fullspec );

  MPI_Barrier ( MPI_COMM_WORLD );

  for ( size_t i = 0; i < NSPEC; ++i ) {
    for ( size_t j = 0; j < SPECSIZE; ++j ) {
      inval = fullspec.Get ( i * SPECSIZE + j, 0 ) * (double)ngang;
      if ( ( ( i >= first_spec ) && ( i < first_spec + nspec_chunk ) ) && ( ( j >= first_lambda ) && ( j < first_lambda + nlambda_chunk ) ) ) {
        inval *= 2.0;
      }
      outval = comp_fullspec.Get ( i * SPECSIZE + j, 0 );
      relerr = fabs ( outval - inval ) / inval;
      if ( relerr > TOL ) {
        cerr << "FAIL on spectrum " << i << ", wavelength " << j << ": " << outval << " != " << inval << endl;
        exit(1);
      }

    }
  }


  if ( myp == 0 ) {
    cerr << "  (PASSED)" << endl;
  }

  return;
}
