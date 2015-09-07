#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include <harp_mpi_test.hpp>


using namespace std;
using namespace harp;


void harp::mpi_test_specslice ( string const & datadir ) {

  // for this test, we pretend that every process in the global communicator is its own "gang"

  boost::mpi::communicator comm;
  int np = comm.size();
  int myp = comm.rank();

  boost::mpi::communicator rcomm = comm.split ( myp, 0 );


  if ( myp == 0 ) {
    cout << "Testing spec slicing..." << endl;
  }

  string specfile = datadir + "/spec_sim.fits.out";

  size_t nspec = 20;
  size_t nlambda = 100;
  size_t overlap_spec = 1;
  size_t overlap_lambda = 10;
  size_t chunk_spec = 2;
  size_t chunk_lambda = 10;

  size_t nchunk_spec = (size_t)( nspec / chunk_spec );
  size_t nchunk_lambda = (size_t)( nlambda / chunk_lambda );
  size_t nchunk = nchunk_spec * nchunk_lambda;

  mpi_spec_slice_p slice ( new mpi_spec_slice ( rcomm, comm, overlap_spec, overlap_lambda, nspec, nlambda, chunk_spec, chunk_lambda, overlap_spec, overlap_lambda ) );

  matrix_float coverage ( slice->full_region().n_spec, slice->full_region().n_lambda );
  coverage.clear();

  std::vector < spec_slice_region > myslice = slice->regions();

  for ( std::vector < spec_slice_region > :: const_iterator sit = myslice.begin(); sit != myslice.end(); ++sit ) {

    for ( size_t i = 0; i < sit->n_good_spec; ++i ) {

      for ( size_t j = 0; j < sit->n_good_lambda; ++j ) {

        size_t global_spec = sit->first_good_spec + i;
        size_t global_lambda = sit->first_good_lambda + j;

        coverage ( global_spec, global_lambda ) += 1.0;

      }

    }

  }

  matrix_float full_coverage ( slice->full_region().n_spec, slice->full_region().n_lambda );
  full_coverage.clear();

  boost::mpi::reduce ( rcomm, coverage, full_coverage, std::plus < matrix_float > (), 0 );

  // check that every spectral point is covered once

  if ( myp == 0 ) {

    for ( size_t i = 0; i < nspec; ++i ) {

      for ( size_t j = 0; j < nlambda; ++j ) {

        size_t sp = overlap_spec + i;
        size_t la = overlap_lambda + j;

        if ( ( full_coverage(sp,la) < 0.5 ) || ( full_coverage(sp,la) > 1.5 ) ) {
          ostringstream o;
          o << "FAIL:  spectrum " << sp << ", lambda " << la << " was assigned to " << full_coverage(sp,la) << " workers";
          cerr << o.str() << endl;
          exit(1);
        }

      }
    
    }

    cout << "  (PASSED)" << endl;

  }

  return;
}



