#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include <harp_test.hpp>


using namespace std;
using namespace harp;


void harp::test_specslice ( string const & datadir ) {

  cout << "Testing spec slicing..." << endl;

  string specfile = datadir + "/spec_sim.fits.out";

  size_t nspec = 20;
  size_t nlambda = 101;
  size_t overlap_spec = 1;
  size_t overlap_lambda = 10;
  size_t chunk_spec = 3;
  size_t chunk_lambda = 13;

  size_t nchunk_spec = (size_t)( nspec / chunk_spec );
  if ( ! ( nspec % chunk_spec == 0 ) ) {
    ++nchunk_spec;
  }
  size_t nchunk_lambda = (size_t)( nlambda / chunk_lambda );
  if ( ! ( nlambda % chunk_lambda == 0 ) ) {
    ++nchunk_lambda;
  }
  size_t nchunk = nchunk_spec * nchunk_lambda;

  size_t procs = 100;

  spec_slice_p slice_one ( new spec_slice ( 1, nspec, nlambda, chunk_spec, chunk_lambda, overlap_spec, overlap_lambda ) );

  spec_slice_p slice ( new spec_slice ( procs, nspec, nlambda, chunk_spec, chunk_lambda, overlap_spec, overlap_lambda ) );

  // immediately serialize and restore, so that any issues with that process will impact the code that follows

  string serialpath = datadir + "/specslice_serialize.xml.out";
  {
    ofstream ofs ( serialpath.c_str() );
    boost::archive::xml_oarchive oa ( ofs );
    oa << BOOST_SERIALIZATION_NVP(slice_one);
    oa << BOOST_SERIALIZATION_NVP(slice);
  }
  {
    ifstream ifs ( serialpath.c_str() );
    boost::archive::xml_iarchive ia ( ifs );
    ia >> BOOST_SERIALIZATION_NVP(slice_one);
    ia >> BOOST_SERIALIZATION_NVP(slice);
  }

  std::vector < spec_slice_region > checkone = slice_one->regions(0);

  if ( checkone.size() != nchunk ) {
    cerr << "FAIL:  spec slicing number of regions (" << checkone.size() << ") does not match correct value (" << nchunk << ")" << endl;
    exit(1);
  }

  matrix_float coverage ( nspec, nlambda );
  coverage.clear();

  for ( size_t p = 0; p < procs; ++p ) {
    std::vector < spec_slice_region > procslice = slice->regions ( p );

    for ( std::vector < spec_slice_region > :: const_iterator sit = procslice.begin(); sit != procslice.end(); ++sit ) {

      for ( size_t i = 0; i < sit->n_good_spec; ++i ) {

        for ( size_t j = 0; j < sit->n_good_lambda; ++j ) {

          size_t global_spec = sit->first_good_spec + i;
          size_t global_lambda = sit->first_good_lambda + j;

          coverage ( global_spec, global_lambda ) += 1.0;

        }

      }

    }
    
  }

  // check that every spectral point is covered once

  for ( size_t i = 0; i < nspec; ++i ) {
    for ( size_t j = 0; j < nlambda; ++j ) {

      if ( ( coverage(i,j) < 0.5 ) || ( coverage(i,j) > 1.5 ) ) {
        ostringstream o;
        o << "FAIL:  spectrum " << i << ", lambda " << j << " was assigned to " << coverage(i,j) << " workers";
        cerr << o.str() << endl;
        exit(1);
      }

    }
  }

  cout << "  (PASSED)" << endl;
     
  return;
}



