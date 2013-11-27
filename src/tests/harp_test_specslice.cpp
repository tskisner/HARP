#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include <harp_test.hpp>


using namespace std;
using namespace harp;


void harp::test_specslice ( string const & datadir ) {

  cerr << "Testing spec slicing..." << endl;

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

  for ( size_t i = 0; i < nchunk; ++i ) {
    size_t lambda_rank = (size_t)( i / nchunk_spec );
    size_t spec_rank = i - lambda_rank * nchunk_spec;

    



  }


  cerr << "  (PASSED)" << endl;
     
  return;
}



