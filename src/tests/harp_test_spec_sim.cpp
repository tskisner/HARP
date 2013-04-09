#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include <harp_test.hpp>


extern "C" {
#include <unistd.h>
}

using namespace std;
using namespace harp;


void harp::test_spec_sim ( string const & datadir ) {

  int np;
  int myp;
  int ret;

  MPI_Comm_size ( MPI_COMM_WORLD, &np );
  MPI_Comm_rank ( MPI_COMM_WORLD, &myp );

  if ( myp == 0 ) {
    cerr << "Testing simulated spec generation..." << endl;
  }

  size_t nlambda = 50;
  size_t nspec = 100;
  size_t first_lambda = 8000.0;
  size_t last_lambda = 8001.1;

  size_t nbins = nspec * nlambda;

  boost::property_tree::ptree spec_props;
  spec_props.clear();
  spec_props.put ( "format", "sim" );
  spec_props.put ( "nspec", nspec );
  spec_props.put ( "nlambda", nlambda );
  spec_props.put ( "first_lambda", first_lambda );
  spec_props.put ( "last_lambda", last_lambda );
  spec_props.put ( "back", 10.0 );
  spec_props.put ( "atm", 500.0 );
  spec_props.put ( "obj", 80.0 );
  spec_props.put ( "atmspace", 12 );
  spec_props.put ( "skymod", 25 );
  spec_p testspec ( spec::create ( spec_props ) );

  if ( nspec != testspec->nspec() ) {
    cerr << "simulated spec nspec (" << testspec->nspec() << ") does not match input (" << nspec << ")" << endl;
  }

  if ( nlambda != testspec->nlambda() ) {
    cerr << "simulated spec nlambda (" << testspec->nlambda() << ") does not match input (" << nlambda << ")" << endl;
  }

  matrix_dist data ( nbins, 1 );
  vector < double > lambda;
  vector < bool > sky;

  testspec->read( data, lambda, sky );

  string outfile = datadir + "/spec_sim_data.out";
  data.Write( outfile );

  fstream fout;
  fout.precision(16);

  outfile = datadir + "/spec_sim_lambda.out";
  fout.open ( outfile.c_str(), ios::out );
  for ( size_t i = 0; i < nlambda; ++i ) {
    fout << lambda[i] << endl;
  }
  fout.close();

  outfile = datadir + "/spec_sim_sky.out";
  fout.open ( outfile.c_str(), ios::out );
  for ( size_t i = 0; i < nspec; ++i ) {
    fout << (int)sky[i] << endl;
  }
  fout.close();

  if ( myp == 0 ) {
    cerr << "  (PASSED)" << endl;
  }
     
  return;
}



