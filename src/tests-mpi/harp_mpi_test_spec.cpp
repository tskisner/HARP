#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include <harp_mpi_test.hpp>


extern "C" {
#include <unistd.h>
}

using namespace std;
using namespace harp;


void harp::mpi_test_spec ( string const & datadir ) {

  boost::mpi::communicator comm;

  int np = comm.size();
  int myp = comm.rank();

  if ( myp == 0 ) {
    cout << "Testing spec reading..." << endl;
  }

  // file from serial tests that we are reading in

  string specfile = datadir + "/spec_sim.fits.out";

  boost::property_tree::ptree spec_props;
  spec_props.clear();
  spec_props.put ( "path", specfile );

  // instantiate a serial version and read data

  size_t serial_nspec;
  size_t serial_nlambda;
  vector_double serial_data;
  vector_double serial_lambda;

  if ( myp == 0 ) {

    spec_fits serial_spec ( spec_props );

    serial_nspec = serial_spec.n_spec();
    serial_nlambda = serial_spec.n_lambda();
    serial_spec.values ( serial_data );
    serial_spec.lambda ( serial_lambda );

  }

  boost::mpi::broadcast ( comm, serial_nspec, 0 );
  boost::mpi::broadcast ( comm, serial_nlambda, 0 );
  boost::mpi::broadcast ( comm, serial_data, 0 );
  boost::mpi::broadcast ( comm, serial_lambda, 0 );

  // instantiate an MPI version and read

  mpi_spec_p par_spec ( new mpi_spec ( comm, "fits", spec_props ) );

  size_t par_nspec = par_spec->n_spec();
  size_t par_nlambda = par_spec->n_lambda();

  vector_double par_data;
  vector_double par_lambda;

  par_spec->values ( par_data );
  par_spec->lambda ( par_lambda );

  // verify

  if ( serial_nspec != par_nspec ) {
    cerr << "FAIL:  MPI nspec (" << par_nspec << ") does not match serial value (" << serial_nspec << ")" << endl;
    exit(1);
  }

  if ( serial_nlambda != par_nlambda ) {
    cerr << "FAIL:  MPI nlambda (" << par_nlambda << ") does not match serial value (" << serial_nlambda << ")" << endl;
    exit(1);
  }

  for ( size_t i = 0; i < serial_lambda.size(); ++i ) {
    if ( fabs ( ( serial_lambda[i] - par_lambda[i] ) / serial_lambda[i] ) > std::numeric_limits < float > :: epsilon() ) {
      cerr << "FAIL:  MPI lambda element " << i << " does not match serial value" << endl;
      exit(1);
    }
  }

  for ( size_t i = 0; i < serial_nspec; ++i ) {
    for ( size_t j = 0; j < serial_nlambda; ++j ) {
      double ser = serial_data[ i * serial_nlambda + j ];
      double par = par_data[ i * par_nlambda + j ];

      if ( fabs ( ( par - ser ) / ser ) > std::numeric_limits < float > :: epsilon() ) {
        cerr << "FAIL:  MPI spec " << i << ", bin " << j << " does not match serial value" << endl;
        exit(1);
      }
    }
  }

  if ( myp == 0 ) {
    cout << "  (PASSED)" << endl;
  }

     
  return;
}



