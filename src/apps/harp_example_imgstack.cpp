// @COPYRIGHT@

/*

  This toy example implements an image "calculator" using a stack (think RPN calculator).  The input JSON
  document specifies images and operators.  After running "make check", there will be some test data
  files created in HARP/src/tests/testdata/.  If you "cd" into that directory, you can test this example
  by running:

  $>  harp_example_imgstack imgstack.json

  NOTE:  this is just an example, and does not have several checks on the input JSON document that would
  be required for production code.  For example, this code does not check that the first operation is a
  "PUSH" or that the stack has a single image remaining at the end.

*/


#include <iostream>
#include <cstdio>

#include <boost/program_options.hpp>

#include <harp.hpp>


using namespace std;
using namespace harp;

namespace popts = boost::program_options;


int main ( int argc, char *argv[] ) {

  double tstart;
  double tstop;

  cout.precision ( 12 );
  cerr.precision ( 12 );

  string jsonpar = "";
  string outfile = "imgstack.fits.out";
  
  // Parse commandline options
  
  popts::options_description desc ( "Allowed Options" );
  
  desc.add_options()
  ( "help,h", "display usage information" )
  ( "out", popts::value<string>( &outfile ), "output image file" )
  ( "par", popts::value<string>( &jsonpar ), "JSON parameter file" )
  ;

  popts::variables_map vm;

  popts::store(popts::command_line_parser( argc, argv ).options(desc).run(), vm);
  
  popts::notify(vm);

  if ( ( argc < 2 ) || vm.count( "help" ) || ( ! vm.count( "par" ) ) ) {
    cerr << endl;
    cerr << desc << endl;
    return 0;
  }

  // Read JSON into a property tree
  
  boost::property_tree::ptree params;
  boost::property_tree::json_parser::read_json ( jsonpar, params );

  // This is our stack

  std::deque < vector_double > stack;

  size_t checkrows = 0;
  size_t checkcols = 0;

  // iterate over the input images and operators

  boost::property_tree::ptree::const_iterator v = params.begin();

  while ( v != params.end() ) {

    if ( v->first == "PUSH" ) {
      // we are pushing a new image onto the stack.  we don't care what format it is
      // so we use the factory method to instantiate the image from parameters specified
      // in the JSON.

      image_p img ( image::create ( v->second ) );

      if ( checkrows == 0 ) {
        // this is the first image

        checkrows = img->n_rows();
        checkcols = img->n_cols();
      } else {
        // verify that the image dimensions are consistent

        if ( ( checkrows != img->n_rows() ) || ( checkcols != img->n_cols() ) ) {
          HARP_THROW( "inconsistent image dimensions on the stack" );
        }
      }

      stack.push_front ( vector_double() );
      stack[0].resize ( checkrows * checkcols );

      img->values ( stack[0] );

      // img goes out of scope here and the raw pointer it wraps is deleted.

    } else if ( v->first == "ADD" ) {

      stack[1] += stack[0];
      stack.pop_front();

    } else if ( v->first == "SUB" ) {

      stack[1] -= stack[0];
      stack.pop_front();

    } else if ( v->first == "MUL" ) {

      stack[1] = boost::numeric::ublas::element_prod ( stack[1], stack[0] );
      stack.pop_front();

    } else if ( v->first == "DIV" ) {

      stack[1] = boost::numeric::ublas::element_div ( stack[1], stack[0] );
      stack.pop_front();

    } else {

      HARP_THROW( "undefined stack operator" );

    }

    ++v;

  }

  // write out the output image.  when writing data, we know exactly what format we are
  // writing and do not use the factory technique- we just instantiate the class directly.

  boost::property_tree::ptree img_props;
  img_props.put ( "format", "fits" );
  img_props.put < size_t > ( "rows", checkrows );
  img_props.put < size_t > ( "cols", checkcols );

  image_fits outimg ( img_props );

  // for this example, we don't care about the inverse variance of each image (we ignored it
  // above).  here we just write fake data.

  vector_double fake_invvar ( checkrows * checkcols );
  for ( size_t i = 0; i < fake_invvar.size(); ++i ) {
    fake_invvar[i] = 1.0;
  }
  std::vector < bool > fake_sky(1);

  outimg.write ( outfile, stack[0], fake_invvar, fake_sky );

  return 0;
}



