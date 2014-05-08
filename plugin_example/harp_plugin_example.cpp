// @COPYRIGHT@

// This header must be included by external plugins
#ifdef USE_MPI
#  include <harp/mpi_plugin.hpp>
#else
#  include <harp/plugin.hpp>
#endif

#include <harp_example.hpp>

using namespace std;
using namespace harp;


// example spec
//-----------------

harp::spec_example::spec_example ( boost::property_tree::ptree const & props ) : spec ( "example", props ) {
  nspec_ = 1;
  nlambda_ = 1;  
}

// This export statement allows use of the plugin with serialization operations in HARP.

BOOST_CLASS_EXPORT(harp::spec_example)

// Define the plugin creation function here.

spec * harp::spec_example_create ( boost::property_tree::ptree const & props ) {
  return new spec_example ( props );
}


// example psf
//-----------------

harp::psf_example::psf_example ( boost::property_tree::ptree const & props ) : psf ( "example", props ) {
  nspec_ = 1;
  nlambda_ = 1;
  rows_ = 1;
  cols_ = 1;
}

// This export statement allows use of the plugin with serialization operations in HARP.

BOOST_CLASS_EXPORT(harp::psf_example)

// Define the plugin creation function here.

psf * harp::psf_example_create ( boost::property_tree::ptree const & props ) {
  return new psf_example ( props );
}


// example image
//-----------------

harp::image_example::image_example ( boost::property_tree::ptree const & props ) : image ( "example", props ) {
  rows_ = 1;
  cols_ = 1;
}

// This export statement allows use of the plugin with serialization operations in HARP.

BOOST_CLASS_EXPORT(harp::image_example)

// Define the plugin creation function here.

image * harp::image_example_create ( boost::property_tree::ptree const & props ) {
  return new image_example ( props );
}


// Define the initialize function that will be called by the plugin registry.  There must
// be only one such function per *.so file, but you can register multiple plugins within
// this function call.

void initialize ( void * registry ) {

  plugin_registry * reg = static_cast < plugin_registry * > ( registry );

  string const & version = harp::source_version();

  // register spec plugin
  reg->register_spec ( "example", spec_example_create, version );

  // register psf plugin
  reg->register_psf ( "example", psf_example_create, version );

  // register image plugin
  reg->register_image ( "example", image_example_create, version );

  return;
}



