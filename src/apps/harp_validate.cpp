// @COPYRIGHT@

#include <iostream>
#include <cstdio>

extern "C" {
  #include <unistd.h>
}

#include <boost/program_options.hpp>

#include <harp.hpp>


using namespace std;
using namespace harp;

namespace popts = boost::program_options;


string prefix = "harp:  ";


void parse_opt ( string const & optstr, string & key, string & value ) {
  size_t pos = 0;
  size_t sep = optstr.find ( '=', pos );

  if ( pos == string::npos ) {
    cerr << "cannot parse key / value" << endl;
    exit(1);
  }

  size_t len = sep - pos;

  key = optstr.substr ( pos, len );

  pos = sep + 1;
  len = optstr.size() - pos;

  value = optstr.substr ( pos, len );

  return;
}


void parse_propstr ( string const & propstr, boost::property_tree::ptree & proptree, string & jsonpath, string & jsonclass ) {

  proptree.clear();

  boost::property_tree::ptree objprops;
  objprops.clear();

  size_t pos = 0;

  size_t sep = propstr.find ( ':', pos );

  if ( pos == string::npos ) {
    cerr << "cannot parse object class" << endl;
    exit(1);
  }

  size_t len = sep - pos;

  jsonclass = propstr.substr ( pos, len );
  jsonpath = jsonclass;

  pos = sep + 1;

  sep = propstr.find ( ':', pos );

  if ( pos == string::npos ) {
    cerr << "cannot parse object type" << endl;
    exit(1);
  }

  len = sep - pos;

  string jsontype = propstr.substr ( pos, len );

  pos = sep + 1;

  len = propstr.size() - pos;

  string parstr = propstr.substr ( pos, len );

  string typestr = jsonclass + "_type";
  proptree.put ( typestr, jsontype );

  if ( parstr.size() != 0 ) {

    pos = 0;
    sep = parstr.find ( ',', pos );

    string optstr;
    string key;
    string val;

    while ( sep != string::npos ) {
      len = sep - pos;

      optstr = parstr.substr ( pos, len );
      parse_opt ( optstr, key, val );
      objprops.put ( key, val );

      pos = sep + 1;
      sep = parstr.find ( ',', pos );
    }

    len = parstr.size() - pos;

    optstr = parstr.substr ( pos, len );
    parse_opt ( optstr, key, val );
    objprops.put ( key, val );

  }

  proptree.put_child ( jsonclass, objprops );

  return;
}


void check_psf ( string const & type, boost::property_tree::ptree & proptree ) {

  cout << prefix << "Creating PSF of type " << type << endl;

  plugin_registry & reg = plugin_registry::get();

  psf_p design ( reg.create_psf ( type, proptree ) );

  size_t imgrows = design->img_rows();
  size_t imgcols = design->img_cols();
  size_t npix = imgrows * imgcols;

  size_t nspec = design->n_spec();
  size_t nlambda = design->n_lambda();
  size_t nbins = nspec * nlambda;

  vector_double lambda = design->lambda();

  cout << prefix << "  " << nspec << " spectra each with " << nlambda << " wavelength points" << endl;
  cout << prefix << "    (" << nbins << " total)" << endl;
  cout << prefix << "    lambda = " << lambda[0] << " ... " << lambda[ lambda.size() - 1 ] << endl;
  cout << prefix << "  image dimensions = " << imgcols << " x " << imgrows << " pixels" << endl;
  cout << prefix << "    (" << npix << " total)" << endl;

  return;
}


void check_spec ( string const & type, boost::property_tree::ptree & proptree ) {

  cout << prefix << "Creating spec of type " << type << endl;

  plugin_registry & reg = plugin_registry::get();

  spec_p sp ( reg.create_spec ( type, proptree ) );

  size_t nspec = sp->n_spec();
  size_t nlambda = sp->n_lambda();
  size_t nbins = nspec * nlambda;

  vector_double lambda;
  sp->lambda( lambda );

  cout << prefix << "  " << nspec << " spectra each with " << nlambda << " wavelength points" << endl;
  cout << prefix << "    (" << nbins << " total)" << endl;
  cout << prefix << "    lambda = " << lambda[0] << " ... " << lambda[ lambda.size() - 1 ] << endl;

  return;
}


void check_image ( string const & type, boost::property_tree::ptree & proptree ) {

  cout << prefix << "Creating image of type " << type << endl;

  plugin_registry & reg = plugin_registry::get();

  image_p img ( reg.create_image ( type, proptree ) );

  size_t imgrows = img->n_rows();
  size_t imgcols = img->n_cols();
  size_t npix = imgrows * imgcols;

  cout << prefix << "  dimensions = " << imgcols << " x " << imgrows << " pixels" << endl;
  cout << prefix << "    (" << npix << " total)" << endl;

  return;
}


void check_targets ( string const & type, boost::property_tree::ptree & proptree ) {

  cout << prefix << "Creating targets of type " << type << endl;

  plugin_registry & reg = plugin_registry::get();

  targets_p trg ( reg.create_targets ( type, proptree ) );

  size_t nobj = trg->n_objects();

  vector < object_p > objs = trg->objects();

  size_t nsky = 0;

  for ( size_t i = 0; i < nobj; ++i ) {
    if ( objs[i]->is_sky() ) {
      ++nsky;
    }
  }

  cout << prefix << "  " << nobj << " total objects, " << nsky << " are sky" << endl;

  return;
}



int main ( int argc, char *argv[] ) {

  cout.precision ( 10 );
  cerr.precision ( 10 );

  string jsonpar = "";
  string jsonpath = "";
  string jsonclass = "";
  string propstr = "";
  
  // Declare options
  
  popts::options_description desc ( "Allowed Options" );
  
  desc.add_options()
  ( "help,h", "display usage information" )
  ( "debug,d", "write out intermediate data products for debugging" )
  ( "props", popts::value < string > ( &propstr ), "property string (see below)" )
  ( "json", popts::value < string > ( &jsonpar ), "JSON file" )
  ( "jsonpath", popts::value < string > ( &jsonpath ), "JSON object to check, as a path (foo/bar/my_object)" )
  ( "jsonclass", popts::value < string > ( &jsonclass ), "JSON object class (psf, image, spec, or targets)" )
  ;

  popts::variables_map vm;

  popts::store(popts::command_line_parser( argc, argv ).options(desc).run(), vm);
  
  popts::notify(vm);

  if ( ( argc < 2 ) || vm.count( "help" ) ) {
    cerr << endl;
    cerr << desc << endl;
    cerr << "  You can specify either a JSON document and the object in that document" << endl;
    cerr << "  to check, or you can manually specify properties of a single object" << endl;
    cerr << "  to validate.  For a JSON object named \"my_object\", we expect to find" << endl;
    cerr << "  an object called \"my_object_type\" at the same level of the document." << endl;
    cerr << "  In the case of specifying a property string, the format is:" << endl << endl;

    cerr << "  <object class>:<type>:<prop1>=<val1>,<prop2>=<val2>,<prop2>=<val2>" << endl << endl;
    
    cerr << "  where valid object classes are psf, spec, image, targets.  The type" << endl;
    cerr << "  is any valid built-in or dynamically loaded plugin name of that" << endl;
    cerr << "  object class." << endl << endl;

    cerr << "  NOTE: for objects with nested properties (such as the simulated image" << endl;
    cerr << "  type), you must specify this as a JSON document." << endl << endl;
    return 0;
  }
  
  // Get plugin registry

  plugin_registry & reg = plugin_registry::get();

  // Read metadata
  
  boost::property_tree::ptree proptree;

  if ( vm.count( "json" ) ) {
    cout << prefix << "Loading JSON file..." << endl;
    boost::property_tree::json_parser::read_json ( jsonpar, proptree );
  } else {
    cout << prefix << "Parsing parameter string..." << endl;
    parse_propstr ( propstr, proptree, jsonpath, jsonclass );
  }

  // seek to our object of interest

  size_t pos = 0;
  size_t len;
  size_t sep = jsonpath.find ( '/', pos );
  boost::property_tree::ptree parent = proptree;
  boost::property_tree::ptree cursor;

  while ( sep != string::npos ) {
    len = sep - pos;

    string substr = jsonpath.substr ( pos, len );

    cursor = parent.get_child ( substr );

    pos = sep + 1;
    sep = jsonpath.find ( '/', pos );

    parent = cursor;
  }

  // now we can get our object and its type

  len = jsonpath.size() - pos;
  string name = jsonpath.substr ( pos, len );
  string typestr = name + "_type";

  string type = parent.get < string > ( typestr );

  cursor = parent.get_child ( name );

  if ( jsonclass == "psf" ) {
    check_psf ( type, cursor );
  } else if ( jsonclass == "image" ) {
    check_image ( type, cursor );
  } else if ( jsonclass == "spec" ) {
    check_spec ( type, cursor );
  } else if ( jsonclass == "targets" ) {
    check_targets ( type, cursor );
  } else {
    cerr << "Invalid object class \"" << jsonclass << "\"" << endl;
  }

  return 0;
}


