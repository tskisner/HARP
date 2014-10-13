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


string top_prefix = "harp:  ";



void check_psf ( string const & prefix, psf_p chk_psf ) {

  size_t imgrows = chk_psf->img_rows();
  size_t imgcols = chk_psf->img_cols();
  size_t npix = imgrows * imgcols;

  size_t nspec = chk_psf->n_spec();
  size_t nlambda = chk_psf->n_lambda();
  size_t nbins = nspec * nlambda;

  vector_double lambda = chk_psf->lambda();

  cout << prefix << "Checking PSF of type " << chk_psf->type() << endl;
  cout << prefix << "  " << nspec << " spectra each with " << nlambda << " wavelength points" << endl;
  cout << prefix << "    (" << nbins << " total)" << endl;
  cout << prefix << "    lambda = " << lambda[0] << " ... " << lambda[ lambda.size() - 1 ] << endl;
  cout << prefix << "  image dimensions = " << imgcols << " x " << imgrows << " pixels" << endl;
  cout << prefix << "    (" << npix << " total)" << endl;

  return;
}


void check_spec ( string const & prefix, spec_p chk_spec ) {

  size_t nspec = chk_spec->n_spec();
  size_t nlambda = chk_spec->n_lambda();
  size_t nbins = nspec * nlambda;

  vector_double lambda;
  chk_spec->lambda( lambda );

  cout << prefix << "Checking spec of type " << chk_spec->type() << endl;
  cout << prefix << "  " << nspec << " spectra each with " << nlambda << " wavelength points" << endl;
  cout << prefix << "    (" << nbins << " total)" << endl;
  cout << prefix << "    lambda = " << lambda[0] << " ... " << lambda[ lambda.size() - 1 ] << endl;

  return;
}


void check_image ( string const & prefix, image_p chk_image ) {

  size_t imgrows = chk_image->n_rows();
  size_t imgcols = chk_image->n_cols();
  size_t npix = imgrows * imgcols;

  cout << prefix << "Checking image of type " << chk_image->type() << endl;
  cout << prefix << "  dimensions = " << imgcols << " x " << imgrows << " pixels" << endl;
  cout << prefix << "    (" << npix << " total)" << endl;

  return;
}


void check_targets ( string const & prefix, targets_p chk_targets ) {

  size_t nobj = chk_targets->n_objects();

  vector < object_p > objs = chk_targets->objects();

  size_t nsky = 0;

  for ( size_t i = 0; i < nobj; ++i ) {
    if ( objs[i]->is_sky() ) {
      ++nsky;
    }
  }

  cout << prefix << "Checking targets of type " << chk_targets->type() << endl;
  cout << prefix << "  " << nobj << " total objects, " << nsky << " are sky" << endl;

  return;
}


void check_group ( string const & prefix, group_p chk_group ) {

  cout << prefix << "Checking group" << endl;

  string new_prefix = prefix + "  ";

  check_psf ( new_prefix, chk_group->psf() );

  for ( list < image_p > :: iterator it = chk_group->images().begin(); it != chk_group->images().end(); ++it ) {
    check_image ( new_prefix, (*it) );
  }

  return;
}




int main ( int argc, char *argv[] ) {

  cout.precision ( 10 );
  cerr.precision ( 10 );

  string jsonpar = "";
  
  // Declare options
  
  popts::options_description desc ( "Allowed Options" );
  
  desc.add_options()
  ( "help,h", "display usage information" )
  ( "debug,d", "write out intermediate data products for debugging" )
  ;

  popts::options_description hidden ( "Hidden Options" );

  hidden.add_options()
    ( "json", popts::value < string > ( &jsonpar ) )
    ;

  popts::options_description all ( "All Options" );

  all.add(desc).add(hidden);

  popts::positional_options_description p;
  p.add("json", -1);

  popts::variables_map vm;
  popts::store(popts::command_line_parser( argc, argv ).options(all).positional(p).run(), vm);
  popts::notify(vm);

  if ( ( argc < 2 ) || vm.count( "help" ) ) {
    cerr << endl;
    cerr << "Usage: " << argv[0] << " [options] [JSON file]" << endl;
    cerr << endl;
    cerr << desc << endl;
    return 0;
  }
  
  // Get plugin registry

  plugin_registry & reg = plugin_registry::get();

  // Read metadata
  
  boost::property_tree::ptree tree;

  if ( vm.count( "json" ) ) {
    cout << top_prefix << "Loading JSON file..." << endl;
    boost::property_tree::json_parser::read_json ( jsonpar, tree );
  } else {
    cerr << "you must specify a JSON document to validate" << endl;
    exit(0);
  }

  // iterate over all recognized JSON keys and check

  boost::property_tree::ptree::const_iterator v = tree.begin();

  while ( v != tree.end() ) {

    if ( v->first == "psf" ) {

      psf_p chk_psf = load_psf ( v->second );
      check_psf ( top_prefix, chk_psf );

    } else if ( v->first == "image" ) {

      image_p chk_image = load_image ( v->second );
      check_image ( top_prefix, chk_image );

    } else if ( v->first == "spec" ) {

      spec_p chk_spec = load_spec ( v->second );
      check_spec ( top_prefix, chk_spec );
    
    } else if ( v->first == "targets" ) {

      targets_p chk_targets = load_targets ( v->second );
      check_targets ( top_prefix, chk_targets );
    
    } else if ( v->first == "group" ) {

      group_p chk_group ( new group ( v->second ) );
      check_group ( top_prefix, chk_group );
    
    } else {
      cout << "Skipping unrecognized object \"" << v->first << "\"" << endl;
    }

    ++v;
  }

  return 0;
}


