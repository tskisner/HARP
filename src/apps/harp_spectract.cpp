// @COPYRIGHT@

#include <iostream>
#include <cstdio>

#include <boost/program_options.hpp>

#include <harp.hpp>


using namespace std;
using namespace harp;


namespace popts = boost::program_options;


void parse_format ( string & input, string & format, map < string, string > & params ) {
  
  params.clear();
  
  size_t stroff = 0;
  size_t strend;
  
  string formatsep = ":";
  string paramsep = ",";
  string keysep = "=";
  
  // read format type
  
  strend = input.find ( formatsep, stroff );
  if ( strend == string::npos ) {
    format = input;
    cerr << "dbg:  format = " << format << endl;
    return;
  } 
  
  format = input.substr ( stroff, strend - stroff );
  stroff = strend + 1;
  
  if ( input.find ( formatsep, stroff ) != string::npos ) {
    MOAT_THROW( "only one format type allowed in format string" );
  }
  
  // read format params
  
  size_t parend;
  
  strend = 0;
  while ( strend != string::npos ) {
    strend = input.find ( paramsep, stroff );
    if ( strend == string::npos ) {
      parend = input.size();
    } else {
      parend = strend;
    }
    
    // parse this key/value
    
    string key;
    string val;
    string par = input.substr ( stroff, parend - stroff );
    
    cerr << "dbg:  parsing parameter \"" << par << "\"" << endl;
    
    size_t keyend = par.find ( keysep );
    
    if ( keyend == string::npos ) {
      key = par;
      val = "TRUE";
    } else {
      key = par.substr ( 0, keyend );
      val = par.substr ( keyend + 1, par.size() - (keyend + 1) );
    }
    
    params[ key ] = val;
    
    stroff = parend + 1;
  }
  
  return;
}


int main ( int argc, char *argv[] ) {
  
  cout.precision ( 16 );
  cerr.precision ( 16 );
  
  string imageformat = "";
  string psfformat = "";
  string psffile = "";
  string outfile = "";
  vector < string > imagefiles;
  
  bool quiet = false;  
  
  
  // Declare options
  
  popts::options_description desc ( "Allowed Options" );
  
  desc.add_options()
  ( "help,h", "Display usage information" )
  ( "psfformat", popts::value < string > ( &psfformat ), "format parameters of the PSF file" )
  ( "psf", popts::value < string > ( &psffile ), "path to PSF file" )
  ( "imageformat", popts::value < string > ( &imageformat ), "format parameters of the image files" )
  ( "output", popts::value < string > ( &outfile ), "path to output spectra file" )
  ( "quiet,q", "supress information printing" )
  ( "image-files", popts::value < vector < string > >(&imagefiles), "input image files" )
  ;
  
  popts::options_description hidden;
  hidden.add_options()
  ("positional", popts::value < vector < string > >())
  ;

  popts::options_description all_options;
  all_options.add(desc).add(hidden);
  
  popts::positional_options_description posargs;
  posargs.add("positional", -1);
  
  popts::variables_map vm;

  popts::store(popts::command_line_parser(argc, argv).options(all_options).positional(posargs).run(), vm);
  popts::notify(vm);


  if ( ( argc < 2 ) || vm.count( "help" ) ) {
    cerr << endl;
    cerr << desc << endl;
    return 0;
  }
  
  if ( vm.count( "positional" ) ) {
    imagefiles = vm["positional"].as < vector < string > > ();
  }
  
  if ( vm.count( "quiet" ) ) {
    quiet = true;
  }
  
  
  
  if ( ! quiet ) {
    cerr << "harp_spectract:  Reading input PSF data...           ";
  }
  
  string psftype;
  map < string, string > psfparams;
  parse_format ( psfformat, psftype, psfparams );
  
  psfparams[ "path" ] = psffile;
  
  //psf_p resp ( psf::create ( psfformat, params ) );
  
  //size_t nspec = resp->nspec();
  //size_t specbins = resp->specsize(0);
  //size_t nbins = nspec * specbins;


  string imagetype;
  map < string, string > imageparams;
  parse_format ( imageformat, imagetype, imageparams );

  
  return 0;
}


