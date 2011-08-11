// @COPYRIGHT@

#include <iostream>
#include <cstdio>

#include <boost/program_options.hpp>

#include <harp.hpp>


using namespace std;
using namespace harp;


namespace popts = boost::program_options;


static bool harp_specstract_quiet = false;


void harp_specstract_prec ( data_vec & in, data_vec & out, int_vec & flags, void * data ) {
  data_vec * prec = (data_vec *) data;
  
  size_t vit;
  size_t n = in.size();
  
  #ifdef _OPENMP
  #pragma omp parallel for default(none) private(vit) shared(n, flags, in, out, prec) schedule(static)
  #endif
  for ( vit = 0; vit < n; ++vit ) {
    if ( flags[ vit ] == 0 ) {
      out[ vit ] = (*prec)[ vit ] * in[ vit ];
    } else {
      out[ vit ] = 0.0;
    }
  }
  
  return;
}


void harp_specstract_report ( double const & norm, double const & deltazero, int const & iter, double const & alpha, double const & beta, double const & delta, double const & epsilon ) {
  
  double relerr = sqrt ( delta / deltazero );
  
  if ( ! harp_specstract_quiet ) {
    cout << "harp_specstract:  PCG iter " << iter << ": epsilon = " << epsilon << " relative err = " << relerr << endl;
  }
  
  return;
}

void harp_specstract_profile ( string const & name, string const & desc, double & totaltime, double & opencltime, map < string, long long int > & papi ) {
  
  cout << "harp_specstract:   " << desc << ":  " << totaltime << " seconds" << endl;
  
  return;
}


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
    //cerr << "dbg:  format = " << format << endl;
    return;
  } 
  
  format = input.substr ( stroff, strend - stroff );
  //cerr << "dbg:  format = " << format << endl;
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
    
    //cerr << "dbg:  parsing parameter \"" << par << "\"" << endl;
    
    size_t keyend = par.find ( keysep );
    
    if ( keyend == string::npos ) {
      key = par;
      val = "TRUE";
    } else {
      key = par.substr ( 0, keyend );
      val = par.substr ( keyend + 1, par.size() - (keyend + 1) );
    }
    
    params[ key ] = val;
    //cerr << "dbg:  key " << key << " = " << params[key] << endl;
    
    stroff = parend + 1;
  }
  
  return;
}


typedef struct {
  image_p handle;
  string path;
  string format;
  string type;
  map < string, string > params;
  size_t rows;
  size_t cols;
  size_t npix;
} harp_image_props;


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
  ;
  
  popts::options_description hidden;
  hidden.add_options()
  ("positional", popts::value < vector < string > >(&imagefiles))
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
    cerr << "  <image file 1> [ <image file 2> ... ]" << endl << endl; 
    return 0;
  }
  
  if ( vm.count( "positional" ) ) {
    imagefiles = vm["positional"].as < vector < string > > ();
  } else {
    cerr << "you must specify at least one image file" << endl;
    return 0;
  }
  
  if ( vm.count( "quiet" ) ) {
    quiet = true;
    harp_specstract_quiet = true;
  }
  
  moat::profile * prof = moat::profile::get ( );
  
  if ( ! quiet ) {

    prof->reg ( "HARP_PSF", "compute projection matrix" );
    prof->reg ( "HARP_REMAP", "remap projection matrix" );
    prof->reg ( "HARP_PRECOND", "compute preconditioner" );
    prof->reg ( "HARP_READ", "read data" );
    prof->reg ( "HARP_WRITE", "write data" );

    prof->reg ( "HARP_PCG_PREC", "applying preconditioner" );
    prof->reg ( "HARP_PCG_PMV", "projection matrix-vector multiply" );
    prof->reg ( "HARP_PCG_NMV", "N^-1 matrix-vector multiply" );
    prof->reg ( "HARP_PCG_VEC", "vector ops time" );
    prof->reg ( "HARP_PCG_TOT", "Total PCG time" );
    
  }
  
  
  size_t nimages = imagefiles.size();
  
  vector < harp_image_props > images ( nimages );
  
  
  if ( ! quiet ) {
    cout << "harp_specstract:  Reading input PSF properties...              ";
  }
  
  string psftype;
  map < string, string > psfparams;
  
  if ( psfformat == "" ) {
    psftype = "boss";
  } else {
    parse_format ( psfformat, psftype, psfparams );
  }
  psfparams[ "path" ] = psffile;
  
  psf_p resp ( psf::create ( psftype, psfparams ) );
  
  size_t nspec = resp->nspec();
  size_t specbins = resp->specsize(0);
  size_t nbins = nspec * specbins;
  
  if ( ! quiet ) {
    cout << "DONE" << endl;
  }

  
  if ( ! quiet ) {
    cout << "harp_specstract:  Reading input image properties...            ";
  }

  string imagetype;
  map < string, string > imageparams;
  parse_format ( imageformat, imagetype, imageparams );
  
  for ( size_t i = 0; i < nimages; ++i ) {
    images[i].format = imageformat;
    images[i].type = imagetype;
    images[i].params = imageparams;
    
    images[i].path = imagefiles[i];
    
    images[i].params[ "path" ] = images[i].path;
    
    images[i].handle = image_p ( image::create ( images[i].type, images[i].params ) );
    
    images[i].rows = images[i].handle->rows();
    images[i].cols = images[i].handle->cols();
    images[i].npix = images[i].rows * images[i].cols;
  }
  
  if ( ! quiet ) {
    cout << "DONE" << endl;
  }
  
  
  if ( ! quiet ) {
    cout << "harp_specstract:  Reading input image and noise covariance...  ";
    prof->start ( "HARP_READ" );
  }

  // FIXME:  for now, use only the first image
  
  data_vec imgdata ( images[0].npix );
  data_vec_view imgdataview ( imgdata, mv_range ( 0, images[0].npix ) );
  
  data_vec imgnoise ( images[0].npix );
  data_vec_view imgnoiseview ( imgnoise, mv_range ( 0, images[0].npix ) );
  
  images[0].handle->read ( imgdataview );
  images[0].handle->read_noise ( imgnoiseview );
  
  
  if ( ! quiet ) {
    prof->stop ( "HARP_READ" );
    cout << "DONE" << endl;
  }
  
  
  if ( ! quiet ) {
    cout << "harp_specstract:  Computing PSF...                             ";
  }
  
  // FIXME:  probably we will have pairs of images and PSFs...
  
  comp_rowmat projmat ( images[0].npix, nbins );
  
  if ( quiet ) {
    resp->projection ( string(""), string(""), 0, nspec - 1, 0, specbins - 1, 0, images[0].cols - 1, 0, images[0].rows - 1, projmat );
  } else {
    resp->projection ( string("HARP_PSF"), string("HARP_REMAP"), 0, nspec - 1, 0, specbins - 1, 0, images[0].cols - 1, 0, images[0].rows - 1, projmat );
  }
  
  if ( ! quiet ) {
    cout << "DONE" << endl;
  }
  
  
  if ( ! quiet ) {
    cout << "harp_specstract:  Computing preconditioner...                  ";
    prof->start ( "HARP_PRECOND" );
  }
  
  data_vec outspec ( nbins );
  int_vec flags ( nbins );
  
  for ( size_t b = 0; b < nbins; ++b ) {
    flags[b] = 0;
    outspec[b] = 0.0;
  }
  
  comp_rowmat invnoise ( images[0].npix, images[0].npix );
  data_vec precdata ( nbins );
  
  for ( size_t i = 0; i < images[0].npix; ++i ) {
    invnoise( i, i ) = imgnoise[i];
  }
  
  for ( size_t i = 0; i < nbins; ++i ) {
    precdata[i] = 0.0;
    for ( size_t j = 0; j < images[0].npix; ++j ) {
      precdata[i] += projmat( j, i ) * projmat( j, i ) * invnoise( j, j );
    }
  }
  
  for ( size_t i = 0; i < nbins; ++i ) {
    precdata[i] = 1.0 / precdata[i];
  }
  
  if ( ! quiet ) {
    cout << "DONE" << endl;
    prof->stop ( "HARP_PRECOND" );
  }
  

  if ( ! quiet ) {
    cout << "harp_specstract:  Solving PCG..." << endl;
  }
  
  data_vec rhs ( nbins );
  data_vec q ( nbins );
  data_vec r ( nbins );
  data_vec s ( nbins );
  data_vec d ( nbins );
  
  double err;
  
  if ( quiet ) {
    err = moat::la::pcg_mle < comp_rowmat, comp_rowmat, data_vec, int_vec > ( true, true, projmat, invnoise, imgdata, outspec, q, r, s, d, flags, rhs, 100, 1.0e-12, harp_specstract_prec, (void*)&precdata, harp_specstract_report, "", "", "", "", "" );
  } else {
    err = moat::la::pcg_mle < comp_rowmat, comp_rowmat, data_vec, int_vec > ( true, true, projmat, invnoise, imgdata, outspec, q, r, s, d, flags, rhs, 100, 1.0e-12, harp_specstract_prec, (void*)&precdata, harp_specstract_report, "HARP_PCG_TOT", "HARP_PCG_VEC", "HARP_PCG_PMV", "HARP_PCG_NMV", "HARP_PCG_PREC" );
  }
  
  
  if ( ! quiet ) {
    cout << "harp_specstract:  Writing output solved spectra...             ";
    prof->start ( "HARP_WRITE" );
  }
  
  string rmcom = "rm -f " + outfile;
  
  system( rmcom.c_str() );
  
  map < string, string > specparams;
  specparams[ "hdu" ] = "1";
  specparams[ "nspec" ] = nspec;
  specparams[ "specsize" ] = specbins;
  
  spec_p solvespec ( spec::create ( string("boss"), specparams ) );
  data_vec_view solvespecview ( outspec, mv_range ( 0, nspec * specbins ) );
  
  solvespec->write ( outfile, solvespecview );
  
  if ( ! quiet ) {
    cout << "DONE" << endl;
    prof->stop ( "HARP_WRITE" );
  }

  if ( ! quiet ) {
    
    cout << "harp_specstract:  Profiling information:" << endl;
    
    prof->stop_all();
    
    prof->query ( harp_specstract_profile );

    prof->unreg ( "HARP_PSF" );
    prof->unreg ( "HARP_REMAP" );
    prof->unreg ( "HARP_PRECOND" );
    prof->unreg ( "HARP_READ" );
    prof->unreg ( "HARP_WRITE" );

    prof->unreg ( "HARP_PCG_PREC" );
    prof->unreg ( "HARP_PCG_PMV" );
    prof->unreg ( "HARP_PCG_NMV" );
    prof->unreg ( "HARP_PCG_VEC" );
    prof->unreg ( "HARP_PCG_TOT" );
    
  }
  
  return 0;
}


