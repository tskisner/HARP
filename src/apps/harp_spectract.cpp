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
  }
  
  size_t nimages = imagefiles.size();
  
  vector < harp_image_props > images ( nimages );
  
  
  if ( ! quiet ) {
    cerr << "harp_spectract:  Reading input PSF properties...           ";
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
    cerr << "DONE" << endl;;
  }

  
  if ( ! quiet ) {
    cerr << "harp_spectract:  Reading input image properties...           ";
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
    cerr << "DONE" << endl;;
  }
  
  
  /*
  cerr << "  Reading pure signal image..." << endl;
  
  params.clear();
  
  params[ "path" ] = datadir + "/realtoy_input_image.fits";
  params[ "hdu" ] = "3";
  
  image_p sigimg ( image::create ( string("toy"), params ) );
  
  dense_rowmat sig ( sigimg->rows(), sigimg->cols() );
  
  dense_rowmat_view sigview ( sig, mv_range ( 0, sig.size1() ), mv_range ( 0, sig.size2() ) );
  
  sigimg->read ( 0, 0, sigview );
  
  data_vec rms ( npix );
  data_vec measured ( npix );
  
  size_t pixoff = 0;
  for ( size_t i = 0; i < rows; ++i ) {
    for ( size_t j = 0; j < cols; ++j ) {
      rms[i*cols + j] = sqrt( 16.0 + sig(i,j) );
      measured[pixoff] = img(i,j);
      ++pixoff;
    }
  }
  
  
  cerr << "  Reading noise covariance image..." << endl;
  
  params.clear();
  
  params[ "path" ] = datadir + "/realtoy_input_image.fits";
  params[ "hdu" ] = "2";
  
  image_p nseimg ( image::create ( string("toy"), params ) );
  
  dense_rowmat nse ( sigimg->rows(), sigimg->cols() );
  
  dense_rowmat_view nseview ( nse, mv_range ( 0, nse.size1() ), mv_range ( 0, nse.size2() ) );
  
  nseimg->read ( 0, 0, nseview );
  
  data_vec trueinv ( npix );
  
  pixoff = 0;
  for ( size_t i = 0; i < rows; ++i ) {
    for ( size_t j = 0; j < cols; ++j ) {
      trueinv[pixoff] = nse(i,j);
      ++pixoff;
    }
  }
  
  
  
  
  cerr << "  Testing PSF convolution..." << endl;
  
  data_vec compare ( npix );
  
  string outtruthfile = datadir + "/realtoy_compare_truth.out";
  
  out.open ( outtruthfile.c_str(), ios::out );
  
  for ( size_t i = 0; i < truth.size(); ++i ) {
    out << i << " " << truth(i) << endl;
  }

  out.close();
  
  boost::numeric::ublas::axpy_prod ( projmat, truth, compare, true );
  
  dense_rowmat tempmat ( rows, cols );
  
  dense_rowmat_view outview ( tempmat, mv_range ( 0, rows ), mv_range ( 0, cols ) );
  
  pixoff = 0;
  for ( size_t i = 0; i < rows; ++i ) {
    for ( size_t j = 0; j < cols; ++j ) {
      outview( i, j ) = compare[pixoff];
      ++pixoff;
    }
  }
  
  params.clear();
  
  params[ "hdu" ] = "1";
  
  ostringstream o;
  o << rows;
  params[ "rows" ] = o.str();
  
  o.str("");
  o << cols;
  params[ "cols" ] = o.str();
  
  image_p outsigimage ( image::create ( string("toy"), params ) );
  
  outsigimage->write ( "!" + datadir + "/realtoy_psf_convolved.fits.out", 0, 0, outview );
  
  string outimgfile = datadir + "/realtoy_compare_sigmap.out";
  
  out.open ( outimgfile.c_str(), ios::out );
  
  double ratio;
  pixoff = 0;
  for ( size_t i = 0; i < rows; ++i ) {
    for ( size_t j = 0; j < cols; ++j ) {
      ratio = compare[pixoff] - sig(i,j);
      out << pixoff << " " << i << " " << j << " " << sig(i,j) << " " << compare[pixoff] << " " << ratio << endl;
      ++pixoff;
    }
  }

  out.close();
  
  
  
  cerr << "  Computing preconditioner..." << endl;
  
  prof->start ( "PCG_PRECALC" );
  
  data_vec outspec ( nbins );
  int_vec flags ( nbins );
  
  for ( size_t b = 0; b < nbins; ++b ) {
    flags[b] = 0;
    outspec[b] = 0.0;
  }
  
  comp_rowmat invnoise ( npix, npix );
  data_vec precdata ( nbins );
  
  for ( size_t i = 0; i < npix; ++i ) {
    invnoise( i, i ) = trueinv[i];
    //invnoise( i, i ) = 1.0 / ( rms[i] * rms[i] );
    //cout << measured[i] << endl;
  }
  
  for ( size_t i = 0; i < nbins; ++i ) {
    precdata[i] = 0.0;
    for ( size_t j = 0; j < npix; ++j ) {
      precdata[i] += projmat( j, i ) * projmat( j, i ) * invnoise( j, j );
    }
  }
  
  for ( size_t i = 0; i < nbins; ++i ) {
    precdata[i] = 1.0 / precdata[i];
    //cout << precdata[i] << endl;
  }
  
  prof->stop ( "PCG_PRECALC" );
  
  
  cerr << "  Solving PCG..." << endl;
  
  data_vec rhs ( nbins );
  data_vec q ( nbins );
  data_vec r ( nbins );
  data_vec s ( nbins );
  data_vec d ( nbins );
  
  double err = moat::la::pcg_mle < comp_rowmat, comp_rowmat, data_vec, int_vec > ( true, true, projmat, invnoise, measured, outspec, q, r, s, d, flags, rhs, 100, 1.0e-12, realtoy_pcgmle_prec, (void*)&precdata, realtoy_pcgmle_report, "PCG_TOT", "PCG_VEC", "PCG_PMV", "PCG_NMV", "PCG_PREC" );
  
  prof->stop_all();
  
  prof->query ( realtoy_pcgmle_profile );
  
  prof->unreg ( "PCG_PREC" );
  prof->unreg ( "PCG_PMV" );
  prof->unreg ( "PCG_NMV" );
  prof->unreg ( "PCG_VEC" );
  prof->unreg ( "PCG_TOT" );
  prof->unreg ( "PCG_PSF" );
  prof->unreg ( "PCG_REMAP" );
  prof->unreg ( "PCG_PRECALC" );
  
  
  string outspecfile = datadir + "/realtoy_solved_spectra.fits.out";
  
  string rmcom = "rm -f " + outspecfile;
  
  system( rmcom.c_str() );
  
  params.clear();
  
  o.str("");
  o << specbins;

  fitsfile * fp;
  fits::open_readwrite ( fp, outspecfile );
  
  fits::img_append ( fp, nspec, specbins );
  
  fits::close ( fp );


  //params[ "path" ] = outspecfile;
  params[ "hdu" ] = "1";
  params[ "size" ] = o.str();
  
  for ( size_t s = 0; s < nspec; ++s ) {
  
    o.str("");
    o << s;
    params[ "pos" ] = o.str();
  
    spectrum_p solvespec ( spectrum::create ( string("toy"), params ) );
  
    data_vec_view solvespecview ( outspec, mv_range ( s * specbins, (s+1) * specbins ) );
  
    solvespec->write ( outspecfile, solvespecview );
  
  }
  
  
  outspecfile = datadir + "/realtoy_solved_spectra.out";
  
  out.open ( outspecfile.c_str(), ios::out );
  
  for ( size_t i = 0; i < nbins; ++i ) {
    if ( flags[i] == 0 ) {
      out << i << " " << truth[i] << " " << outspec[i] << " " << sqrt(precdata[i]) << endl;
    }
  }
  
  out.close();
*/
  
  return 0;
}


