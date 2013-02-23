// @COPYRIGHT@

#include <iostream>
#include <cstdio>

#include <boost/program_options.hpp>


#include <harp.hpp>


using namespace std;
using namespace harp;


namespace popts = boost::program_options;



int main ( int argc, char *argv[] ) {

  cliq::Initialize( argc, argv );

  int np;
  int myp;
  int ret;

  MPI_Comm_size ( MPI_COMM_WORLD, &np );
  MPI_Comm_rank ( MPI_COMM_WORLD, &myp );

  double tstart;
  double tstop;

  cout.precision ( 16 );
  cerr.precision ( 16 );

  size_t lambda_width = 3;
  size_t lambda_overlap = 1;
  
  string jsonconf = "";

  string prefix = "harp:  ";
  
  bool quiet = false;
  
  
  // Declare options
  
  popts::options_description desc ( "Allowed Options" );
  
  desc.add_options()
  ( "help,h", "Display usage information" )
  ( "quiet,q", "supress information printing" )
  ( "lambda_width", popts::value < size_t > ( &lambda_width ), "maximum wavelength points to process simultaneously" )
  ( "lambda_overlap", popts::value < size_t > ( &lambda_overlap ), "minimum wavelength points to overlap" )
  ( "conf", popts::value<string>(&jsonconf), "JSON configuration file" )
  ;

  popts::variables_map vm;
  
  popts::positional_options_description posargs;
  posargs.add("conf", -1);

  popts::store(popts::command_line_parser( argc, argv ).options(desc).positional(posargs).run(), vm);
  
  popts::notify(vm);

  if ( ( argc < 2 ) || vm.count( "help" ) || ( ! vm.count( "conf" ) ) ) {
    cerr << endl;
    cerr << desc << endl;
    cliq::Finalize();
    return 0;
  }

  if ( 2 * lambda_overlap > lambda_width ) {
    cerr << prefix << "wavelength overlap is more than half the wavelength band!" << endl;
    ret = MPI_Abort ( MPI_COMM_WORLD, 1 );
  }

  size_t lambda_core = lambda_width - 2 * lambda_overlap;
  
  if ( vm.count( "quiet" ) ) {
    quiet = true;
  }
  
  boost::property_tree::ptree conf;

  if ( ( myp == 0 ) && ( ! quiet ) ) {
    cout << prefix << "Loading JSON config file..." << endl;
  }

  boost::property_tree::json_parser::read_json ( jsonconf, conf );
  
  if ( ( myp == 0 ) && ( ! quiet ) ) {
    cout << prefix << "Reading input image..." << endl;
  }

  tstart = MPI_Wtime();
  
  boost::property_tree::ptree img_props;
  img_props = conf.get_child ( "image" );

  image_p img ( image::create ( img_props ) );

  size_t imgrows = img->rows();
  size_t imgcols = img->cols();
  size_t npix = imgrows * imgcols;

  matrix_local measured ( npix, 1 );  
  local_matrix_zero ( measured );

  matrix_local invnoise ( npix, 1 );
  local_matrix_zero ( invnoise );

  img->read ( measured );
  img->read_noise ( invnoise );

  tstop = MPI_Wtime();

  if ( ( myp == 0 ) && ( ! quiet ) ) {
    cout << prefix << "  time = " << tstop-tstart << " seconds" << endl;
    cout << prefix << "  dimensions = " << imgrows << " x " << imgcols << endl;
    cout << prefix << "Creating PSF..." << endl;
  }

  tstart = MPI_Wtime();

  boost::property_tree::ptree psf_props;
  psf_props = conf.get_child ( "psf" );

  psf_p epsf ( psf::create ( psf_props ) );

  size_t psf_imgrows = epsf->pixrows();
  size_t psf_imgcols = epsf->pixcols();
  size_t psf_npix = psf_imgrows * psf_imgcols;

  if ( psf_npix != npix ) {
    cerr << prefix << "PSF image dimensions (" << psf_imgrows << " x " << psf_imgcols << ") do not match image size" << endl;
    ret = MPI_Abort ( MPI_COMM_WORLD, 1 );
  }

  size_t nspec = epsf->nspec();
  size_t nlambda = epsf->nlambda();

  vector < double > lambda = epsf->lambda();

  size_t nbins = nspec * nlambda;

  tstop = MPI_Wtime();

  if ( ( myp == 0 ) && ( ! quiet ) ) {
    cout << prefix << "  time = " << tstop-tstart << " seconds" << endl;
    cout << prefix << "  image dimensions = " << psf_imgrows << " x " << psf_imgcols << endl;
    cout << prefix << "  " << nspec << " spectra with " << nlambda << " wavelength points each:" << endl;
    cout << prefix << "  lambda = " << lambda[0] << " ... " << lambda[nlambda - 1] << endl;

    cout << prefix << "Processing " << lambda_width << " wavelength points with overlap of " << lambda_overlap << " :" << endl;
  }

  vector < size_t > band_start;
  vector < size_t > band_stop;

  size_t offset = 0;

  while ( offset + lambda_width < nlambda ) {
    band_start.push_back( offset );
    band_stop.push_back( offset + lambda_width );
    offset += lambda_core;
  }

  if ( offset < nlambda ) {
    band_start.push_back( offset );
    band_stop.push_back( nlambda - 1 );
  }

  size_t nband = band_start.size();

  matrix_dist fullRf ( nbins, 1 );
  matrix_dist fullcov ( nbins, 1 );

  for ( size_t band = 0; band < nband; ++band ) {

    if ( ( myp == 0 ) && ( ! quiet ) ) {
      cout << prefix << "  Wavelength band " << band << "/" << nband << " (" << band_start[band] << " - " << band_stop[band] << endl;
    }

    tstart = MPI_Wtime();

    double tsubstart = tstart;
    double tsubstop;

    size_t bandsize = band_stop[ band ] - band_start[ band ] + 1;

    size_t nbins_band = nspec * bandsize;

    matrix_sparse design;

    epsf->projection ( band_start[ band ], band_stop[ band ], design );

    tsubstop = MPI_Wtime();

    if ( ( myp == 0 ) && ( ! quiet ) ) {
      cout << prefix << "    computing A^T = " << tsubstop-tsubstart << " seconds" << endl;
    }

    tsubstart = MPI_Wtime();

    elem::Grid grid ( elem::mpi::COMM_WORLD );
    
    matrix_dist inv ( nbins_band, nbins_band, grid );

    inverse_covariance ( design, invnoise, inv );

    tsubstop = MPI_Wtime();

    if ( ( myp == 0 ) && ( ! quiet ) ) {
      cout << prefix << "    building inverse covariance = " << tsubstop-tsubstart << " seconds" << endl;
    }

    tsubstart = MPI_Wtime();

    matrix_dist W ( nbins_band, nbins_band, grid );
    
    matrix_dist D ( nbins_band, 1, grid );

    eigen_decompose ( inv, D, W );

    tsubstop = MPI_Wtime();

    if ( ( myp == 0 ) && ( ! quiet ) ) {
      cout << prefix << "    eigendecompose inverse covariance = " << tsubstop-tsubstart << " seconds" << endl;
    }

    tsubstart = MPI_Wtime();

    matrix_dist S ( nbins_band, 1, grid );

    norm ( D, W, S );

    tsubstop = MPI_Wtime();

    if ( ( myp == 0 ) && ( ! quiet ) ) {
      cout << prefix << "    compute matrix column norm = " << tsubstop-tsubstart << " seconds" << endl;
    }

    tsubstart = MPI_Wtime();

    matrix_dist z ( nbins_band, 1 );
  
    matrix_dist Rf ( nbins_band, 1 );

    noise_weighted_spec ( design, invnoise, measured, z );

    tsubstop = MPI_Wtime();

    if ( ( myp == 0 ) && ( ! quiet ) ) {
      cout << prefix << "    compute noise weighted spec = " << tsubstop-tsubstart << " seconds" << endl;
    }

    tsubstart = MPI_Wtime();
  
    extract ( D, W, S, z, Rf );

    // FIXME: stitch together into full Rf and covariance. check for agreement within some overlap.



  
    tsubstop = MPI_Wtime();

    tstop = tsubstop;

    if ( ( myp == 0 ) && ( ! quiet ) ) {
      cout << prefix << "    extraction = " << tsubstop-tsubstart << " seconds" << endl;
      cout << prefix << "    total band time = " << tstop-tstart << " seconds" << endl;
    }

  }

  cliq::Finalize();

  return 0;
}


