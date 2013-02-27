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

  cout.precision ( 10 );
  cerr.precision ( 10 );

  size_t lambda_width = 3;
  size_t lambda_overlap = 1;

  size_t spec_first = 0;
  size_t spec_last = 1000000000;
  
  string jsonconf = "";

  string prefix = "harp:  ";

  bool quiet = false;
  bool debug = false;
  bool dosky = false;
  
  fitsfile * fp;
  
  // Declare options
  
  popts::options_description desc ( "Allowed Options" );
  
  desc.add_options()
  ( "help,h", "display usage information" )
  ( "quiet,q", "supress information printing" )
  ( "debug,d", "write out intermediate data products for debugging" )
  ( "skysub", "simultaneously remove common sky spectrum" )
  ( "spec_first", popts::value < size_t > ( &spec_first ), "first spectrum to extract (default 0)" )
  ( "spec_last", popts::value < size_t > ( &spec_last ), "last spectrum to extract (default all spectra)" )
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

  if ( vm.count( "debug" ) ) {
    debug = true;
  }

  if ( vm.count( "skysub" ) ) {
    dosky = true;
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

  double global_start = tstart;
  
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

  vector < bool > is_sky;

  img->read ( measured );
  img->read_noise ( invnoise );
  is_sky = img->sky ();

  tstop = MPI_Wtime();

  if ( debug ) {
    if ( myp == 0 ) {
      string outimg = "specstract_input_image.fits";
      fits::create ( fp, outimg );
      fits::img_append ( fp, imgrows, imgcols );
      fits::write_key ( fp, "EXTNAME", "Data", "harp_specstract input image" );
      fits::img_write ( fp, measured );
      fits::img_append ( fp, imgrows, imgcols );
      fits::write_key ( fp, "EXTNAME", "InvCov", "harp_specstract inverse image covariance" );
      fits::img_write ( fp, invnoise );
      fits::close ( fp );
    }
  }

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

  size_t psf_nspec = epsf->nspec();
  size_t nlambda = epsf->nlambda();

  if ( spec_first > psf_nspec - 1 ) {
    spec_first = psf_nspec - 1;
  }

  if ( spec_last > psf_nspec - 1 ) {
    spec_last = psf_nspec - 1;
  }

  size_t nspec = spec_last - spec_first + 1;

  vector < double > lambda = epsf->lambda();

  size_t nbins = nspec * nlambda;

  tstop = MPI_Wtime();

  if ( ( myp == 0 ) && ( ! quiet ) ) {
    cout << prefix << "  time = " << tstop-tstart << " seconds" << endl;
    cout << prefix << "  image dimensions = " << psf_imgrows << " x " << psf_imgcols << endl;
    cout << prefix << "  " << nspec << " spectra with " << nlambda << " wavelength points each:" << endl;
    cout << prefix << "  lambda = " << lambda[0] << " ... " << lambda[nlambda - 1] << endl;
    cout << prefix << "Extracting spectra " << spec_first << " - " << spec_last << endl;
    cout << prefix << "Processing " << lambda_width << " wavelength points with overlap of " << lambda_overlap << " :" << endl;
  }

  vector < size_t > band_start;
  vector < size_t > band_stop;

  size_t offset = 0;

  while ( offset + lambda_width < nlambda ) {
    band_start.push_back( offset );
    band_stop.push_back( offset + lambda_width - 1 );
    offset += lambda_core;
  }

  if ( offset < nlambda ) {
    band_start.push_back( offset );
    band_stop.push_back( nlambda - 1 );
  }

  size_t nband = band_start.size();

  matrix_dist fullRf ( nbins, 1 );
  matrix_dist fullcov ( nbins, 1 );
  matrix_dist fulltruth ( nbins, 1 );
  matrix_dist fullRtruth ( nbins, 1 );
  vector < bool > fulltruth_sky;

  bool dotruth = false;

  if ( conf.count ( "truth" ) > 0 ) {

    boost::property_tree::ptree truth_props = conf.get_child ( "truth" );

    spec_p truth_spec ( spec::create ( truth_props ) );
    size_t truth_nspec = truth_spec->nspec();
    size_t truth_nlambda = truth_spec->nlambda();

    size_t truth_nbins = truth_nspec * truth_nlambda;

    if ( truth_nbins != psf_nspec * nlambda ) {
      ostringstream o;
      o << "truth spectrum has " << truth_nbins << " spectral bins, but PSF has " << psf_nspec * nlambda << " bins";
      cerr << o.str() << endl;
      MPI_Abort ( MPI_COMM_WORLD, 1 );
    }

    dist_matrix_zero ( fulltruth );
    vector < double > truth_lambda;

    truth_spec->read ( fulltruth, truth_lambda, fulltruth_sky );

    for ( size_t i = 0; i < nlambda; ++i ) {
      if ( fabs ( truth_lambda[i] - lambda[i] ) / lambda[i] > 0.001 ) {
        ostringstream o;
        o << "wavelength point " << i << " does not match between PSF (" << lambda[i] << ") and truth spec (" << truth_lambda[i] << ")";
        HARP_THROW( o.str().c_str() );
      }
    }

    dotruth = true;

  }

  // Process all bands

  //for ( size_t band = 0; band < nband; ++band ) {
  for ( size_t band = 0; band < 1; ++band ) {

    if ( ( myp == 0 ) && ( ! quiet ) ) {
      cout << prefix << "  Wavelength band " << band << "/" << nband << " (" << band_start[band] << " - " << band_stop[band] << ")" << endl;
    }

    tstart = MPI_Wtime();

    double tsubstart = tstart;
    double tsubstop;

    size_t bandsize = band_stop[ band ] - band_start[ band ] + 1;

    size_t nbins_band = nspec * bandsize;

    matrix_sparse design;

    epsf->projection ( spec_first, spec_last, band_start[ band ], band_stop[ band ], design );

    matrix_sparse design_sky;

    if ( dosky ) {
      sky_design ( design, is_sky, design_sky );
    }

    tsubstop = MPI_Wtime();

    if ( ( myp == 0 ) && ( ! quiet ) ) {
      cout << prefix << "    computing A^T = " << tsubstop-tsubstart << " seconds" << endl;
    }

    tsubstart = MPI_Wtime();

    elem::Grid grid ( elem::mpi::COMM_WORLD );
    
    matrix_dist inv ( nbins_band, nbins_band, grid );

    if ( dosky ) {
      inverse_covariance ( design_sky, invnoise, inv );
    } else {
      inverse_covariance ( design, invnoise, inv );
    }

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

    if ( dotruth ) {

      // since we need the explicit resolution matrix anyway, compute it
      // here along with the normalization vector

      matrix_dist Rtruth ( nbins_band, 1 );

      matrix_dist truth_band ( nbins_band, 1 );
      dist_matrix_zero ( truth_band );

      sub_block ( fulltruth, band_start[ band ], 0, nbins_band, 1, truth_band );

      dist_matrix_zero ( Rtruth );

      matrix_dist R;

      resolution ( D, W, S, R );

      elem::Gemv ( elem::NORMAL, 1.0, R, truth_band, 0.0, Rtruth );

      Rtruth.Write ( "Rtruth.txt" );

      // accumulate to global resolution convolved truth



      // FIXME: do this at the end...

      if ( debug ) {

        matrix_local truth_image ( npix, 1 );
        local_matrix_zero ( truth_image );

        spec_project ( design, Rtruth, truth_image );

        if ( myp == 0 ) {
          string outimg = "specstract_Rtruth_projection.fits";
          fits::create ( fp, outimg );
          fits::img_append ( fp, imgrows, imgcols );
          fits::write_key ( fp, "EXTNAME", "RTruth", "harp_specstract resolution convolved truth" );
          fits::img_write ( fp, truth_image );
          fits::close ( fp );
        }

      }

    } else {

      // we just need the norm

      norm ( D, W, S );

    }    

    tsubstop = MPI_Wtime();

    if ( ( myp == 0 ) && ( ! quiet ) ) {
      cout << prefix << "    compute matrix column norm = " << tsubstop-tsubstart << " seconds" << endl;
    }

    tsubstart = MPI_Wtime();

    matrix_dist z ( nbins_band, 1 );
  
    matrix_dist Rf ( nbins_band, 1 );

    matrix_dist Rf_err ( nbins_band, 1 );    

    noise_weighted_spec ( design, invnoise, measured, z );

    tsubstop = MPI_Wtime();

    if ( ( myp == 0 ) && ( ! quiet ) ) {
      cout << prefix << "    compute noise weighted spec = " << tsubstop-tsubstart << " seconds" << endl;
    }

    tsubstart = MPI_Wtime();
  
    extract ( D, W, S, z, Rf, Rf_err );



    // FIXME: stitch together into full Rf and covariance. check for agreement within some overlap.

    Rf.Write( "Rf.txt" );
    Rf_err.Write( "Rf_err.txt" );


    // FIXME: put this at the end of the program and dump the full projected image

    matrix_local solution_image ( npix, 1 );
    local_matrix_zero ( solution_image );

    spec_project ( design, Rf, solution_image );

    if ( myp == 0 ) {
      string outimg = "specstract_Rf_projection.fits";
      fits::create ( fp, outimg );
      fits::img_append ( fp, imgrows, imgcols );
      fits::write_key ( fp, "EXTNAME", "Rf", "harp_specstract solution" );
      fits::img_write ( fp, solution_image );
      fits::close ( fp );
    }


  
    tsubstop = MPI_Wtime();

    tstop = tsubstop;

    if ( ( myp == 0 ) && ( ! quiet ) ) {
      cout << prefix << "    extraction = " << tsubstop-tsubstart << " seconds" << endl;
      cout << prefix << "    total band time = " << tstop-tstart << " seconds" << endl;
    }

  }

  double global_stop = MPI_Wtime();

  if ( ( myp == 0 ) && ( ! quiet ) ) {
    cout << prefix << "Total run time = " << global_stop-global_start << " seconds" << endl;
  }

  cliq::Finalize();

  return 0;
}


