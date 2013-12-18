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

  size_t spec_width = 0;
  
  string jsonconf = "";

  string prefix = "harp:  ";

  bool quiet = false;
  bool debug = false;
  bool dosky = false;
  bool doskysub = false;
  int planproc = 0;

  int gangsize = np;
  
  fitsfile * fp;
  
  // Declare options
  
  popts::options_description desc ( "Allowed Options" );
  
  desc.add_options()
  ( "help,h", "display usage information" )
  ( "quiet,q", "supress information printing" )
  ( "debug,d", "write out intermediate data products for debugging" )
  ( "sky", "simultaneously remove common sky spectrum" )
  ( "skysub", "solve for common sky and subtract" )
  ( "gangsize", popts::value < int > ( &gangsize ), "number of processes per gang (choose perfect square number if possible)" )
  ( "spec_width", popts::value < size_t > ( &spec_width ), "number of spectra to process at once (default = all)" )
  ( "lambda_width", popts::value < size_t > ( &lambda_width ), "maximum wavelength points to process simultaneously" )
  ( "lambda_overlap", popts::value < size_t > ( &lambda_overlap ), "minimum wavelength points to overlap" )
  ( "plan", popts::value < int > ( &planproc ), "# procs to simulate. Check inputs and compute distribution." )
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
  
  if ( vm.count( "quiet" ) ) {
    quiet = true;
  }

  if ( vm.count( "debug" ) ) {
    debug = true;
  }

  if ( vm.count( "skysub" ) ) {
    doskysub = true;
  }

  if ( vm.count( "sky" ) ) {
    dosky = true;
  }

  bool skystuff = dosky || doskysub;

  // create global process grid

  elem::Grid grid ( MPI_COMM_WORLD );

  // Read metadata
  
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

  if ( ( myp == 0 ) && ( ! quiet ) ) {
    cout << prefix << "  time = " << tstop-tstart << " seconds" << endl;
    cout << prefix << "  dimensions = " << imgrows << " x " << imgcols << endl;
  }

  if ( debug && ( planproc < 1 ) ) {

    tstart = MPI_Wtime();

    if ( myp == 0 ) {
      cout << prefix << "(debug mode) Dumping input image to disk..." << endl;

      string outimg = "harp_image_input.fits";
      fits::create ( fp, outimg );
      fits::img_append ( fp, imgrows, imgcols );
      fits::write_key ( fp, "EXTNAME", "Data", "input image" );
      fits::img_write ( fp, measured );
      fits::img_append ( fp, imgrows, imgcols );
      fits::write_key ( fp, "EXTNAME", "InvCov", "inverse image covariance" );
      fits::img_write ( fp, invnoise );
      fits::close ( fp );
    }

    tstop = MPI_Wtime();

    if ( ( myp == 0 ) && ( ! quiet ) ) {
      cout << prefix << "  time = " << tstop-tstart << " seconds" << endl;
    }

  }

  if ( ( myp == 0 ) && ( ! quiet ) ) {
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
  vector < double > lambda = epsf->lambda();
  size_t psf_nbins = psf_nspec * nlambda;

  if ( spec_width == 0 ) {
    spec_width = psf_nspec;
  }

  // If we are doing sky subtraction, force user to solve for bands of all spectra

  if ( skystuff && ( spec_width != psf_nspec ) ) {
    cerr << prefix << "when doing sky subtraction, spec_width must be all spectra" << endl;
    ret = MPI_Abort ( MPI_COMM_WORLD, 1 );
  }

  tstop = MPI_Wtime();

  if ( ( myp == 0 ) && ( ! quiet ) ) {
    cout << prefix << "  time = " << tstop-tstart << " seconds" << endl;
    cout << prefix << "  image dimensions = " << psf_imgrows << " x " << psf_imgcols << endl;
    cout << prefix << "  " << psf_nspec << " spectra with " << nlambda << " wavelength points each:" << endl;
    cout << prefix << "  lambda = " << lambda[0] << " ... " << lambda[nlambda - 1] << endl;
  }

  // Determine reduced problem size for sky subtraction

  size_t psf_nspec_obj = 0;

  for ( size_t i = 0; i < psf_nspec; ++i ) {
    if ( ! is_sky[i] ) {
      ++psf_nspec_obj;
    }
  }

  //size_t psf_nspec_sky = psf_nspec + 1;
  //size_t psf_nbins_sky = psf_nspec_sky * nlambda;

  size_t psf_nspec_sky = psf_nspec_obj + 1;
  size_t psf_nbins_sky = psf_nspec_sky * nlambda;

  vector < bool > is_sky_solve ( psf_nspec_sky );
  vector < int > sky_lookup ( psf_nspec );
  int sky_off = 0;

  for ( size_t i = 0; i < psf_nspec; ++i ) {
    if ( is_sky[i] ) {
      sky_lookup[i] = psf_nspec_sky - 1;
    } else {
      sky_lookup[i] = sky_off;
      is_sky_solve[ sky_off ] = false;
      ++sky_off;
    }
  }
  is_sky_solve[ sky_off ] = true;

  size_t psf_nspec_solve;
  size_t psf_nbins_solve;

  if ( skystuff ) {
    psf_nspec_solve = psf_nspec_sky;
    psf_nbins_solve = psf_nbins_sky;
  } else {
    psf_nspec_solve = psf_nspec;
    psf_nbins_solve = psf_nbins;
  }

  // Determine spectral chunks

  vector < size_t > spec_start;
  vector < size_t > spec_stop;

  size_t nspec_chunk = (size_t) ( psf_nspec / spec_width );

  for ( size_t i = 0; i < nspec_chunk; ++i ) {
    spec_start.push_back ( i * spec_width );
    spec_stop.push_back ( (i + 1) * spec_width - 1 );
  }

  if ( psf_nspec > nspec_chunk * spec_width ) {
    spec_start.push_back ( nspec_chunk * spec_width );
    spec_stop.push_back ( psf_nspec - 1 );
    ++nspec_chunk;
  } 

  // Determine lambda chunks

  size_t lambda_core = lambda_width - 2 * lambda_overlap;

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

  // select output block from solved spectra

  vector < size_t > band_out;
  vector < size_t > band_write;

  if ( nband == 1 ) {

    // there is only one band...
    band_out.push_back( 0 );
    band_write.push_back( band_stop[0] - band_start[0] + 1 );

  } else {

    band_out.push_back( 0 );
    band_write.push_back( lambda_overlap + lambda_core );

    for ( size_t i = 1; i < ( nband - 1 ); ++i ) {
      band_out.push_back( band_start[i] + lambda_overlap );
      band_write.push_back( lambda_core );
    }

    band_out.push_back( band_start[ nband - 1 ] + lambda_overlap );
    band_write.push_back( band_stop[ nband - 1 ] + 1 - ( band_start[ nband - 1 ] + lambda_overlap ) );

  }

  // Compute gang distribution

  int ngang;
  if ( planproc > 0 ) {
    ngang = (int)( planproc / gangsize );
  } else {
    ngang = (int)( np / gangsize );
  }
  if ( ( myp == 0 ) && ( ngang < 1 ) ) {
    cerr << "Process gang size is larger than the total number of processes!" << endl;
    cliq::Finalize();
    return 0;
  }

  int gangtot = ngang * gangsize;
  if ( ( myp == 0 ) && ( ! quiet ) ) {
    cout << prefix << "Using " << ngang << " gangs of " << gangsize << " processes each" << endl;
  }

  if ( planproc > 0 ) {
    if ( gangtot < planproc ) {
      if ( ( myp == 0 ) && ( ! quiet ) ) {
        cout << prefix << "WARNING: " << (planproc-gangtot) << " processes are idle" << endl;
      }
    }
  } else {
    if ( gangtot < np ) {
      if ( ( myp == 0 ) && ( ! quiet ) ) {
        cout << prefix << "WARNING: " << (np-gangtot) << " processes are idle" << endl;
      }
    }
  }

  vector < int > fake_gang ( planproc );
  vector < int > fake_grank ( planproc );

  int gang = MPI_UNDEFINED;
  int grank = MPI_UNDEFINED;

  MPI_Comm gcomm;
  MPI_Comm grootcomm;

  if ( planproc > 0 ) {
    fake_gang.resize ( planproc );
    fake_grank.resize ( planproc );
    for ( int fakep = 0; fakep < planproc; ++fakep ) {
      fake_gang[ fakep ] = (int)( fakep / gangsize );
      fake_grank[ fakep ] = fakep % gangsize;
    }
  } else {
    gang = (int)( myp / gangsize );
    grank = myp % gangsize;

    ret = MPI_Comm_split ( MPI_COMM_WORLD, gang, grank, &gcomm );
    mpi_check ( MPI_COMM_WORLD, ret );

    if ( grank == 0 ) {
      ret = MPI_Comm_split ( MPI_COMM_WORLD, 0, gang, &grootcomm );
    } else {
      ret = MPI_Comm_split ( MPI_COMM_WORLD, MPI_UNDEFINED, 0, &grootcomm );
    }
    mpi_check ( MPI_COMM_WORLD, ret );
  }

  // distribute spectral chunks among the gangs

  size_t total_chunks = nspec_chunk * nband;

  if ( ( myp == 0 ) && ( total_chunks < ngang ) ) {
    cerr << "There are more process gangs than chunks: " << endl;
    cerr << "   use fewer total processes, larger gangs, or smaller chunks!" << endl;
    cliq::Finalize();
    return 0;
  }

  size_t gang_nchunk = (size_t) ( total_chunks / ngang );
  size_t gang_offset;

  size_t leftover = total_chunks % ngang;

  vector < size_t > fake_gang_nchunk ( planproc );
  vector < size_t > fake_gang_offset ( planproc );

  if ( planproc > 0 ) {

    for ( int fakep = 0; fakep < planproc; ++fakep ) {
      fake_gang_nchunk[ fakep ] = gang_nchunk;
      if ( fake_gang[ fakep ] < leftover ) {
        ++fake_gang_nchunk[ fakep ];
        fake_gang_offset[ fakep ] = fake_gang[ fakep ] * fake_gang_nchunk[ fakep ];
      } else {
        fake_gang_offset[ fakep ] = ( (fake_gang_nchunk[ fakep ] + 1) * leftover ) + ( fake_gang_nchunk[ fakep ] * (fake_gang[ fakep ] - leftover) );
      }

    }

  } else {

    if ( gang < leftover ) {
      ++gang_nchunk;
      gang_offset = gang * gang_nchunk;
    } else {
      gang_offset = ( (gang_nchunk + 1) * leftover ) + ( gang_nchunk * (gang - leftover) );
    }

  }

  // Print out "planning information"

  if ( ( myp == 0 ) && ( ! quiet ) ) {
    cout << prefix << "Extracting " << total_chunks << " spectral chunks, each with " << spec_width << " spectra and " << lambda_width << " lambda points ( overlap = " << lambda_overlap << " )" << endl;
  }

  if ( planproc > 0 ) {

    if ( myp == 0 ) {
      for ( int g = 0; g < ngang; ++g ) {
        for ( int fakep = 0; fakep < planproc; ++fakep ) {
          if ( ( fake_gang[ fakep ] == g ) && ( fake_grank[ fakep ] == 0 ) ) {
            cout << prefix << "  gang " << g << ": assigned spectral chunks " << fake_gang_offset[ fakep ] << " - " << (fake_gang_offset[ fakep ] + fake_gang_nchunk[ fakep ] - 1) << endl;
          }
        }
      }
    }

    MPI_Barrier ( MPI_COMM_WORLD );

    cliq::Finalize();
    return 0;

  } else {

    for ( int g = 0; g < ngang; ++g ) {
      if ( ( g == gang ) && ( grank == 0 ) ) {
        cout << prefix << "  gang " << g << ": assigned spectral chunks " << gang_offset << " - " << (gang_offset + gang_nchunk - 1) << endl;
      }
      MPI_Barrier ( MPI_COMM_WORLD );
    }
    
  }

  // create gang process grid

  elem::Grid gang_grid ( gcomm );

  // global distributed spectral products

  matrix_dist fullf ( psf_nbins_solve, 1, grid );
  dist_matrix_zero ( fullf );

  matrix_dist fullRf ( psf_nbins_solve, 1, grid );
  dist_matrix_zero ( fullRf );

  matrix_dist fullerr ( psf_nbins_solve, 1, grid );
  dist_matrix_zero ( fullerr );

  matrix_dist fulltruth ( psf_nbins_solve, 1, grid );
  dist_matrix_zero ( fulltruth );

  matrix_dist fullRtruth ( psf_nbins_solve, 1, grid );
  dist_matrix_zero ( fullRtruth );

  // gang distributed spectral products

  matrix_dist gang_fullf ( psf_nbins_solve, 1, gang_grid );
  dist_matrix_zero ( gang_fullf );

  matrix_dist gang_fullRf ( psf_nbins_solve, 1, gang_grid );
  dist_matrix_zero ( gang_fullRf );

  matrix_dist gang_fullerr ( psf_nbins_solve, 1, gang_grid );
  dist_matrix_zero ( gang_fullerr );

  matrix_dist gang_fulltruth ( psf_nbins_solve, 1, gang_grid );
  dist_matrix_zero ( gang_fulltruth );

  matrix_dist gang_fullRtruth ( psf_nbins_solve, 1, gang_grid );
  dist_matrix_zero ( gang_fullRtruth );

  vector < bool > fulltruth_sky;

  bool dotruth = false;

  if ( conf.count ( "truth" ) > 0 ) {

    if ( ( myp == 0 ) && ( ! quiet ) ) {
      cout << prefix << "Reading input truth spectra..." << endl;
    }

    tstart = MPI_Wtime();

    boost::property_tree::ptree truth_props = conf.get_child ( "truth" );

    spec_p truth_spec ( spec::create ( truth_props ) );
    size_t truth_nspec = truth_spec->nspec();
    size_t truth_nlambda = truth_spec->nlambda();

    size_t truth_nbins = truth_nspec * truth_nlambda;

    if ( truth_nbins != psf_nspec_solve * nlambda ) {
      ostringstream o;
      o << "truth spectrum has " << truth_nbins << " spectral bins, but we are solving for " << psf_nspec_solve * nlambda << " bins";
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

    if ( skystuff ) {
      for ( size_t i = 0; i < psf_nspec_sky - 1; ++i ) {
        if ( fulltruth_sky[i] ) {
          ostringstream o;
          o << "doing sky subtraction, but truth spec " << i << " is a sky fiber!";
          HARP_THROW( o.str().c_str() );
        }
      }
      if ( ! fulltruth_sky[ psf_nspec_solve - 1 ] ) {
        ostringstream o;
        o << "doing sky subtraction, but last truth spec is not a sky fiber!";
        HARP_THROW( o.str().c_str() );
      }
    } else {
      for ( size_t i = 0; i < psf_nspec; ++i ) {
        if ( fulltruth_sky[i] != is_sky[i] ) {
          ostringstream o;
          o << "truth sky point " << i << " (" << fulltruth_sky[i] << ") does not match input image (" << is_sky[i] << ")";
          HARP_THROW( o.str().c_str() );
        }
      }
    }

    dotruth = true;

    gang_distribute ( fulltruth, gang_fulltruth );

    tstop = MPI_Wtime();

    if ( ( myp == 0 ) && ( ! quiet ) ) {
      cout << prefix << "  time = " << tstop-tstart << " seconds" << endl;
    }

  }

  // aggregate timing

  double tot_design = 0.0;
  double time_design;
  double tot_inverse = 0.0;
  double time_inverse;
  double tot_eigen = 0.0;
  double time_eigen;
  double tot_norm = 0.0;
  double time_norm;
  double tot_nsespec = 0.0;
  double time_nsespec;
  double tot_extract = 0.0;
  double time_extract;
  double time_chunk;

  // root gang processes implement RMA lock for gangs to take turns writing to stdout

  int * root_lockbuf;

  MPI_Win printlock;

  if ( grank == 0 ) {
    if ( gang == 0 ) {
      ret = MPI_Alloc_mem ( sizeof(int), MPI_INFO_NULL, (void*)&root_lockbuf );
      mpi_check ( grootcomm, ret );

      ret = MPI_Win_create ( (void*)(root_lockbuf), sizeof(int), sizeof(int), MPI_INFO_NULL, grootcomm, &printlock );
    } else {
      ret = MPI_Win_create ( NULL, 0, sizeof(int), MPI_INFO_NULL, grootcomm, &printlock );
    }
    mpi_check ( grootcomm, ret );
  }

  if ( ( myp == 0 ) && ( ! quiet ) ) {
    cout << prefix << "Processing..." << endl;
  }

  // Process all chunks for this gang

  for ( size_t gchunk = 0; gchunk < gang_nchunk; ++gchunk ) {

    size_t band = (size_t) ( (gang_offset + gchunk) / nspec_chunk );

    size_t spec = (gang_offset + gchunk) % nspec_chunk;
    
    size_t bandsize = band_stop[ band ] - band_start[ band ] + 1;

    tstart = MPI_Wtime();

    double tsubstart = tstart;
    double tsubstop;

    size_t nspec = spec_stop[ spec ] - spec_start[ spec ] + 1;

    size_t spec_start_solve;
    size_t spec_stop_solve;

    if ( skystuff ) {
      // we are doing the full range of spectra (enforced by checks on spec_width)
      spec_start_solve = 0;
      spec_stop_solve = psf_nspec_sky - 1;
    } else {
      spec_start_solve = spec_start[ spec ];
      spec_stop_solve = spec_stop[ spec ];
    }

    size_t nspec_solve = spec_stop_solve - spec_start_solve + 1;
    
    size_t nbins = nspec * bandsize;
    size_t nbins_solve = nspec_solve * bandsize;

    matrix_sparse design ( gcomm );

    epsf->projection ( spec_start[ spec ], spec_stop[ spec ], band_start[ band ], band_stop[ band ], design );

    matrix_sparse design_sky ( gcomm );

    if ( skystuff ) {
      //sky_design2 ( design, is_sky, design_sky );
      sky_design ( design, is_sky, design_sky, dosky );
    }

    tsubstop = MPI_Wtime();
    time_design = ( tsubstop - tsubstart );

    tsubstart = MPI_Wtime();
    
    matrix_dist inv ( nbins_solve, nbins_solve, gang_grid );

    if ( skystuff ) {
      inverse_covariance ( design_sky, invnoise, inv );
    } else {
      inverse_covariance ( design, invnoise, inv );
    }

    tsubstop = MPI_Wtime();
    time_inverse = ( tsubstop - tsubstart );

    tsubstart = MPI_Wtime();

    matrix_dist W ( nbins_solve, nbins_solve, gang_grid );
    
    matrix_dist D ( nbins_solve, 1, gang_grid );

    eigen_decompose ( inv, D, W );

    elem::Write ( D, "eigen", "dbg_eigevals.txt" );

    tsubstop = MPI_Wtime();
    time_eigen = ( tsubstop - tsubstart );

    matrix_dist S ( nbins_solve, 1, gang_grid );

    matrix_dist out_spec ( nspec_solve * band_write [ band ], 1, gang_grid );

    if ( dotruth ) {

      tsubstart = MPI_Wtime();

      // since we need the explicit resolution matrix anyway, compute it
      // here along with the normalization vector

      matrix_dist Rtruth ( nbins_solve, 1, gang_grid );
      dist_matrix_zero ( Rtruth );

      matrix_dist truth_band ( nbins_solve, 1, gang_grid );
      dist_matrix_zero ( truth_band );

      sub_spec ( gang_fulltruth, psf_nspec_solve, spec_start_solve, nspec_solve, band_start[ band ], bandsize, truth_band );

      matrix_dist R ( W );

      resolution ( D, W, S, R );

      elem::Gemv ( elem::NORMAL, 1.0, R, truth_band, 0.0, Rtruth );
      
      // accumulate to global resolution convolved truth

      sub_spec ( Rtruth, nspec_solve, 0, nspec_solve, band_out[ band ] - band_start[ band ], band_write[ band ], out_spec );

      accum_spec ( gang_fullRtruth, psf_nspec_solve, spec_start_solve, nspec_solve, band_out[ band ], band_write[ band ], out_spec );

      tsubstop = MPI_Wtime();
      time_norm = ( tsubstop - tsubstart );

    } else {

      // we just need the norm

      tsubstart = MPI_Wtime();

      norm ( D, W, S );

      tsubstop = MPI_Wtime();
      time_norm = ( tsubstop - tsubstart );

    }

    matrix_dist z ( nbins_solve, 1, gang_grid );
    dist_matrix_zero ( z );

    matrix_dist Rf ( nbins_solve, 1, gang_grid );
    dist_matrix_zero ( Rf );

    matrix_dist f ( nbins_solve, 1, gang_grid );
    dist_matrix_zero ( f );

    tsubstart = MPI_Wtime();

    if ( skystuff ) {
      noise_weighted_spec ( design_sky, invnoise, measured, z );
    } else {
      noise_weighted_spec ( design, invnoise, measured, z );
    }

    tsubstop = MPI_Wtime();
    time_nsespec = ( tsubstop - tsubstart );

    tsubstart = MPI_Wtime();
  
    extract ( D, W, S, z, Rf, f );

    // accumulate to global solution

    sub_spec ( f, nspec_solve, 0, nspec_solve, band_out[ band ] - band_start[ band ], band_write[ band ], out_spec );
    accum_spec ( gang_fullf, psf_nspec_solve, spec_start_solve, nspec_solve, band_out[ band ], band_write[ band ], out_spec );

    MPI_Barrier ( gcomm );

    sub_spec ( Rf, nspec_solve, 0, nspec_solve, band_out[ band ] - band_start[ band ], band_write[ band ], out_spec );
    accum_spec ( gang_fullRf, psf_nspec_solve, spec_start_solve, nspec_solve, band_out[ band ], band_write[ band ], out_spec );

    MPI_Barrier ( gcomm );

    sub_spec ( S, nspec_solve, 0, nspec_solve, band_out[ band ] - band_start[ band ], band_write[ band ], out_spec );
    accum_spec ( gang_fullerr, psf_nspec_solve, spec_start_solve, nspec_solve, band_out[ band ], band_write[ band ], out_spec );

    tsubstop = MPI_Wtime();
    time_extract = ( tsubstop - tsubstart );

    tstop = tsubstop;

    time_chunk = tstop - tstart;

    if ( ( grank == 0 ) && ( ! quiet ) ) {

      ret = MPI_Win_lock ( MPI_LOCK_EXCLUSIVE, 0, 0, printlock );
      mpi_check ( grootcomm, ret );

      cout << prefix << "  gang " << gang << ": finished chunk " << (gang_offset + gchunk) << endl;
      cout << prefix << "    computing A^T = " << time_design << " seconds" << endl;
      cout << prefix << "    building inverse covariance = " << time_inverse << " seconds" << endl;
      cout << prefix << "    eigendecompose inverse covariance = " << time_eigen << " seconds" << endl;
      if ( dotruth ) {
        cout << prefix << "    compute column norm and resolution convolved truth = " << time_norm << " seconds" << endl;
      } else {
        cout << prefix << "    compute column norm = " << time_norm << " seconds" << endl;
      }
      cout << prefix << "    compute noise weighted spec = " << time_nsespec << " seconds" << endl;
      cout << prefix << "    extraction = " << time_extract << " seconds" << endl;
      cout << prefix << "    total band time = " << time_chunk << " seconds" << endl;

      // free the lock
      ret = MPI_Win_unlock ( 0, printlock );
      mpi_check ( grootcomm, ret );

      tot_design += time_design;
      tot_inverse += time_inverse;
      tot_eigen += time_eigen;
      tot_norm += time_norm;
      tot_nsespec += time_nsespec;
      tot_extract += time_extract;

    }

    MPI_Barrier ( gcomm );

  }

  if ( grank == 0 ) {
    MPI_Win_free ( &printlock );
    if ( gang == 0 ) {
      ret = MPI_Free_mem ( (void*)root_lockbuf );
      mpi_check ( grootcomm, ret );
    }
  }
  
  // Merge gang-wise outputs

  gang_accum ( gang_fullf, fullf );
  gang_accum ( gang_fullRf, fullRf );
  gang_accum ( gang_fullerr, fullerr );
  gang_accum ( gang_fullRtruth, fullRtruth );

  // subtract sky if needed

  matrix_dist fullRfsky ( psf_nbins_solve, 1, grid );
  dist_matrix_zero ( fullRfsky );

  matrix_dist fullerrsky ( psf_nbins_solve, 1, grid );
  dist_matrix_zero ( fullerrsky );

  if ( doskysub ) {
    sky_subtract ( psf_nspec_obj, fullRf, fullerr, fullRfsky, fullerrsky );
  }

  // Total timings

  if ( ( myp == 0 ) && ( ! quiet ) ) {
    cout << prefix << "Aggregate Timings:" << endl;
    cout << prefix << "  Build design matrix = " << tot_design << " seconds" << endl;
    cout << prefix << "  Build inverse covariance = " << tot_inverse << " seconds" << endl;
    cout << prefix << "  Eigendecompose inverse = " << tot_eigen << " seconds" << endl;
    cout << prefix << "  Compute column norm = " << tot_norm << " seconds" << endl;
    cout << prefix << "  Compute noise weighted spec = " << tot_nsespec << " seconds" << endl;
    cout << prefix << "  Extract spectra = " << tot_extract << " seconds" << endl;
  }

  // Write outputs

  if ( ( myp == 0 ) && ( ! quiet ) ) {
    cout << prefix << "Writing solution and error..." << endl;
  }

  tstart = MPI_Wtime();

  boost::property_tree::ptree solution_spec_props;
  solution_spec_props.put ( "format", "specter" );
  solution_spec_props.put ( "nspec", psf_nspec_solve );
  solution_spec_props.put ( "nlambda", nlambda );

  spec_p solution_spec ( spec::create ( solution_spec_props ) );

  solution_spec->write ( "harp_spec_Rf.fits", fullRf, lambda, is_sky_solve );
  solution_spec->write ( "harp_spec_Rf-err.fits", fullerr, lambda, is_sky_solve );

  if ( doskysub ) {
    solution_spec->write ( "harp_spec_Rfsky.fits", fullRfsky, lambda, is_sky_solve );
    solution_spec->write ( "harp_spec_Rfsky-err.fits", fullerrsky, lambda, is_sky_solve );
  }

  if ( debug ) {
    elem::Write ( fullf, "spec_f", "harp_spec_f.txt" );
    elem::Write ( fullRf, "spec_Rf", "harp_spec_Rf.txt" );
    elem::Write ( fullerr, "spec_Rf-err", "harp_spec_Rf-err.txt" );
    if ( doskysub ) {
      elem::Write ( fullRfsky, "spec_Rfsky", "harp_spec_Rfsky.txt" );
      elem::Write ( fullerrsky, "spec_Rf-errsky", "harp_spec_Rf-errsky.txt" );
    }
  }

  tstop = MPI_Wtime();

  if ( ( myp == 0 ) && ( ! quiet ) ) {
    cout << prefix << "  time = " << tstop-tstart << " seconds" << endl;
  }

  matrix_local solution_image;

  if ( debug ) {

    if ( ( myp == 0 ) && ( ! quiet ) ) {
      cout << prefix << "(debug mode) Writing projected, deconvolved spectra..." << endl;
    }

    tstart = MPI_Wtime();

    solution_image.ResizeTo ( npix, 1 );
    local_matrix_zero ( solution_image );

    matrix_sparse design ( psf_nspec * nlambda, npix, MPI_COMM_WORLD );

    epsf->projection ( 0, psf_nspec - 1, 0, nlambda - 1, design );

    matrix_sparse design_sky ( gcomm );

    if ( skystuff ) {
      //sky_design2 ( design, is_sky, design_sky );
      sky_design ( design, is_sky, design_sky, dosky );
      spec_project ( design_sky, fullf, solution_image );
    } else {
      spec_project ( design, fullf, solution_image );
    }

    if ( myp == 0 ) {
      string outimg = "harp_image_f-project.fits";
      fits::create ( fp, outimg );
      fits::img_append ( fp, imgrows, imgcols );
      fits::write_key ( fp, "EXTNAME", "Solution", "Projected solution" );
      fits::img_write ( fp, solution_image );
      fits::close ( fp );
    }

    tstop = MPI_Wtime();

    if ( ( myp == 0 ) && ( ! quiet ) ) {
      cout << prefix << "  time = " << tstop-tstart << " seconds" << endl;
    }

  }

  if ( dotruth ) {

    if ( ( myp == 0 ) && ( ! quiet ) ) {
      cout << prefix << "Writing resolution convolved truth..." << endl;
    }

    tstart = MPI_Wtime();

    boost::property_tree::ptree rtruth_spec_props;
    rtruth_spec_props.put ( "format", "specter" );
    rtruth_spec_props.put ( "nspec", psf_nspec_solve );
    rtruth_spec_props.put ( "nlambda", nlambda );

    spec_p rtruth_spec ( spec::create ( rtruth_spec_props ) );

    rtruth_spec->write ( "harp_spec_Rtruth.fits", fullRtruth, lambda, is_sky_solve );

    if ( debug ) {
      elem::Write ( fulltruth, "spec_truth", "harp_spec_truth.txt" );
      elem::Write ( fullRtruth, "spec_Rtruth", "harp_spec_Rtruth.txt" );
    }

    tstop = MPI_Wtime();

    if ( ( myp == 0 ) && ( ! quiet ) ) {
      cout << prefix << "  time = " << tstop-tstart << " seconds" << endl;
    }

    if ( ( myp == 0 ) && ( ! quiet ) ) {
      cout << prefix << "Computing spectral Chi-square..." << endl;
    }

    tstart = MPI_Wtime();

    matrix_local loc_fullRtruth;
    matrix_local loc_fullRf;
    matrix_local loc_fullerr;
    
    elem::AxpyInterface < double > globloc;

    globloc.Attach( elem::GLOBAL_TO_LOCAL, fullRtruth );
    if ( myp == 0 ) {
      loc_fullRtruth.ResizeTo ( fullRtruth.Height(), 1 );
      local_matrix_zero ( loc_fullRtruth );
      globloc.Axpy ( 1.0, loc_fullRtruth, 0, 0 );
    }
    globloc.Detach();

    globloc.Attach( elem::GLOBAL_TO_LOCAL, fullRf );
    if ( myp == 0 ) {
      loc_fullRf.ResizeTo ( fullRf.Height(), 1 );
      local_matrix_zero ( loc_fullRf );
      globloc.Axpy ( 1.0, loc_fullRf, 0, 0 );
    }
    globloc.Detach();

    globloc.Attach( elem::GLOBAL_TO_LOCAL, fullerr );
    if ( myp == 0 ) {
      loc_fullerr.ResizeTo ( fullerr.Height(), 1 );
      local_matrix_zero ( loc_fullerr );
      globloc.Axpy ( 1.0, loc_fullerr, 0, 0 );
    }
    globloc.Detach();

    if ( myp == 0 ) {
      fstream outfile ( "harp_spec_chisq.txt", ios::out );
      outfile.precision(15);

      if ( ! outfile.is_open() ) {
        std::ostringstream o;
        o << "cannot open output spec chisq file";
        HARP_THROW( o.str().c_str() );
      }

      double chisq_reduced = 0.0;
      double val;
      for ( size_t i = 0; i < loc_fullRf.Height(); ++i ) {
        val = ( loc_fullRf.Get ( i, 0 ) - loc_fullRtruth.Get ( i, 0 ) ) / loc_fullerr.Get ( i, 0 );
        val *= val;
        chisq_reduced += val;
        outfile << val << endl;
      }

      chisq_reduced /= (double)( loc_fullRf.Height() - 1 );

      cout << prefix << "  Reduced Chi square = " << chisq_reduced << endl;

      outfile.close();
    }

    tstop = MPI_Wtime();

    if ( ( myp == 0 ) && ( ! quiet ) ) {
      cout << prefix << "  time = " << tstop-tstart << " seconds" << endl;
    }

    matrix_local truth_image;

    if ( debug ) {

      if ( ( myp == 0 ) && ( ! quiet ) ) {
        cout << prefix << "(debug mode) Writing projected truth..." << endl;
      }

      tstart = MPI_Wtime();

      matrix_sparse design ( psf_nspec * nlambda, npix, MPI_COMM_WORLD );

      epsf->projection ( 0, psf_nspec - 1, 0, nlambda - 1, design );

      truth_image.ResizeTo ( npix, 1 );
      local_matrix_zero ( truth_image );

      matrix_sparse design_sky ( psf_nspec_solve * nlambda, npix, MPI_COMM_WORLD );

      if ( skystuff ) {
        //sky_design2 ( design, is_sky, design_sky );
        sky_design ( design, is_sky, design_sky, dosky );
        spec_project ( design_sky, fulltruth, truth_image );
      } else {
        spec_project ( design, fulltruth, truth_image );
      }

      if ( myp == 0 ) {
        string outimg = "harp_image_truth-project.fits";
        fits::create ( fp, outimg );
        fits::img_append ( fp, imgrows, imgcols );
        fits::write_key ( fp, "EXTNAME", "Truth", "Projected truth" );
        fits::img_write ( fp, truth_image );
        fits::close ( fp );
      }

      tstop = MPI_Wtime();

      if ( ( myp == 0 ) && ( ! quiet ) ) {
        cout << prefix << "  time = " << tstop-tstart << " seconds" << endl;
      }

      if ( ( myp == 0 ) && ( ! quiet ) ) {
        cout << prefix << "Computing pixel space Chi-square..." << endl;
      }

      tstart = MPI_Wtime();

      if ( myp == 0 ) {
        fstream outfile ( "harp_image_chisq.txt", ios::out );
        outfile.precision(15);

        if ( ! outfile.is_open() ) {
          std::ostringstream o;
          o << "cannot open output spec chisq file";
          HARP_THROW( o.str().c_str() );
        }

        double chisq_reduced = 0.0;
        double val;
        for ( size_t i = 0; i < npix; ++i ) {
          val = ( solution_image.Get ( i, 0 ) - truth_image.Get ( i, 0 ) ) * sqrt ( invnoise.Get ( i, 0 ) );
          val *= val;
          chisq_reduced += val;
          outfile << val << endl;
        }

        // (same number of degrees as freedom, even though projected onto pixels)
        chisq_reduced /= (double)( loc_fullRf.Height() - 1 );

        cout << prefix << "  Reduced Chi square = " << chisq_reduced << endl;

        outfile.close();
      }

      tstop = MPI_Wtime();

      if ( ( myp == 0 ) && ( ! quiet ) ) {
        cout << prefix << "  time = " << tstop-tstart << " seconds" << endl;
      }

    }

  }

  


  double global_stop = MPI_Wtime();

  if ( ( myp == 0 ) && ( ! quiet ) ) {
    cout << prefix << "Total run time = " << global_stop-global_start << " seconds" << endl;
  }

  cliq::Finalize();

  return 0;
}

