/*
  High Performance Astrophysical Reconstruction and Processing (HARP)
  (c) 2014-2015, The Regents of the University of California, 
  through Lawrence Berkeley National Laboratory.  See top
  level LICENSE file for details.
*/

#include <iostream>
#include <cstdio>

extern "C" {
  #include <unistd.h>
}

#include <boost/program_options.hpp>

#include <harp_mpi.hpp>


using namespace std;
using namespace harp;


namespace popts = boost::program_options;


int main ( int argc, char *argv[] ) {

  El::Initialize ( argc, argv );

  boost::mpi::environment env;

  double tstart;
  double tstop;

  cout.precision ( 10 );
  cerr.precision ( 10 );

  size_t lambda_width = 50;
  double lambda_res = 0.6;
  double lambda_b_min = 3579.0;
  double lambda_b_max = 5939.0;
  double lambda_r_min = 5635.0;
  double lambda_r_max = 7731.0;
  double lambda_z_min = 7445.0;
  double lambda_z_max = 9824.0;
  double lambda_min;
  double lambda_max;
  double lambda_hard_min;
  double lambda_hard_max;

  size_t spec_width = 0;
  size_t spec_overlap = 0;
  size_t spec_min = 0;
  size_t spec_max = 499;

  size_t res_width = 10;

  bool lambda_mask = false;
  
  string outdir = ".";
  string fibermap = "";
  string in_image = "";
  string in_psf = "";
  string in_sim = "";

  string band = "";
  string camera = "";
  string expid = "";

  string prefix = "harp:  ";

  bool quiet = false;
  bool debug = false;

  double global_start = MPI_Wtime();

  int gangsize = 0;

  // global communicator and grid

  boost::mpi::communicator comm;

  int np = comm.size();
  int myp = comm.rank();

  El::Grid grid ( El::mpi::COMM_WORLD );
  
  // Declare options
  
  popts::options_description desc ( "Allowed Options" );
  
  desc.add_options()
  ( "help,h", "display usage information" )
  ( "quiet,q", "supress information printing" )
  ( "debug,d", "write out intermediate data products for debugging" )
  ( "gangsize", popts::value < int > ( &gangsize ), "number of processes per gang (choose whole number of nodes and perfect square number of processes if possible)" )
  ( "image", popts::value < string > ( &in_image ), "input image file" )
  ( "psf", popts::value < string > ( &in_psf ), "input PSF file (XML)" )
  ( "psf_interp", "enable PSF interpolation" )
  ( "fibermap", popts::value < string > ( &fibermap ), "input fibermap" )
  ( "sim", popts::value < string > ( &in_sim ), "(optional) original simulated spectra" )
  ( "out", popts::value < string > ( &outdir ), "output directory (default = \".\")" )
  ( "spec_width", popts::value < size_t > ( &spec_width ), "number of spectra to solve at once" )
  ( "spec_overlap", popts::value < size_t > ( &spec_overlap ), "minimum spectra to overlap" )
  ( "spec_min", popts::value < size_t > ( &spec_min ), "first spectrum to solve (0)" )
  ( "spec_max", popts::value < size_t > ( &spec_max ), "last spectrum to solve (499)" )
  ( "lambda_width", popts::value < size_t > ( &lambda_width ), "number of wavelength points to solve at once" )
  ( "lambda_res", popts::value < double > ( &lambda_res ), "extraction resolution, in Angstroms" )
  ( "lambda_min", popts::value < double > ( &lambda_min ), "minimum extracted wavelength in Angstroms" )
  ( "lambda_max", popts::value < double > ( &lambda_max ), "maximum extracted wavelength in Angstroms" )
  ( "res_width", popts::value < size_t > ( &res_width ), "half-band width of diagonal resolution to keep (10)" )
  ;

  popts::variables_map vm;

  popts::store(popts::command_line_parser( argc, argv ).options(desc).run(), vm);
  
  popts::notify(vm);

  if ( ( ! vm.count("image") ) || ( ! vm.count("psf") ) || vm.count( "help" ) ) {
    if ( myp == 0 ) {
      cerr << endl;
      cerr << desc << endl;
    }
    El::Finalize();
    return 0;
  }
  
  if ( vm.count( "quiet" ) ) {
    quiet = true;
  }

  if ( vm.count( "debug" ) ) {
    debug = true;
  }

  size_t res_band = 2 * res_width + 1;

  // attempt to get metadata from the file name

  size_t baseoff = 0;
  size_t basetmp = in_image.find( '/', baseoff );
  while ( basetmp != string::npos ) {
    baseoff = basetmp + 1;
    basetmp = in_image.find( '/', baseoff );
  }

  size_t off = in_image.find ( '-', baseoff );
  if ( off != string::npos ) {
    off += 1;
    size_t pos = in_image.find ( '-', off );
    if ( pos != string::npos ) {
      band = in_image.substr ( off, 1 );
      camera = in_image.substr ( off, 2 );
      off = pos + 1;
      size_t ext = in_image.find ( '.', off );
      expid = in_image.substr ( off, ext - off );
    }
  }

  if ( ( myp == 0 ) && ( ! quiet ) ) {
    cout << prefix << "Image from camera \"" << camera << "\", exposure \"" << expid << "\"" << endl;
  }

  if ( ! vm.count( "lambda_min" ) ) {
    if ( band == "b" ) {
      lambda_min = lambda_b_min;
    } else if ( band == "r" ) {
      lambda_min = lambda_r_min;
    } else if ( band == "z" ) {
      lambda_min = lambda_z_min;
    } else {
      HARP_MPI_ABORT( myp, "image name does not conform to standard, and lambda min not specified" );
    }
  }

  if ( ! vm.count( "lambda_max" ) ) {
    if ( band == "b" ) {
      lambda_max = lambda_b_max;
    } else if ( band == "r" ) {
      lambda_max = lambda_r_max;
    } else if ( band == "z" ) {
      lambda_max = lambda_z_max;
    } else {
      HARP_MPI_ABORT( myp, "image name does not conform to standard, and lambda max not specified" );
    }
  }

  // we specify an extra overlap range in lambda, to remove
  // edge effects.

  size_t lambda_pre = (res_band * lambda_res);
  size_t lambda_post = (res_band * lambda_res);
  //size_t lambda_post = ((res_band + 1) * lambda_res);
  lambda_hard_min = lambda_min - lambda_pre;
  lambda_hard_max = lambda_max + lambda_post;

  // intra-gang communicator and grid

  if ( gangsize == 0 ) {
    gangsize = np;
  }

  int ngang = (int)( np / gangsize );
  int gangtot = ngang * gangsize;
  if ( ( myp == 0 ) && ( ! quiet ) ) {
    cout << prefix << "Using " << ngang << " gangs of " << gangsize << " processes each" << endl;
  }
  if ( gangtot < np ) {
    if ( myp == 0 ) {
      cout << prefix << "  WARNING: " << (np-gangtot) << " processes are idle" << endl;
    }
  }
  int gang = (int)( myp / gangsize );
  int grank = myp % gangsize;
  if ( gang >= ngang ) {
    gang = MPI_UNDEFINED;
    grank = MPI_UNDEFINED;
  }

  boost::mpi::communicator gcomm = comm.split ( gang, grank );

  El::Grid gang_grid ( (MPI_Comm)gcomm );

  // inter-gang communicator

  boost::mpi::communicator rcomm = comm.split ( grank, gang );


  // Create input globally-distributed PSF from the serial one

  if ( ( myp == 0 ) && ( ! quiet ) ) {
    cout << prefix << "Creating PSF..." << endl;
  }

  tstart = MPI_Wtime();

  boost::property_tree::ptree pathprops;
  boost::property_tree::ptree props;

  pathprops.clear();
  pathprops.put ( "path", in_psf );
  pathprops.put ( "wavemin", lambda_hard_min );
  pathprops.put ( "wavemax", lambda_hard_max );
  pathprops.put ( "wavebin", lambda_res );
  pathprops.put ( "fibermin", 0 );
  pathprops.put ( "fibermax", 499 );
  pathprops.put ( "type", "xml" );
  if ( vm.count( "psf_interp" ) ) {
    pathprops.put ( "interpolation", 1 );
  } else {
    pathprops.put ( "interpolation", 0 );
  }

  props.clear();
  props.put ( "type", "specex" );
  props.put_child ( "props", pathprops );

  mpi_psf_p design = mpi_load_psf ( comm, props );

  size_t psf_imgrows = design->img_rows();
  size_t psf_imgcols = design->img_cols();
  size_t psf_npix = psf_imgrows * psf_imgcols;

  size_t psf_nspec = design->n_spec();
  size_t psf_nlambda = design->n_lambda();
  size_t psf_nbins = psf_nspec * psf_nlambda;

  vector_double psf_lambda = design->lambda();

  size_t nlambda = psf_nlambda - 2*res_band;

  vector_double lambda(nlambda);
  for ( size_t i = 0; i < nlambda; ++i ) {
    lambda[i] = psf_lambda[i + res_band];
  }

  tstop = MPI_Wtime();

  if ( ( myp == 0 ) && ( ! quiet ) ) {
    cout << prefix << "  time = " << tstop-tstart << " seconds" << endl;
    cout << prefix << "  " << psf_nspec << " spectra each with " << psf_nlambda << " wavelength points" << endl;
    cout << prefix << "  image dimensions = " << psf_imgcols << " x " << psf_imgrows << endl;
  }

  // Read input image

  if ( ( myp == 0 ) && ( ! quiet ) ) {
    cout << prefix << "Creating input image..." << endl;
  }

  tstart = MPI_Wtime();

  pathprops.clear();
  pathprops.put ( "path", in_image );

  props.clear();
  props.put ( "type", "desi" );
  props.put_child ( "props", pathprops );

  mpi_image_p img = mpi_load_image ( comm, props );

  size_t imgrows = img->n_rows();
  size_t imgcols = img->n_cols();
  size_t npix = imgrows * imgcols;

  tstop = MPI_Wtime();

  if ( ( myp == 0 ) && ( ! quiet ) ) {
    cout << prefix << "  time = " << tstop-tstart << " seconds" << endl;
    cout << prefix << "  dimensions = " << imgcols << " x " << imgrows << endl;
  }

  if ( ( imgrows != psf_imgrows ) || ( imgcols != psf_imgcols ) ) {
    HARP_MPI_ABORT( myp, "PSF and image have inconsistent dimensions" );
  }

  if ( ( myp == 0 ) && ( ! quiet ) ) {
    cout << prefix << "Reading input image..." << endl;
  }

  tstart = MPI_Wtime();

  elem_matrix_local measured ( npix );
  local_matrix_zero ( measured );

  elem_matrix_local invnoise ( npix );
  local_matrix_zero ( invnoise );

  vector_mask msk;
  msk.clear();

  img->values ( measured );
  img->inv_variance ( invnoise );
  img->mask ( msk );

  tstop = MPI_Wtime();

  if ( ( myp == 0 ) && ( ! quiet ) ) {
    cout << prefix << "  time = " << tstop-tstart << " seconds" << endl;
  }


  // target list

  mpi_targets_p fibers;
  vector < object_p > target_list;

  if ( fibermap != "" ) {

    pathprops.clear();
    pathprops.put ( "path", fibermap );

    props.clear();
    props.put ( "type", "desi" );
    props.put_child ( "props", pathprops );

    fibers = mpi_load_targets ( comm, props );

    target_list = fibers->objects();

  }

  // select the global region of spec and wavelength space that we are
  // solving for.

  size_t global_first_spec = spec_min;
  size_t global_nspec = spec_max - spec_min + 1;

  size_t global_nbins = global_nspec * psf_nlambda;

  if ( ( myp == 0 ) && ( ! quiet ) ) {
    cout << prefix << "Solving globally for " << global_nspec << " spectra and " << psf_nlambda << " wavelength points" << endl;
  }

  // output spectral products

  mpi_matrix data_truth ( 0, 1, grid );
  mpi_matrix data_Rtruth ( 0, 1, grid );
  
  mpi_matrix data_f ( global_nbins, 1, grid );
  mpi_matrix_zero ( data_f );
  
  mpi_matrix data_Rf ( global_nbins, 1, grid );
  mpi_matrix_zero ( data_Rf );
  
  mpi_matrix data_err ( global_nbins, 1, grid );
  mpi_matrix_zero ( data_err );

  mpi_matrix data_Rdiag ( global_nbins, res_band, grid );
  mpi_matrix_zero ( data_Rdiag );


  // read truth spectra, if provided

  bool dotruth = false;

  if ( vm.count ( "sim" ) > 0 ) {

    if ( ( myp == 0 ) && ( ! quiet ) ) {
      cout << prefix << "Reading input truth spectra..." << endl;
    }

    tstart = MPI_Wtime();

    pathprops.clear();
    pathprops.put ( "path", in_sim );

    props.clear();
    props.put ( "type", "desisim" );
    props.put_child ( "props", pathprops );

    mpi_spec_p truth_spec = mpi_load_spec ( comm, props );
    size_t truth_nspec = truth_spec->n_spec();
    size_t truth_nlambda = truth_spec->n_lambda();
    size_t truth_nbins = truth_nspec * truth_nlambda;

    if ( truth_nbins != global_nbins ) {
      ostringstream o;
      o << "truth spectrum has " << truth_nbins << " spectral bins, but the PSF has " << global_nbins << " bins";
      HARP_MPI_ABORT( myp, o.str().c_str() );
    }

    data_truth.Resize ( global_nbins, 1 );
    mpi_matrix_zero ( data_truth );

    data_Rtruth.Resize ( global_nbins, 1 );
    mpi_matrix_zero ( data_Rtruth );

    vector_double truth_lambda;

    truth_spec->values ( data_truth );
    truth_spec->lambda ( truth_lambda );

    for ( size_t i = 0; i < nlambda; ++i ) {
      if ( fabs ( truth_lambda[i] - lambda[i] ) / lambda[i] > 1.0e-5 ) {
        ostringstream o;
        o << "wavelength point " << i << " does not match between PSF (" << lambda[i] << ") and truth spec (" << truth_lambda[i] << ")";
        HARP_MPI_ABORT( myp, o.str().c_str() );
      }
    }

    dotruth = true;

    tstop = MPI_Wtime();

    if ( ( myp == 0 ) && ( ! quiet ) ) {
      cout << prefix << "  time = " << tstop-tstart << " seconds" << endl;
    }

  }


  // Divide spectral data into overlapping regions

  if ( spec_width < 1 ) {
    // do all spectra
    spec_width = global_nspec;
  }

  if ( lambda_width < 1 ) {
    lambda_width = psf_nlambda;
  }

  mpi_spec_slice_p slice ( new mpi_spec_slice ( rcomm, gcomm, global_first_spec, res_band, global_nspec, nlambda, spec_width, lambda_width, spec_overlap, res_band ) );

  size_t gang_nregion = slice->regions().size();
  size_t nregion = 0;

  boost::mpi::reduce ( rcomm, gang_nregion, nregion, std::plus < size_t > (), 0 );

  if ( ( myp == 0 ) && ( ! quiet ) ) {
    cout << prefix << "Extracting " << nregion << " spectral chunks, each with " << spec_width << " spectra ( overlap = " << spec_overlap << " ) and " << lambda_width << " lambda points ( overlap = " << res_band << " )" << endl;
  }

  if ( debug && ( ! quiet ) ) {

    // each gang prints info about assigned slices

    if ( gcomm.rank() == 0 ) {

      vector < spec_slice_region > gslices = slice->regions();

      for ( size_t r = 0; r < rcomm.size(); ++r ) {

        if ( rcomm.rank() == r ) {

          cout << prefix << "  Gang " << r << " processing " << gslices.size() << " chunks:" << endl;
          for ( vector < spec_slice_region > :: const_iterator gsit = gslices.begin(); gsit != gslices.end(); ++gsit ) {
            cout << prefix << "    spec (" << gsit->first_good_spec << " - " << (gsit->first_good_spec + gsit->n_good_spec - 1) << ") wavelength (" << gsit->first_good_lambda << " - " << (gsit->first_good_lambda + gsit->n_good_lambda - 1) << ")" << endl;
          }

          cout.flush();

        }

        rcomm.barrier();

        cout.flush();

        rcomm.barrier();

      }

    }

  }


  // timing results

  map < string, double > timing;


  // Process all spectral slices

  string extract_prefix = "";
  if ( ! quiet ) {
    extract_prefix = prefix;
  }

  mpi_extract_slices ( slice, design, measured, invnoise, data_truth, res_band, data_Rdiag, data_Rf, data_f, data_err, data_Rtruth, timing, lambda_mask, prefix );

  /*
  string outres = outdir + "/debug_res.out";
  ofstream fres;
  fres.open(outres.c_str(), ios::out);
  El::Print ( data_Rdiag, "res_diag", fres );
  fres.close();
  */

  // subtract sky if needed



  // Total timings- get the mean, min and max from all processes

  if ( ( myp == 0 ) && ( ! quiet ) ) {
    cout << prefix << "Aggregate Timings:" << endl;
  }

  map < string, string > counter_desc;
  counter_desc [ "design" ] = "  Build design matrix =";
  counter_desc [ "inverse" ] = "  Build inverse covariance =";
  counter_desc [ "eigen" ] = "  Eigendecompose inverse =";
  counter_desc [ "norm" ] = "  Compute column norm =";
  counter_desc [ "nsespec" ] = "  Compute noise weighted spec =";
  counter_desc [ "extract" ] = "  Extract spectra =";

  double time_mean;
  double time_min;
  double time_max;

  for ( map < string, double > :: const_iterator itcounter = timing.begin(); itcounter != timing.end(); ++itcounter ) {

    time_mean = 0.0;
    time_min = 0.0;
    time_max = 0.0;

    boost::mpi::reduce ( comm, itcounter->second, time_mean, std::plus < double > (), 0 );
    time_mean /= (double)comm.size();

    boost::mpi::reduce ( comm, itcounter->second, time_min, boost::mpi::minimum < double > (), 0 );

    boost::mpi::reduce ( comm, itcounter->second, time_max, boost::mpi::maximum < double > (), 0 );

    if ( ( myp == 0 ) && ( ! quiet ) ) {
      cout << prefix << counter_desc [ itcounter->first ] << " " << time_mean << " ( " << time_min << " low / " << time_max << " high ) seconds" << endl;
    }

  }

  // Write outputs and compute chi-square if input truth is given

  if ( ( myp == 0 ) && ( ! quiet ) ) {
    cout << prefix << "Writing solution and error..." << endl;
  }

  tstart = MPI_Wtime();

  // get local copies.

  elem_matrix_local loc_truth;
  elem_matrix_local loc_Rtruth;
  elem_matrix_local loc_Rf;
  elem_matrix_local loc_f;
  elem_matrix_local loc_err;

  El::AxpyInterface < double > globloc;

  if ( dotruth ) {
    globloc.Attach( El::GLOBAL_TO_LOCAL, data_Rtruth );
    if ( myp == 0 ) {
      loc_Rtruth.Resize ( data_Rtruth.Height(), 1 );
      local_matrix_zero ( loc_Rtruth );
      globloc.Axpy ( 1.0, loc_Rtruth, 0, 0 );
    }
    globloc.Detach();

    globloc.Attach( El::GLOBAL_TO_LOCAL, data_truth );
    if ( myp == 0 ) {
      loc_truth.Resize ( data_truth.Height(), 1 );
      local_matrix_zero ( loc_truth );
      globloc.Axpy ( 1.0, loc_truth, 0, 0 );
    }
    globloc.Detach();
  }

  globloc.Attach( El::GLOBAL_TO_LOCAL, data_Rf );
  if ( myp == 0 ) {
    loc_Rf.Resize ( data_Rf.Height(), 1 );
    local_matrix_zero ( loc_Rf );
    globloc.Axpy ( 1.0, loc_Rf, 0, 0 );
  }
  globloc.Detach();

  globloc.Attach( El::GLOBAL_TO_LOCAL, data_f );
  if ( myp == 0 ) {
    loc_f.Resize ( data_f.Height(), 1 );
    local_matrix_zero ( loc_f );
    globloc.Axpy ( 1.0, loc_f, 0, 0 );
  }
  globloc.Detach();

  globloc.Attach( El::GLOBAL_TO_LOCAL, data_err );
  if ( myp == 0 ) {
    loc_err.Resize ( data_err.Height(), 1 );
    local_matrix_zero ( loc_err );
    globloc.Axpy ( 1.0, loc_err, 0, 0 );
  }
  globloc.Detach();

  boost::property_tree::ptree solution_props;
  string outfile;

  // define the clipping in wavelength to remove the extra overlap
  // region that we added at the beginning.

  size_t clip_nlambda = nlambda;
  //size_t clip_nlambda = psf_nlambda;

  vector_double & clip_lambda = lambda;
  //vector_double & clip_lambda = psf_lambda;
  
  size_t clip_lambda_start = res_band;
  //size_t clip_lambda_start = 0;

  size_t clip_lambda_stop = clip_lambda_start + clip_nlambda - 1;
  size_t clip_nbins = global_nspec * clip_nlambda;


  if ( myp == 0 ) {

    boost::property_tree::ptree child;

    solution_props.clear();

    // initialize to image metadata
    solution_props = boost::dynamic_pointer_cast < image_desi, image > ( img->local() )->meta();

    child.clear();
    child.put ( "VAL", clip_lambda[0] );
    child.put ( "TYPE", "F" );
    child.put ( "COM", "Starting wavelength [Angstroms]" );
    solution_props.put_child ( "CRVAL1", child );

    child.clear();
    child.put ( "VAL", lambda_res );
    child.put ( "TYPE", "F" );
    child.put ( "COM", "Wavelength step [Angstroms]" );
    solution_props.put_child ( "CDELT1", child );

    child.clear();
    child.put ( "VAL", "vac" );
    child.put ( "TYPE", "C" );
    child.put ( "COM", "Vacuum wavelengths" );
    solution_props.put_child ( "AIRORVAC", child );

    child.clear();
    child.put ( "VAL", 0 );
    child.put ( "TYPE", "I" );
    child.put ( "COM", "linear wavelength steps, not log10" );
    solution_props.put_child ( "LOGLAM", child );

    if ( dotruth ) {
      child.clear();
      child.put ( "VAL", in_sim );
      child.put ( "TYPE", "C" );
      child.put ( "COM", "Input simulation file" );
      solution_props.put_child ( "SIMFILE", child );
    }

    child.clear();
    child.put ( "VAL", in_psf );
    child.put ( "TYPE", "C" );
    child.put ( "COM", "Input spectral PSF" );
    solution_props.put_child ( "IN_PSF", child );

    child.clear();
    child.put ( "VAL", in_image );
    child.put ( "TYPE", "C" );
    child.put ( "COM", "Input image" );
    solution_props.put_child ( "IN_IMG", child );

    // write out solution

    vector_double errbuf;
    elem_to_ublas ( loc_err, errbuf );
    // convert from rms to inverse variance
    for ( size_t i = 0; i < errbuf.size(); ++i ) {
      errbuf[i] = 1.0 / ( errbuf[i] * errbuf[i] );
    }

    vector_double ubuf;
    elem_to_ublas ( loc_Rf, ubuf );

    //cerr << "clipping " << global_nspec << " spectra to lambda [0," << (psf_nlambda-1) << "] --> [" << clip_lambda_start << "," << (clip_lambda_stop) << "] (" << clip_nbins << " bins)" << endl; 

    vector_double clip_errbuf( clip_nbins );
    vector_double clip_ubuf( clip_nbins );

    for ( size_t i = 0; i < global_nspec; ++i ) {
      for ( size_t j = 0; j < clip_nlambda; ++j ) {
        clip_errbuf[i * clip_nlambda + j] = errbuf[i * psf_nlambda + clip_lambda_start + j];
        clip_ubuf[i * clip_nlambda + j] = ubuf[i * psf_nlambda + clip_lambda_start + j];
      }
    }

    outfile = outdir + "/frame-" + camera + "-" + expid + ".fits";
    spec_desi::write ( outfile, solution_props, clip_ubuf, clip_errbuf, clip_lambda );

    if ( dotruth ) {

      elem_to_ublas ( loc_Rtruth, ubuf );

      for ( size_t i = 0; i < global_nspec; ++i ) {
        for ( size_t j = 0; j < clip_nlambda; ++j ) {
          clip_ubuf[i * clip_nlambda + j] = ubuf[i * psf_nlambda + clip_lambda_start + j];
        }
      }

      outfile = outdir + "/harp_rsim-" + camera + "-" + expid + ".fits";
      spec_desi::write ( outfile, solution_props, clip_ubuf, clip_errbuf, clip_lambda );
      
      vector_double data_chisq ( clip_nbins );

      double chisq_reduced = 0.0;

      for ( size_t i = 0; i < global_nspec; ++i ) {
        for ( size_t j = 0; j < clip_nlambda; ++j ) {
          size_t w = i * psf_nlambda + clip_lambda_start + j;
          size_t v = i * clip_nlambda + j;
          if ( loc_err.Get(w,0) > std::numeric_limits < double > :: epsilon() ) {
            data_chisq[v] = ( loc_Rf.Get(w,0) - loc_Rtruth.Get(w,0) ) / loc_err.Get(w,0);
            data_chisq[v] *= data_chisq[v];
          } else {
            data_chisq[v] = 0.0;
          }
          chisq_reduced += data_chisq[v];
        }
      }

      chisq_reduced /= (double)( clip_nbins - 1 );

      if ( ! quiet ) {
        cout << prefix << "  Reduced Chi square = " << chisq_reduced << endl;
      }

      outfile = outdir + "/harp_chisq-" + camera + "-" + expid + ".fits";
      spec_desi::write ( outfile, solution_props, data_chisq, clip_errbuf, clip_lambda );

    }

  }

  // Free unused buffers before writing resolution matrix.

  loc_truth.Empty();
  loc_Rtruth.Empty();
  loc_Rf.Empty();
  loc_err.Empty();

  measured.Empty();
  data_Rf.Empty();
  data_Rtruth.Empty();
  data_err.Empty();

  // The resolution matrix is too large to reduce to a single process
  // all at once.  We have no choice but to reduce this in buffers and append
  // to the FITS file...

  fitsfile * fp;
  int ret;
  int status = 0;
  
  long naxes[3];
  naxes[0] = clip_nlambda;
  naxes[1] = res_band;
  naxes[2] = global_nspec;

  int fitstype = fits::ftype < double > :: datatype();
  long fpixel[2];

  elem_matrix_local loc_Rdiag;

  size_t write_spec_chunk = 10;
  size_t write_spec_offset = 0;

  if (write_spec_chunk > global_nspec) {
    write_spec_chunk = global_nspec;
  }

  // write buffer
  vector_double resbuffer;

  if ( myp == 0 ) {
    fits::open_readwrite ( fp, outfile );
    fits::img_seek ( fp, 3 );

    ret = fits_create_img ( fp, fits::ftype< double >::bitpix(), 3, naxes, &status );
    fits::check ( status );

    fits::key_write ( fp, "EXTNAME", string("RESOLUTION"), "" );
  }

  while ( write_spec_offset < global_nspec ) {
    if ( write_spec_offset + write_spec_chunk > global_nspec ) {
      write_spec_chunk = global_nspec - write_spec_offset;
    }

    size_t write_spec_nbuf = write_spec_chunk * res_band * clip_nlambda;

    if ( myp == 0 ) {
      resbuffer.resize(write_spec_nbuf);
      resbuffer.clear();
    }

    // we have to copy the data into the buffer one spec
    // at a time, due to memory layout.
    for ( size_t k = 0; k < write_spec_chunk; ++k ) {

      globloc.Attach( El::GLOBAL_TO_LOCAL, data_Rdiag );
      //cerr << "data_Rdiag has dimensions " << data_Rdiag.Height() << " x " << data_Rdiag.Width() << endl;
      if ( myp == 0 ) {
        loc_Rdiag.Resize ( clip_nlambda, res_band );
        local_matrix_zero ( loc_Rdiag );
        globloc.Axpy ( 1.0, loc_Rdiag, (write_spec_offset+k)*psf_nlambda + clip_lambda_start, 0 );
      }
      globloc.Detach();

      if ( myp == 0 ) {
        size_t boff = k * res_band * clip_nlambda;
        //cerr << "copying resolution spec " << (k + write_spec_offset) << "(" << boff << "-" << (boff + res_band * clip_nlambda) << ") to memory buffer" << endl;

        for ( size_t j = 0; j < res_band; ++j ) {

          for ( size_t i = 0; i < clip_nlambda; ++i ) {
          
            size_t outdiag = (res_band - 1) - j;
            int64_t outlambda = -1;
            int64_t inlambda = i;
            
            int64_t loff = (int64_t)j - (int64_t)res_width;

            outlambda = inlambda + loff;
            if ( outlambda >= (int64_t)clip_nlambda ) {
              outlambda = -1;
            }

            if ( outlambda >= 0 ) {
              double val = loc_Rdiag.Get(inlambda, j);
              if ( fabs(val) < 1.0e-100 ) {
                val = 0.0;
              }
              resbuffer[boff + outdiag * clip_nlambda + outlambda] = val;
            }
          }
        }
      }
    }

    if ( myp == 0 ) {
      fpixel[0] = 1;
      fpixel[1] = 1;
      fpixel[2] = write_spec_offset + 1;
      long npix = (long)( write_spec_nbuf );

      ret = fits_write_pix ( fp, fitstype, fpixel, npix, &(resbuffer[0]), &status );
      fits::check ( status );
    }

    write_spec_offset += write_spec_chunk;
  }

  if ( myp == 0 ) {
    fits::close( fp );
  }

  tstop = MPI_Wtime();

  if ( ( myp == 0 ) && ( ! quiet ) ) {
    cout << prefix << "  time = " << tstop-tstart << " seconds" << endl;
  }

  // free more memory

  loc_Rdiag.Empty();
  data_Rdiag.Empty();
  

  if ( debug ) {

    boost::property_tree::ptree child;
    child.clear();
    child.put ( "VAL", lambda[0] );
    child.put ( "TYPE", "F" );
    child.put ( "COM", "Starting wavelength [Angstroms]" );
    solution_props.erase("CRVAL1");
    solution_props.put_child ( "CRVAL1", child );

    if ( ( myp == 0 ) && ( ! quiet ) ) {
      cout << prefix << "Writing projected, deconvolved spectra..." << endl;
    }

    tstart = MPI_Wtime();

    map < size_t, set < size_t > > speclambda;

    //size_t slbins = 0;
    for ( size_t s = 0; s < global_nspec; ++s ) {
      for ( size_t l = 0; l < clip_nlambda; ++l ) {
        //cerr << "using spec " << (s + spec_min) << ", lambda " << (l + clip_lambda_start) << endl;
        speclambda[ s + spec_min ].insert ( l + clip_lambda_start );
        //++slbins;
      }
    }

    mpi_matrix_sparse AT ( comm, clip_nbins, psf_npix );

    design->project_transpose ( speclambda, AT );

    elem_matrix_local f_projected;

    elem_matrix_local newloc_f;
    newloc_f.Resize ( clip_nbins, 1 );
    local_matrix_zero ( newloc_f );

    for ( size_t i = 0; i < global_nspec; ++i ) {
      for ( size_t j = 0; j < clip_nlambda; ++j ) {
        newloc_f.Set( i * clip_nlambda + j, 0, loc_f.Get( i * psf_nlambda + clip_lambda_start + j, 0 ) );
      }
    }

    mpi_matrix new_f ( clip_nbins, 1, grid );
    mpi_matrix_zero ( new_f );

    globloc.Attach( El::LOCAL_TO_GLOBAL, new_f );
    if ( myp == 0 ) {
      globloc.Axpy ( 1.0, newloc_f, 0, 0 );
    }
    globloc.Detach();

    mpi_sparse_mv_trans ( AT, new_f, f_projected );

    vector_double u_f_projected;
    vector_double u_invnoise;

    if ( myp == 0 ) {

      elem_to_ublas ( f_projected, u_f_projected );
      f_projected.Empty();

      elem_to_ublas ( invnoise, u_invnoise );
      invnoise.Empty();

      outfile = outdir + "/harp_rpix-" + camera + "-" + expid + ".fits";
      image_desi::write ( outfile, solution_props, imgrows, u_f_projected, u_invnoise, msk );

    }

    loc_f.Empty();
    data_f.Empty();

    tstop = MPI_Wtime();

    if ( ( myp == 0 ) && ( ! quiet ) ) {
      cout << prefix << "  time = " << tstop-tstart << " seconds" << endl;
    }

    if ( dotruth ) {

      if ( ( myp == 0 ) && ( ! quiet ) ) {
        cout << prefix << "Writing projected truth spectra..." << endl;
      }

      tstart = MPI_Wtime();

      elem_matrix_local truth_projected;

      mpi_sparse_mv_trans ( AT, data_truth, truth_projected );

      // save more memory before allocating chisquare image

      AT.clear();
      data_truth.Empty();

      vector_double u_truth_projected;

      if ( myp == 0 ) {

        elem_to_ublas ( truth_projected, u_truth_projected );
        truth_projected.Empty();

        // pixel space chi-square

        vector_double data_chisq ( npix );

        double chisq_reduced = 0.0;

        for ( size_t i = 0; i < npix; ++i ) {
          if ( u_invnoise[i] > std::numeric_limits < double > :: epsilon() ) {
            data_chisq[i] = ( u_f_projected[i] - u_truth_projected[i] ) * sqrt( u_invnoise[i] );
            data_chisq[i] *= data_chisq[i];
          } else {
            data_chisq[i] = 0.0;
          }
          chisq_reduced += data_chisq[i];
        }

        chisq_reduced /= (double)( npix - 1 );

        if ( ! quiet ) {
          cout << prefix << "  Pixel space reduced Chi square = " << chisq_reduced << endl;
        }

        outfile = outdir + "/harp_rsimpix-" + camera + "-" + expid + ".fits";
        image_desi::write ( outfile, solution_props, imgrows, u_truth_projected, u_invnoise, msk );

        outfile = outdir + "/harp_chisqpix-" + camera + "-" + expid + ".fits";
        image_desi::write ( outfile, solution_props, imgrows, data_chisq, u_invnoise, msk );

      }

      tstop = MPI_Wtime();

      if ( ( myp == 0 ) && ( ! quiet ) ) {
        cout << prefix << "  time = " << tstop-tstart << " seconds" << endl;
      }

    }

  }

  // final timing

  double global_stop = MPI_Wtime();

  if ( ( myp == 0 ) && ( ! quiet ) ) {
    cout << prefix << "Total run time = " << global_stop-global_start << " seconds" << endl;
  }

  El::Finalize();

  return 0;
}



