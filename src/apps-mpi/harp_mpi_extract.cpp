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

  elem::Initialize ( argc, argv );

  boost::mpi::environment env;

  double tstart;
  double tstop;

  cout.precision ( 10 );
  cerr.precision ( 10 );

  size_t lambda_width = 0;
  size_t lambda_overlap = 0;

  size_t spec_width = 0;
  size_t spec_overlap = 0;

  bool lambda_mask = true;
  
  string outroot = "harp_";

  string jsonpar = "";

  string prefix = "harp:  ";

  bool quiet = false;
  bool debug = false;

  double global_start = MPI_Wtime();

  int gangsize = 0;

  // global communicator and grid

  boost::mpi::communicator comm;

  int np = comm.size();
  int myp = comm.rank();

  elem::Grid grid ( elem::mpi::COMM_WORLD );
  
  // Declare options
  
  popts::options_description desc ( "Allowed Options" );
  
  desc.add_options()
  ( "help,h", "display usage information" )
  ( "quiet,q", "supress information printing" )
  ( "debug,d", "write out intermediate data products for debugging" )
  ( "gangsize", popts::value < int > ( &gangsize ), "number of processes per gang (choose whole number of nodes and perfect square number of processes if possible)" )
  ( "out", popts::value<string>( &outroot ), "root path of output files (default = \"harp_\")" )
  ( "spec_width", popts::value < size_t > ( &spec_width ), "number of spectra to solve at once" )
  ( "spec_overlap", popts::value < size_t > ( &spec_overlap ), "minimum spectra to overlap" )
  ( "lambda_width", popts::value < size_t > ( &lambda_width ), "number of wavelength points to solve at once" )
  ( "lambda_overlap", popts::value < size_t > ( &lambda_overlap ), "minimum wavelength points to overlap" )
  ( "no_mask", "disable pixel masking in lambda overlap region" )
  ( "par", popts::value < string > ( &jsonpar ), "JSON parameter file" )
  ;

  popts::variables_map vm;

  popts::store(popts::command_line_parser( argc, argv ).options(desc).run(), vm);
  
  popts::notify(vm);

  if ( ( argc < 2 ) || vm.count( "help" ) || ( ! vm.count( "par" ) ) ) {
    if ( myp == 0 ) {
      cerr << endl;
      cerr << desc << endl;
    }
    elem::Finalize();
    return 0;
  }
  
  if ( vm.count( "quiet" ) ) {
    quiet = true;
  }

  if ( vm.count( "debug" ) ) {
    debug = true;
  }

  if ( vm.count( "no_mask" ) ) {
    lambda_mask = false;
  }

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

  elem::Grid gang_grid ( (MPI_Comm)gcomm );

  // inter-gang communicator

  boost::mpi::communicator rcomm = comm.split ( grank, gang );

  // Read metadata
  
  boost::property_tree::ptree par;

  if ( ( myp == 0 ) && ( ! quiet ) ) {
    cout << prefix << "Loading JSON parameter file..." << endl;
  }

  if ( myp == 0 ) {
    boost::property_tree::json_parser::read_json ( jsonpar, par );
  }

  boost::mpi::broadcast ( comm, par, 0 );

  // For now, we expect one psf and one image at the lop level of the
  // document.

  // Create input globally-distributed PSF

  if ( ( myp == 0 ) && ( ! quiet ) ) {
    cout << prefix << "Creating PSF..." << endl;
  }

  tstart = MPI_Wtime();

  boost::property_tree::ptree psf_props = par.get_child ( "psf" );

  mpi_psf_p design = mpi_load_psf ( comm, psf_props );

  size_t psf_imgrows = design->img_rows();
  size_t psf_imgcols = design->img_cols();
  size_t psf_npix = psf_imgrows * psf_imgcols;

  size_t psf_nspec = design->n_spec();
  size_t psf_nlambda = design->n_lambda();
  size_t psf_nbins = psf_nspec * psf_nlambda;

  vector_double lambda = design->lambda();

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

  boost::property_tree::ptree img_props = par.get_child ( "image" );

  mpi_image_p img = mpi_load_image ( comm, img_props );

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

  img->values ( measured );
  img->inv_variance ( invnoise );

  tstop = MPI_Wtime();

  if ( ( myp == 0 ) && ( ! quiet ) ) {
    cout << prefix << "  time = " << tstop-tstart << " seconds" << endl;
  }

  if ( debug ) {

    if ( ( myp == 0 ) && ( ! quiet ) ) {
      cout << prefix << "(debug mode) Dumping input image to disk..." << endl;
    }

    tstart = MPI_Wtime();

    if ( myp == 0 ) {

      string out_img_path = outroot + "debug_input_image.fits";

      vector_double temp_measured;
      vector_double temp_invnoise;

      elem_to_ublas ( measured, temp_measured );
      elem_to_ublas ( invnoise, temp_invnoise );

      image_fits::write ( out_img_path, imgrows, temp_measured, temp_invnoise );

    }

    tstop = MPI_Wtime();

    if ( ( myp == 0 ) && ( ! quiet ) ) {
      cout << prefix << "  time = " << tstop-tstart << " seconds" << endl;
    }

  }


  // target list

  // eventually we will read a target list and this will be used for sky subtraction and will
  // be written to the output spectral files.  For now, we create a fake list here.


  vector < object_p > target_list ( psf_nspec );

  for ( size_t i = 0; i < psf_nspec; ++i ) {
    target_list[i].reset ( new object ( OBJECT_UNKNOWN, "" ) );
  }


  // output spectral products

  mpi_matrix data_truth ( psf_nbins, 1, grid );
  mpi_matrix_zero ( data_truth );

  mpi_matrix data_Rtruth ( psf_nbins, 1, grid );
  mpi_matrix_zero ( data_Rtruth );
  
  mpi_matrix data_f ( psf_nbins, 1, grid );
  mpi_matrix_zero ( data_f );
  
  mpi_matrix data_Rf ( psf_nbins, 1, grid );
  mpi_matrix_zero ( data_Rf );
  
  mpi_matrix data_err ( psf_nbins, 1, grid );
  mpi_matrix_zero ( data_err );


  // read truth spectra, if provided

  bool dotruth = false;

  if ( par.count ( "truth" ) > 0 ) {

    if ( ( myp == 0 ) && ( ! quiet ) ) {
      cout << prefix << "Reading input truth spectra..." << endl;
    }

    tstart = MPI_Wtime();

    boost::property_tree::ptree truth_props = par.get_child ( "truth" );

    mpi_spec_p truth_spec = mpi_load_spec ( comm, truth_props );
    size_t truth_nspec = truth_spec->n_spec();
    size_t truth_nlambda = truth_spec->n_lambda();
    size_t truth_nbins = truth_nspec * truth_nlambda;

    if ( truth_nbins != psf_nbins ) {
      ostringstream o;
      o << "truth spectrum has " << truth_nbins << " spectral bins, but the PSF has " << psf_nbins << " bins";
      HARP_MPI_ABORT( myp, o.str().c_str() );
    }

    vector_double truth_lambda;

    truth_spec->values ( data_truth );
    truth_spec->lambda ( truth_lambda );

    for ( size_t i = 0; i < psf_nlambda; ++i ) {
      if ( fabs ( truth_lambda[i] - lambda[i] ) / lambda[i] > 1.0e-5 ) {
        ostringstream o;
        o << "wavelength point " << i << " does not match between PSF (" << lambda[i] << ") and truth spec (" << truth_lambda[i] << ")";
        HARP_MPI_ABORT( myp, o.str().c_str() );
      }
    }

    // FIXME: verify that target list matches input, after we actually read the input
    // target list.

    dotruth = true;

    tstop = MPI_Wtime();

    if ( ( myp == 0 ) && ( ! quiet ) ) {
      cout << prefix << "  time = " << tstop-tstart << " seconds" << endl;
    }

  }


  // Divide spectral data into overlapping regions

  if ( spec_width < 1 ) {
    // do all spectra
    spec_width = psf_nspec;
  }

  if ( lambda_width < 1 ) {
    lambda_width = psf_nlambda;
  }

  mpi_spec_slice_p slice ( new mpi_spec_slice ( rcomm, gcomm, psf_nspec, psf_nlambda, spec_width, lambda_width, spec_overlap, lambda_overlap ) );

  size_t gang_nregion = slice->regions().size();
  size_t nregion = 0;

  boost::mpi::reduce ( rcomm, gang_nregion, nregion, std::plus < size_t > (), 0 );

  if ( ( myp == 0 ) && ( ! quiet ) ) {
    cout << prefix << "Extracting " << nregion << " spectral chunks, each with " << spec_width << " spectra ( overlap = " << spec_overlap << " ) and " << lambda_width << " lambda points ( overlap = " << lambda_overlap << " )" << endl;
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

  mpi_extract_slices ( slice, design, measured, invnoise, data_truth, data_Rf, data_f, data_err, data_Rtruth, timing, lambda_mask, prefix );


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

  // get local copies

  elem_matrix_local loc_truth;
  elem_matrix_local loc_Rtruth;
  elem_matrix_local loc_Rf;
  elem_matrix_local loc_f;
  elem_matrix_local loc_err;

  elem::AxpyInterface < double > globloc;

  globloc.Attach( elem::GLOBAL_TO_LOCAL, data_Rtruth );
  if ( myp == 0 ) {
    loc_Rtruth.Resize ( data_Rtruth.Height(), 1 );
    local_matrix_zero ( loc_Rtruth );
    globloc.Axpy ( 1.0, loc_Rtruth, 0, 0 );
  }
  globloc.Detach();

  globloc.Attach( elem::GLOBAL_TO_LOCAL, data_truth );
  if ( myp == 0 ) {
    loc_truth.Resize ( data_truth.Height(), 1 );
    local_matrix_zero ( loc_truth );
    globloc.Axpy ( 1.0, loc_truth, 0, 0 );
  }
  globloc.Detach();

  globloc.Attach( elem::GLOBAL_TO_LOCAL, data_Rf );
  if ( myp == 0 ) {
    loc_Rf.Resize ( data_Rf.Height(), 1 );
    local_matrix_zero ( loc_Rf );
    globloc.Axpy ( 1.0, loc_Rf, 0, 0 );
  }
  globloc.Detach();

  globloc.Attach( elem::GLOBAL_TO_LOCAL, data_f );
  if ( myp == 0 ) {
    loc_f.Resize ( data_f.Height(), 1 );
    local_matrix_zero ( loc_f );
    globloc.Axpy ( 1.0, loc_f, 0, 0 );
  }
  globloc.Detach();

  globloc.Attach( elem::GLOBAL_TO_LOCAL, data_err );
  if ( myp == 0 ) {
    loc_err.Resize ( data_err.Height(), 1 );
    local_matrix_zero ( loc_err );
    globloc.Axpy ( 1.0, loc_err, 0, 0 );
  }
  globloc.Detach();

  boost::property_tree::ptree solution_props;
  string outfile;

  if ( myp == 0 ) {

    vector_double errbuf;
    elem_to_ublas ( loc_err, errbuf );

    vector_double ubuf;

    elem_to_ublas ( loc_Rf, ubuf );
    outfile = outroot + "spec_Rf.fits";
    spec_fits::write ( outfile, ubuf, errbuf, lambda );

    vector_double fake_err ( ubuf.size() );
    fake_err.clear();

    elem_to_ublas ( loc_f, ubuf );
    outfile = outroot + "spec_f.fits";
    spec_fits::write ( outfile, ubuf, fake_err, lambda );

    if ( dotruth ) {

      elem_to_ublas ( loc_Rtruth, ubuf );
      
      vector_double data_chisq ( psf_nbins );

      double chisq_reduced = 0.0;

      for ( size_t i = 0; i < psf_nbins; ++i ) {
        if ( loc_err.Get(i,0) > std::numeric_limits < double > :: epsilon() ) {
          data_chisq[i] = ( loc_Rf.Get(i,0) - loc_Rtruth.Get(i,0) ) / loc_err.Get(i,0);
          data_chisq[i] *= data_chisq[i];
        } else {
          data_chisq[i] = 0.0;
        }
        chisq_reduced += data_chisq[i];
      }

      chisq_reduced /= (double)( psf_nbins - 1 );

      if ( ! quiet ) {
        cout << prefix << "  Reduced Chi square = " << chisq_reduced << endl;
      }

      outfile = outroot + "spec_Rtruth.fits";
      spec_fits::write ( outfile, ubuf, data_chisq, lambda );

    }

  }

  tstop = MPI_Wtime();

  if ( ( myp == 0 ) && ( ! quiet ) ) {
    cout << prefix << "  time = " << tstop-tstart << " seconds" << endl;
  }

  // this next bit of code might (if debug == dotruth == true) generate 3 new images and the full PSF for
  // debugging purposes.  If we are tight on RAM, this could push us over the edge.  Explicitly deallocate
  // un-needed data objects at this point.

  loc_truth.Empty();
  loc_Rtruth.Empty();
  loc_Rf.Empty();
  loc_f.Empty();
  loc_err.Empty();

  measured.Empty();
  data_Rf.Empty();
  data_Rtruth.Empty();
  data_err.Empty();

  if ( debug ) {

    if ( ( myp == 0 ) && ( ! quiet ) ) {
      cout << prefix << "Writing projected, deconvolved spectra..." << endl;
    }

    tstart = MPI_Wtime();

    mpi_matrix_sparse AT ( comm, psf_nbins, npix );

    design->project_transpose ( AT );

    elem_matrix_local f_projected;

    mpi_sparse_mv_trans ( AT, data_f, f_projected );

    vector_double u_f_projected;
    vector_double u_invnoise;

    if ( myp == 0 ) {

      elem_to_ublas ( f_projected, u_f_projected );
      f_projected.Empty();

      elem_to_ublas ( invnoise, u_invnoise );
      invnoise.Empty();

      outfile = outroot + "image_f-project.fits";
      image_fits::write ( outfile, imgrows, u_f_projected, u_invnoise );

    }

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

        outfile = outroot + "image_truth-project.fits";
        image_fits::write ( outfile, imgrows, u_truth_projected, data_chisq );

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

  elem::Finalize();

  return 0;
}



