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

  double tstart;
  double tstop;

  cout.precision ( 10 );
  cerr.precision ( 10 );

  size_t lambda_width = 1;
  size_t lambda_overlap = 1;

  size_t spec_width = 1;
  size_t spec_overlap = 1;
  
  string outfile = "extracted.fits";

  string jsonpar = "";

  string prefix = "harp:  ";

  bool quiet = false;
  bool debug = false;
  
  fitsfile * fp;
  
  // Declare options
  
  popts::options_description desc ( "Allowed Options" );
  
  desc.add_options()
  ( "help,h", "display usage information" )
  ( "quiet,q", "supress information printing" )
  ( "debug,d", "write out intermediate data products for debugging" )
  ( "out", popts::value<string>( &outfile ), "output spectra file" )
  ( "spec_width", popts::value < size_t > ( &spec_width ), "number of spectra to solve at once" )
  ( "spec_overlap", popts::value < size_t > ( &spec_overlap ), "minimum spectra to overlap" )
  ( "lambda_width", popts::value < size_t > ( &lambda_width ), "number of wavelength points to solve at once" )
  ( "lambda_overlap", popts::value < size_t > ( &lambda_overlap ), "minimum wavelength points to overlap" )
  ( "par", popts::value < string > ( &jsonpar ), "JSON parameter file" )
  ;

  popts::variables_map vm;

  popts::store(popts::command_line_parser( argc, argv ).options(desc).run(), vm);
  
  popts::notify(vm);

  if ( ( argc < 2 ) || vm.count( "help" ) || ( ! vm.count( "par" ) ) ) {
    cerr << endl;
    cerr << desc << endl;
    return 0;
  }

  if ( 2 * lambda_overlap > lambda_width ) {
    cerr << prefix << "wavelength overlap is more than half the wavelength band!" << endl;
    exit(1);
  }
  
  if ( vm.count( "quiet" ) ) {
    quiet = true;
  }

  if ( vm.count( "debug" ) ) {
    debug = true;
  }


  // Read metadata
  
  boost::property_tree::ptree par;

  if ( ! quiet ) {
    cout << prefix << "Loading JSON parameter file..." << endl;
  }

  boost::property_tree::json_parser::read_json ( jsonpar, par );


  // Create input PSF

  if ( ! quiet ) {
    cout << prefix << "Creating PSF..." << endl;
  }

  tstart = wtime();

  boost::property_tree::ptree psf_props;
  psf_props = par.get_child ( "psf" );

  psf_p epsf ( psf::create ( psf_props ) );

  size_t psf_imgrows = epsf->img_rows();
  size_t psf_imgcols = epsf->img_cols();
  size_t psf_npix = psf_imgrows * psf_imgcols;

  size_t psf_nspec = epsf->n_spec();
  size_t psf_nlambda = epsf->n_lambda();
  size_t psf_nbins = psf_nspec * psf_nlambda;

  vector < double > lambda = epsf->lambda();

  tstop = wtime();

  if ( ! quiet ) {
    cout << prefix << "  time = " << tstop-tstart << " seconds" << endl;
    cout << prefix << "  " << psf_nspec << " spectra each with " << psf_nlambda << " wavelength points" << endl;
    cout << prefix << "  image dimensions = " << psf_imgcols << " x " << psf_imgrows << endl;
  }


  // Read input image

  if ( ! quiet ) {
    cout << prefix << "Creating input image..." << endl;
  }

  tstart = wtime();

  boost::property_tree::ptree img_props;
  img_props = conf.get_child ( "image" );

  image_p img ( image::create ( img_props ) );

  size_t imgrows = img->n_rows();
  size_t imgcols = img->n_cols();
  size_t npix = imgrows * imgcols;

  tstop = wtime();

  if ( ! quiet ) {
    cout << prefix << "  time = " << tstop-tstart << " seconds" << endl;
    cout << prefix << "  dimensions = " << imgcols << " x " << imgrows << endl;
  }

  if ( ( imgrows != psf_imgrows ) || ( imgcols != psf_imgcols ) ) {
    HARP_THROW( "PSF and image have inconsistent dimensions" );
    exit(1);
  }

  if ( ! quiet ) {
    cout << prefix << "Reading input image..." << endl;
  }

  tstart = wtime();

  vector_double measured ( npix );
  measured.clear();

  vector_double invnoise ( npix );
  invnoise.clear();

  img->values ( measured );
  img->inv_variance ( invnoise );

  tstop = wtime();

  if ( ! quiet ) {
    cout << prefix << "  time = " << tstop-tstart << " seconds" << endl;
  }

  if ( debug ) {

    if ( ! quiet ) {
      cout << prefix << "(debug mode) Dumping input image to disk..." << endl;
    }

    tstart = wtime();

    string out_img_path = "harp_debug_input_image.fits";

    boost::property_tree::ptree out_img_props;

    out_img_props.clear();
    out_img_props.put ( "format", "fits" );
    out_img_props.put ( "rows", imgrows );
    out_img_props.put ( "cols", imgcols );

    image_fits outimg ( out_img_props );

    outimg.write ( out_img_path, measured, invnoise );

    tstop = wtime();

    if ( ! quiet ) {
      cout << prefix << "  time = " << tstop-tstart << " seconds" << endl;
    }

  }

  // Divide spectral data into overlapping regions

  spec_slice_p slice ( new spec_slice ( 1, psf_nspec, psf_nlambda, spec_width, lambda_width, spec_overlap, lambda_overlap ) );

  vector < spec_slice_region > regions = slice->regions ( 0 );

  if ( ! quiet ) {
    cout << prefix << "Extracting " << regions.size() << " spectral chunks, each with " << spec_width << " spectra and " << lambda_width << " lambda points ( overlap = " << lambda_overlap << " )" << endl;
  }

  // output spectral products

  vector_double full_f ( psf_nbins );
  full_f.clear();

  vector_double full_Rf ( psf_nbins );
  full_Rf.clear();

  vector_double full_err ( psf_nbins );
  full_err.clear();

  vector_double full_truth ( psf_nbins );
  full_truth.clear();

  vector_double full_Rtruth ( psf_nbins );
  full_Rtruth.clear();

  // read truth spectra, if provided

  bool dotruth = false;

  if ( conf.count ( "truth" ) > 0 ) {

    if ( ! quiet ) {
      cout << prefix << "Reading input truth spectra..." << endl;
    }

    tstart = wtime();

    boost::property_tree::ptree truth_props = conf.get_child ( "truth" );

    spec_p truth_spec ( spec::create ( truth_props ) );
    size_t truth_nspec = truth_spec->n_spec();
    size_t truth_nlambda = truth_spec->n_lambda();
    size_t truth_nbins = truth_nspec * truth_nlambda;

    if ( truth_nbins != psf_nbins ) {
      ostringstream o;
      o << "truth spectrum has " << truth_nbins << " spectral bins, but the PSF has " << psf_nbins << " bins";
      HARP_THROW( o.str().c_str() );
    }

    vector_double truth_lambda;
    vector < target > truth_target_list;

    testspec->values ( full_truth );
    testspec->lambda ( truth_lambda );
    testspec->targets ( truth_target_list );

    for ( size_t i = 0; i < psf_nlambda; ++i ) {
      if ( fabs ( truth_lambda[i] - lambda[i] ) / lambda[i] > 1.0e-5 ) {
        ostringstream o;
        o << "wavelength point " << i << " does not match between PSF (" << lambda[i] << ") and truth spec (" << truth_lambda[i] << ")";
        HARP_THROW( o.str().c_str() );
      }
    }

    // FIXME: verify that target list matches input, after we actually read the input
    // target list.

    dotruth = true;

    tstop = wtime();

    if ( ! quiet ) {
      cout << prefix << "  time = " << tstop-tstart << " seconds" << endl;
    }

  }

  // timing

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

  // Process all spectral slices

  for ( vector < spec_slice_region > :: iterator regit = regions.begin(); regit != regions.end(); ++regit ) {

    tstart = wtime();
    double tsubstart = tstart;

    size_t nbins = regit->n_spec * regit->n_lambda;

    // build the list of spectral points we want for the projection

    std::map < size_t, std::set < size_t > > speclambda;

    for ( size_t s = 0; s < regit->n_spec; ++s ) {
      for ( size_t l = 0; l < regit->n_lambda; ++l ) {
        speclambda[ s + regit->first_spec ].insert ( l + regit->first_lambda );
      }
    }

    // get the projection matrix for this slice

    matrix_double_sparse AT;

    epsf->project_transpose ( speclambda, AT );

    double tsubstop = wtime();
    time_design = ( tsubstop - tsubstart );

    // build the inverse spectral covariance for this slice

    matrix_double inv ( nbins, nbins );

    




    size_t overlap_spec;
      size_t overlap_lambda;
      size_t first_spec;
      size_t first_lambda;
      size_t first_good_spec;
      size_t first_good_lambda;
      size_t n_spec;
      size_t n_lambda;
      size_t n_good_spec;
      size_t n_good_lambda;



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


