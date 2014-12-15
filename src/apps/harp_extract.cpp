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

#include <harp.hpp>


using namespace std;
using namespace harp;


namespace popts = boost::program_options;


int main ( int argc, char *argv[] ) {

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
  
  fitsfile * fp;

  double global_start = wtime();
  
  // Declare options
  
  popts::options_description desc ( "Allowed Options" );
  
  desc.add_options()
  ( "help,h", "display usage information" )
  ( "quiet,q", "supress information printing" )
  ( "debug,d", "write out intermediate data products for debugging" )
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
    cerr << endl;
    cerr << desc << endl;
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

  // Read metadata
  
  boost::property_tree::ptree par;

  if ( ! quiet ) {
    cout << prefix << "Loading JSON parameter file..." << endl;
  }

  boost::property_tree::json_parser::read_json ( jsonpar, par );

  // For now, we expect one psf and one image at the lop level of the
  // document.

  // Create input PSF

  if ( ! quiet ) {
    cout << prefix << "Creating PSF..." << endl;
  }

  tstart = wtime();

  boost::property_tree::ptree psf_props = par.get_child ( "psf" );

  psf_p design = load_psf ( psf_props );

  size_t psf_imgrows = design->img_rows();
  size_t psf_imgcols = design->img_cols();
  size_t psf_npix = psf_imgrows * psf_imgcols;

  size_t psf_nspec = design->n_spec();
  size_t psf_nlambda = design->n_lambda();
  size_t psf_nbins = psf_nspec * psf_nlambda;

  vector_double lambda = design->lambda();

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

  boost::property_tree::ptree img_props = par.get_child ( "image" );

  image_p img = load_image ( img_props );

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

    string out_img_path = outroot + "debug_input_image.fits";

    image_fits::write ( out_img_path, imgrows, measured, invnoise );

    tstop = wtime();

    if ( ! quiet ) {
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

  vector_double data_f ( psf_nbins );
  data_f.clear();

  vector_double data_Rf ( psf_nbins );
  data_Rf.clear();

  vector_double data_err ( psf_nbins );
  data_err.clear();

  vector_double data_truth ( psf_nbins );
  data_truth.clear();

  vector_double data_Rtruth ( psf_nbins );
  data_Rtruth.clear();


  // read truth spectra, if provided

  bool dotruth = false;

  if ( par.count ( "truth" ) > 0 ) {

    if ( ! quiet ) {
      cout << prefix << "Reading input truth spectra..." << endl;
    }

    tstart = wtime();

    boost::property_tree::ptree truth_props = par.get_child ( "truth" );

    spec_p truth_spec = load_spec ( truth_props );
    size_t truth_nspec = truth_spec->n_spec();
    size_t truth_nlambda = truth_spec->n_lambda();
    size_t truth_nbins = truth_nspec * truth_nlambda;

    if ( truth_nbins != psf_nbins ) {
      ostringstream o;
      o << "truth spectrum has " << truth_nbins << " spectral bins, but the PSF has " << psf_nbins << " bins";
      HARP_THROW( o.str().c_str() );
    }

    vector_double truth_lambda;

    truth_spec->values ( data_truth );
    truth_spec->lambda ( truth_lambda );

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


  // Divide spectral data into overlapping regions

  if ( spec_width < 1 ) {
    // do all spectra
    spec_width = psf_nspec;
  }

  if ( lambda_width < 1 ) {
    lambda_width = psf_nlambda;
  }

  spec_slice_p slice ( new spec_slice ( 1, psf_nspec, psf_nlambda, spec_width, lambda_width, spec_overlap, lambda_overlap ) );

  vector < spec_slice_region > regions = slice->regions ( 0 );

  if ( ! quiet ) {
    cout << prefix << "Extracting " << regions.size() << " spectral chunks, each with " << spec_width << " spectra ( overlap = " << spec_overlap << " ) and " << lambda_width << " lambda points ( overlap = " << lambda_overlap << " )" << endl;
  }

  if ( debug && ( ! quiet ) ) {

    for ( vector < spec_slice_region > :: const_iterator gsit = regions.begin(); gsit != regions.end(); ++gsit ) {
      cout << prefix << "    spec (" << gsit->first_good_spec << " - " << (gsit->first_good_spec + gsit->n_good_spec - 1) << ") wavelength (" << gsit->first_good_lambda << " - " << (gsit->first_good_lambda + gsit->n_good_lambda - 1) << ")" << endl;
    }

    cout.flush();

  }


  // timing results

  map < string, double > timing;


  // Process all spectral slices

  string extract_prefix = "";
  if ( ! quiet ) {
    extract_prefix = prefix;
  }

  extract_slices ( slice, design, measured, invnoise, data_truth, data_Rf, data_f, data_err, data_Rtruth, timing, false, lambda_mask, prefix );


  // subtract sky if needed



  // Total timings

  if ( ! quiet ) {
    cout << prefix << "Aggregate Timings:" << endl;
    cout << prefix << "  Build design matrix = " << timing["design"] << " seconds" << endl;
    cout << prefix << "  Build inverse covariance = " << timing["inverse"] << " seconds" << endl;
    cout << prefix << "  Eigendecompose inverse = " << timing["eigen"] << " seconds" << endl;
    cout << prefix << "  Compute column norm = " << timing["norm"] << " seconds" << endl;
    cout << prefix << "  Compute noise weighted spec = " << timing["nsespec"] << " seconds" << endl;
    cout << prefix << "  Extract spectra = " << timing["extract"] << " seconds" << endl;
  }


  // Write outputs and compute chi-square if input truth is given

  if ( ! quiet ) {
    cout << prefix << "Writing solution and error..." << endl;
  }

  tstart = wtime();

  string outfile = outroot + "spec_Rf.fits";
  spec_fits::write ( outfile, data_Rf, data_err, lambda );

  vector_double fake_err ( psf_nbins );
  fake_err.clear();

  outfile = outroot + "spec_f.fits";
  spec_fits::write ( outfile, data_f, fake_err, lambda );

  if ( dotruth ) {
    
    vector_double data_chisq ( psf_nbins );

    double chisq_reduced = 0.0;

    for ( size_t i = 0; i < psf_nbins; ++i ) {
      if ( data_err[i] > std::numeric_limits < double > :: epsilon() ) {
        data_chisq[i] = ( data_Rf[i] - data_Rtruth[i] ) / data_err[i];
        data_chisq[i] *= data_chisq[i];
      } else {
        data_chisq[i] = 0.0;
      }
      chisq_reduced += data_chisq[i];
    }

    chisq_reduced /= (double)( psf_nbins );

    if ( ! quiet ) {
      cout << prefix << "  Reduced Chi square = " << chisq_reduced << endl;
    }

    outfile = outroot + "spec_Rtruth.fits";
    spec_fits::write ( outfile, data_Rtruth, data_chisq, lambda );

  }

  tstop = wtime();

  if ( ! quiet ) {
    cout << prefix << "  time = " << tstop-tstart << " seconds" << endl;
  }

  // this next bit of code might (if debug == dotruth == true) generate 3 new images and the full PSF for
  // debugging purposes.  If we are tight on RAM, this could push us over the edge.  Explicitly deallocate
  // un-needed data objects at this point.

  measured.resize(0);
  data_Rf.resize(0);
  data_Rtruth.resize(0);
  data_err.resize(0);

  if ( debug ) {

    if ( ! quiet ) {
      cout << prefix << "Writing projected, deconvolved spectra..." << endl;
    }

    tstart = wtime();

    matrix_double_sparse AT;
    design->project_transpose ( AT );

    vector_double f_projected;
    sparse_mv_trans ( AT, data_f, f_projected );

    data_f.resize(0);

    outfile = outroot + "image_f-project.fits";
    image_fits::write ( outfile, imgrows, f_projected, invnoise );

    tstop = wtime();

    if ( ! quiet ) {
      cout << prefix << "  time = " << tstop-tstart << " seconds" << endl;
    }

    if ( dotruth ) {

      if ( ! quiet ) {
        cout << prefix << "Writing projected truth spectra..." << endl;
      }

      tstart = wtime();

      vector_double truth_projected;
      sparse_mv_trans ( AT, data_truth, truth_projected );

      // save more memory before allocating chisquare image

      AT.resize ( 0, 0, false );
      data_truth.resize(0);

      // pixel space chi-square

      vector_double data_chisq ( npix );

      double chisq_reduced = 0.0;

      for ( size_t i = 0; i < npix; ++i ) {
        if ( invnoise[i] > std::numeric_limits < double > :: epsilon() ) {
          data_chisq[i] = ( f_projected[i] - truth_projected[i] ) * sqrt( invnoise[i] );
          data_chisq[i] *= data_chisq[i];
        } else {
          data_chisq[i] = 0.0;
        }
        chisq_reduced += data_chisq[i];
      }

      chisq_reduced /= (double)( npix );

      if ( ! quiet ) {
        cout << prefix << "  Pixel space reduced Chi square = " << chisq_reduced << endl;
      }

      outfile = outroot + "image_truth-project.fits";
      image_fits::write ( outfile, imgrows, truth_projected, data_chisq );

      tstop = wtime();

      if ( ! quiet ) {
        cout << prefix << "  time = " << tstop-tstart << " seconds" << endl;
      }

    }

  }

  // final timing

  double global_stop = wtime();

  if ( ! quiet ) {
    cout << prefix << "Total run time = " << global_stop-global_start << " seconds" << endl;
  }

  return 0;
}


