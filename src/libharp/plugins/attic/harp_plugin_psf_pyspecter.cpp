/*
  High Performance Astrophysical Reconstruction and Processing (HARP)
  (c) 2014-2015, The Regents of the University of California, 
  through Lawrence Berkeley National Laboratory.  See top
  level LICENSE file for details.
*/


#include <harp/static_plugins.hpp>

using namespace std;
using namespace harp;


#include <config.h>

// parameters

static const char * psf_pyspecter_key_path = "path";
static const char * psf_pyspecter_key_wmin = "wavemin";
static const char * psf_pyspecter_key_wmax = "wavemax";
static const char * psf_pyspecter_key_wbin = "wavebin";


#ifdef HAVE_BOOST_PYTHON_HPP

#include <boost/python.hpp>
#include <boost/numpy.hpp>

namespace py = boost::python;


harp::psf_pyspecter::psf_pyspecter ( ) : psf () {
}


harp::psf_pyspecter::psf_pyspecter ( boost::property_tree::ptree const & props ) : psf ( "pyspecter", props ) {

  path_ = props.get < string > ( psf_pyspecter_key_path );

  // instantiate the specter psf, get the pixel dimensions,
  // and pickle it for later use.

  Py_Initialize();

  // Initialize NumPy
  boost::numpy::initialize();

  try {

    py::object main_module = py::import("__main__");
    py::object main_namespace = main_module.attr("__dict__");

    py::exec ( "import specter.psf as specpsf", main_namespace );

    ostringstream com;
    com.str("");
    com << "pypsf = specpsf.load_psf(\"" << path_ << "\")";

    // instantiate psf
    py::exec ( com.str().c_str(), main_namespace );

    // number of spectra
    nspec_ = py::extract < size_t > ( py::eval ( "pypsf.nspec", main_namespace ) );

    // get image dims
    imgrows_ = py::extract < size_t > ( py::eval ( "pypsf.npix_y", main_namespace ) );
    imgcols_ = py::extract < size_t > ( py::eval ( "pypsf.npix_x", main_namespace ) );

    // get maximum wavelength extent that is valid for all spectra
    lambda_min_ = py::extract < double > ( py::eval ( "pypsf._wmin", main_namespace ) );
    lambda_max_ = py::extract < double > ( py::eval ( "pypsf._wmax", main_namespace ) );

    // to get an estimate of the NNZ, we just use the first point.

    py::exec ( 
      "xmin, xmax, ymin, ymax = pypsf.xyrange(0, pypsf._wmin)\n"
      "nnzapprox = (xmax-xmin)*(ymax-ymin)\n"
      , main_namespace );

    nnz_ = py::extract < size_t > ( py::eval ( "nnzapprox", main_namespace ) );

    // pickle this and return a string to c++
    py::exec ( "import pickle", main_namespace );
    pickled_ = py::extract < string > ( py::eval ( "pickle.dumps( pypsf, -1 )", main_namespace ) );

  } catch ( py::error_already_set ) {

    string error_str = parse_python_exception();
    throw std::runtime_error ( error_str );

    //PyErr_Print();
  }

  wavemin_ = props.get < double > ( psf_pyspecter_key_wmin, lambda_min_ );
  wavemax_ = props.get < double > ( psf_pyspecter_key_wmax, lambda_max_ );
  wavebin_ = props.get < double > ( psf_pyspecter_key_wbin, 1.0 );

  nlambda_ = (size_t)( ( wavemax_ - wavemin_ ) / wavebin_ ) + 1;
  
  lambda_.resize ( nlambda_ );
    
  for ( size_t i = 0; i < nlambda_; ++i ) {
    lambda_[i] = wavemin_ + wavebin_ * (double)i;
  }

  // global dimensions

  nglobal_ = nspec_ * nlambda_;

  npix_ = imgrows_ * imgcols_;

}


harp::psf_pyspecter::~psf_pyspecter ( ) {
  
}


size_t harp::psf_pyspecter::n_spec ( ) const {
  return nspec_;
}


size_t harp::psf_pyspecter::n_lambda ( ) const {
  return nlambda_;
}


size_t harp::psf_pyspecter::img_rows ( ) const {
  return imgrows_;
}


size_t harp::psf_pyspecter::img_cols ( ) const {
  return imgcols_;
}


vector_double harp::psf_pyspecter::lambda ( ) const {
  return lambda_;
}


size_t harp::psf_pyspecter::response_nnz_estimate ( ) const {
  return nnz_;
}


// For all subsequent functions, unpickle local version of psf, and use it for operations

void harp::psf_pyspecter::extent ( size_t spec, size_t lambda, size_t & x_offset, size_t & y_offset, size_t & n_x, size_t & n_y ) const {

  HARP_THROW( "psf_pyspecter::extent should never be called" );

  return;
}


void harp::psf_pyspecter::response ( size_t spec, size_t lambda, size_t & x_offset, size_t & y_offset, matrix_double & patch ) const {

  HARP_THROW( "psf_pyspecter::response should never be called" );

  return;
}


void harp::psf_pyspecter::extent_multi ( std::map < size_t, std::set < size_t > > const & speclambda, std::vector < size_t > & x_offset, std::vector < size_t > & y_offset, std::vector < size_t > & n_x, std::vector < size_t > & n_y ) const {

  size_t total = psf::total_bins ( speclambda );

  x_offset.resize ( total );
  y_offset.resize ( total );
  n_x.resize ( total );
  n_y.resize ( total );

  Py_Initialize();

  try {

    py::object main_module = py::import("__main__");
    py::object main_namespace = main_module.attr("__dict__");

    py::exec ( "import pickle", main_namespace );
    py::exec ( "import specter.psf as specpsf", main_namespace );

    // unpickle a new instance in python
    main_module.attr ( "picklestr" ) = pickled_;
    py::exec ( "pypsf = pickle.loads( picklestr )", main_namespace );

    size_t cur = 0;

    for ( std::map < size_t, std::set < size_t > > :: const_iterator itspec = speclambda.begin(); itspec != speclambda.end(); ++itspec ) {

      for ( std::set < size_t > :: const_iterator itlambda = itspec->second.begin(); itlambda != itspec->second.end(); ++itlambda ) {

        double wave = lambda_[ (*itlambda) ];
        main_module.attr ( "ext_spec" ) = itspec->first;
        main_module.attr ( "ext_wave" ) = wave;

        py::exec ( "xmin, xmax, ymin, ymax = pypsf.xyrange(ext_spec, ext_wave)", main_namespace );

        size_t xmin = py::extract < size_t > ( py::eval ( "xmin", main_namespace ) );
        size_t xmax = py::extract < size_t > ( py::eval ( "xmax", main_namespace ) );
        size_t ymin = py::extract < size_t > ( py::eval ( "ymin", main_namespace ) );
        size_t ymax = py::extract < size_t > ( py::eval ( "ymax", main_namespace ) );

        size_t nx = xmax - xmin + 1;
        size_t ny = ymax - ymin + 1;

        x_offset[ cur ] = xmin;
        y_offset[ cur ] = ymin;
        n_x[ cur ] = nx;
        n_y[ cur ] = ny;

        ++cur;

      }

    }

  } catch ( py::error_already_set ) {

    string error_str = parse_python_exception();
    throw std::runtime_error ( error_str );

    //PyErr_Print();
  }

  return;
}


void harp::psf_pyspecter::project ( std::map < size_t, std::set < size_t > > const & speclambda, matrix_double & A ) const {

  size_t total = psf::total_bins ( speclambda );

  // resize output to correct dimensions

  A.resize ( npix_, total, false );
  A.clear();

  // iterate over spectral bins and populate the matrix elements

  Py_Initialize();

  try {

    py::object main_module = py::import("__main__");
    py::object main_namespace = main_module.attr("__dict__");

    py::exec ( "import pickle", main_namespace );
    py::exec ( "import numpy", main_namespace );
    py::exec ( "import specter.psf as specpsf", main_namespace );

    // unpickle a new instance in python
    main_module.attr ( "picklestr" ) = pickled_;
    py::exec ( "pypsf = pickle.loads( picklestr )", main_namespace );

    matrix_double patch;
    size_t xoff;
    size_t yoff;

    size_t col = 0;
    size_t row;

    for ( std::map < size_t, std::set < size_t > > :: const_iterator itspec = speclambda.begin(); itspec != speclambda.end(); ++itspec ) {

      for ( std::set < size_t > :: const_iterator itlambda = itspec->second.begin(); itlambda != itspec->second.end(); ++itlambda ) {

        double wave = lambda_[ (*itlambda) ];
        main_module.attr ( "ext_spec" ) = itspec->first;
        main_module.attr ( "ext_wave" ) = wave;

        py::exec ( "xmin, xmax, ymin, ymax = pypsf.xyrange(ext_spec, ext_wave)", main_namespace );
        py::exec ( "xslice, yslice, pixels = pypsf.xypix(ext_spec, ext_wave, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)", main_namespace );

        // pass back to C++

        xoff = py::extract < size_t > ( py::eval ( "xslice[0]", main_namespace ) );
        yoff = py::extract < size_t > ( py::eval ( "yslice[0]", main_namespace ) );

        boost::numpy::ndarray patch = py::extract < boost::numpy::ndarray > ( py::eval ( "pixels", main_namespace ) );

        int patchdims = patch.get_nd();

        if ( patchdims != 2 ) {
          HARP_THROW( "PSF patch is not a 2D array" );
          return;
        }

        Py_intptr_t const * patchshape = patch.get_shape();

        double * pydata = reinterpret_cast < double * > ( patch.get_data() );

        for ( size_t patch_col = 0; patch_col < patchshape[1]; ++patch_col ) {

          if ( xoff + patch_col < imgcols_ ) {
            // this column is within the image dimensions

            for ( size_t patch_row = 0; patch_row < patchshape[0]; ++patch_row ) {

              if ( yoff + patch_row < imgrows_ ) {
                // this row is within the image dimensions

                row = ( xoff + patch_col ) * imgrows_ + yoff + patch_row;

                A ( row, col ) = pydata[ patch_row * patchshape[1] + patch_col ];

              }

            }

          }

        }

        ++col;

      }

    }

  } catch ( py::error_already_set ) {

    string error_str = parse_python_exception();
    throw std::runtime_error ( error_str );

    //PyErr_Print();
  }

  return;
}


void harp::psf_pyspecter::project_transpose ( std::map < size_t, std::set < size_t > > const & speclambda, matrix_double_sparse & AT ) const {

  size_t total = total_bins ( speclambda );

  // resize output to correct dimensions

  AT.resize ( total, npix_, false );
  AT.clear();

  // we use the estimated number of non-zeros per row times the number 
  // of rows to reserve space in the sparse matrix.

  size_t nnz_estimate = response_nnz_estimate();
  nnz_estimate *= total;

  AT.reserve ( nnz_estimate, false );

  // iterate over spectral bins and populate the matrix elements

  Py_Initialize();

  try {

    py::object main_module = py::import("__main__");
    py::object main_namespace = main_module.attr("__dict__");

    py::exec ( "import pickle", main_namespace );
    py::exec ( "import specter.psf as specpsf", main_namespace );

    // unpickle a new instance in python
    main_module.attr ( "picklestr" ) = pickled_;
    py::exec ( "pypsf = pickle.loads( picklestr )", main_namespace );

    matrix_double patch;
    size_t xoff;
    size_t yoff;

    size_t row = 0;
    size_t col;

    for ( std::map < size_t, std::set < size_t > > :: const_iterator itspec = speclambda.begin(); itspec != speclambda.end(); ++itspec ) {

      for ( std::set < size_t > :: const_iterator itlambda = itspec->second.begin(); itlambda != itspec->second.end(); ++itlambda ) {

        double wave = lambda_[ (*itlambda) ];
        main_module.attr ( "ext_spec" ) = itspec->first;
        main_module.attr ( "ext_wave" ) = wave;

        py::exec ( "xmin, xmax, ymin, ymax = pypsf.xyrange(ext_spec, ext_wave)", main_namespace );
        py::exec ( "xslice, yslice, pixels = pypsf.xypix(ext_spec, ext_wave, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)", main_namespace );

        // pass back to C++

        xoff = py::extract < size_t > ( py::eval ( "xslice[0]", main_namespace ) );
        yoff = py::extract < size_t > ( py::eval ( "yslice[0]", main_namespace ) );

        boost::numpy::ndarray patch = py::extract < boost::numpy::ndarray > ( py::eval ( "pixels", main_namespace ) );

        int patchdims = patch.get_nd();

        if ( patchdims != 2 ) {
          HARP_THROW( "PSF patch is not a 2D array" );
          return;
        }

        Py_intptr_t const * patchshape = patch.get_shape();

        double * pydata = reinterpret_cast < double * > ( patch.get_data() );

        for ( size_t patch_col = 0; patch_col < patchshape[1]; ++patch_col ) {

          if ( xoff + patch_col < imgcols_ ) {
            // this column is within the image dimensions

            for ( size_t patch_row = 0; patch_row < patchshape[0]; ++patch_row ) {

              if ( yoff + patch_row < imgrows_ ) {
                // this row is within the image dimensions

                col = ( xoff + patch_col ) * imgrows_ + yoff + patch_row;

                AT ( row, col ) = pydata[ patch_row * patchshape[1] + patch_col ];

              }

            }

          }

        }

        ++row;

      }

    }

  } catch ( py::error_already_set ) {

    string error_str = parse_python_exception();
    throw std::runtime_error ( error_str );

    //PyErr_Print();
  }


  return;
}


#else

harp::psf_pyspecter::psf_pyspecter ( boost::property_tree::ptree const & props ) : psf ( "pyspecter", props ) {
  HARP_THROW( "HARP not compiled with python support, cannot instantiate pyspecter PSF" );
}


harp::psf_pyspecter::~psf_pyspecter ( ) { 
}

void harp::psf_pyspecter::extent ( size_t spec, size_t lambda, size_t & x_offset, size_t & y_offset, size_t & n_x, size_t & n_y ) const {
  return;
}


void harp::psf_pyspecter::response ( size_t spec, size_t lambda, size_t & x_offset, size_t & y_offset, matrix_double & patch ) const {
  return;
}

void harp::psf_pyspecter::extent_multi ( std::map < size_t, std::set < size_t > > const & speclambda, std::vector < size_t > & x_offset, std::vector < size_t > & y_offset, std::vector < size_t > & n_x, std::vector < size_t > & n_y ) const {
  return;
}

void harp::psf_pyspecter::project ( std::map < size_t, std::set < size_t > > const & speclambda, matrix_double & A ) const {
  return;
}


void harp::psf_pyspecter::project_transpose ( std::map < size_t, std::set < size_t > > const & speclambda, matrix_double_sparse & AT ) const {
  return;
}

#endif


BOOST_CLASS_EXPORT(harp::psf_pyspecter)

psf * harp::psf_pyspecter_create ( boost::property_tree::ptree const & props ) {
  return new psf_pyspecter ( props );
}


