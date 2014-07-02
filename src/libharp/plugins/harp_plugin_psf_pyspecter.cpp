// @COPYRIGHT@


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


harp::psf_pyspecter::psf_pyspecter ( boost::property_tree::ptree const & props ) : psf ( "pyspecter", props ) {

  path_ = props.get < string > ( psf_pyspecter_key_path );

  // instantiate the specter psf, get the pixel dimensions,
  // and pickle it for later use.

  Py_Initialize();

  try {

    py::object main_module = py::import("__main__");
    py::object main_namespace = main_module.attr("__dict__");

    py::exec ( "import specter", main_namespace );

    ostringstream com;
    com.str("");
    com << "pypsf = specter.load_psf(\"" << path_ << "\")";

    // instantiate psf
    py::exec ( com.str().c_str(), main_namespace );

    // number of spectra
    nspec_ = py::extract < size_t > ( py::eval ( "pypsf.nspec", main_namespace ) );

    // get image dims
    imgrows_ = py::extract < size_t > ( py::eval ( "pypsf.npix_y", main_namespace ) );
    imgcols_ = py::extract < size_t > ( py::eval ( "pypsf.npix_x", main_namespace ) );

    // get maximum wavelength extent that is valid for all spectra
    lambda_min_ = py::extract < size_t > ( py::eval ( "pypsf._wmin", main_namespace ) );
    lambda_max_ = py::extract < size_t > ( py::eval ( "pypsf._wmax", main_namespace ) );

    // to get an estimate of the NNZ, we just use the first point.

    py::exec ( 
      "xmin, xmax, ymin, ymax = pypsf.xyrange(0, 0)\n"
      "nnzapprox = (xmax-xmin)*(ymax-ymin)\n"
      , main_namespace );

    nnz_ = py::extract < size_t > ( py::eval ( "nnzapprox", main_namespace ) );

    // pickle this and return a string to c++
    py::exec ( "import pickle", main_namespace );
    pickled_ = py::extract < string > ( py::eval ( "pickle.dumps( pypsf, -1 )", main_namespace ) );

  } catch ( py::error_already_set ) {
    PyErr_Print();
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


// For all subsequent functions, unpickle local version of psf, and use it for operations

void harp::psf_pyspecter::extent ( size_t spec, size_t lambda, size_t & x_offset, size_t & y_offset, size_t & n_x, size_t & n_y ) const {



  return;
}


void harp::psf_pyspecter::response ( size_t spec, size_t lambda, size_t & x_offset, size_t & y_offset, matrix_double & patch ) const {

  HARP_THROW( "psf_pyspecter::response should never be called" );

  return;
}


void harp::psf_pyspecter::project ( std::map < size_t, std::set < size_t > > const & speclambda, matrix_double & A ) const {

  Py_Initialize();

  try {

    py::object main_module = py::import("__main__");
    py::object main_namespace = main_module.attr("__dict__");

    py::exec ( "import pickle", main_namespace );
    py::exec ( "import specter", main_namespace );

    // unpickle a new instance in python
    main_module.attr ( "picklestr" ) = pickled_;
    py::exec ( "pypsf = pickle.loads( picklestr )", main_namespace );






  } catch ( py::error_already_set ) {
    PyErr_Print();
  }


  return;
}


void harp::psf_pyspecter::project_transpose ( std::map < size_t, std::set < size_t > > const & speclambda, matrix_double_sparse & AT ) const {

  Py_Initialize();

  try {

    py::object main_module = py::import("__main__");
    py::object main_namespace = main_module.attr("__dict__");

    py::exec ( "import pickle", main_namespace );
    py::exec ( "import specter", main_namespace );

    // unpickle a new instance in python
    main_module.attr ( "picklestr" ) = pickled_;
    py::exec ( "pypsf = pickle.loads( picklestr )", main_namespace );




    

  } catch ( py::error_already_set ) {
    PyErr_Print();
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


