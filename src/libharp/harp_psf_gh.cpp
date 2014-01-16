// @COPYRIGHT@

#include <harp_internal.hpp>

#include <boost/optional.hpp>

#include <cmath>

using namespace std;
using namespace harp;




// parameters

static const char * psf_gh_key_corr = "corr";
static const char * psf_gh_key_path = "path";


harp::psf_gh::psf_gh ( boost::property_tree::ptree const & props ) : psf ( props ) {

  pixcorr_ = props.get < int > ( psf_gh_key_corr );

  path_ = props.get < string > ( psf_gh_key_path, "" );

  if ( path_ == "" ) {

    // we are solving for the PSF

    // parse other options here to instantiate the calibration images and other information.
    // from this, determine the number of spectra, image size and wavelength solution.
    // Solve for the PSF.

    // NOTE:  we can pass in properties needed to instantiate these images as child nodes of
    // props.  See the constructor in harp_image_sim.cpp as an example.

    // nspec_ = ???;
    // nlambda_ = ???;
    // rows_ = ???;
    // cols_ = ???;
    // 
    // lambda_.resize ( nlambda_ );
    // ... set elements of lambda_ ...
    //


  } else {

    // we are loading the PSF from a file





  }

  // global dimensions

  nglobal_ = nspec_ * nlambda_;

  npix_ = rows_ * cols_;

}


harp::psf_gh::~psf_gh ( ) {
  
}


boost::property_tree::ptree harp::psf_gh::metadata ( ) const {

  return boost::property_tree::ptree();
}


size_t harp::psf_gh::response_nnz_estimate ( ) const {
  size_t est = 2 * pixcorr_ + 1;
  est *= est;
  return est;
}


void harp::psf_gh::response ( size_t spec, size_t lambda, size_t & x_offset, size_t & y_offset, matrix_double & patch ) const {

  size_t bin = spec * nlambda_ + lambda;

  // for this spectral bin, use pixcorr_ to determine the size of the patch in pixel space that should
  // be sampled.  resize patch and set x_offset and y_offset.  populate the elements of patch.



  return;
}


void harp::psf_gh::write ( std::string const & path ) {

  fitsfile *fp;

  fits::create ( fp, path );

  // write out the HDUs with information needed to read it back in...




  fits::close ( fp );

  return;
}





