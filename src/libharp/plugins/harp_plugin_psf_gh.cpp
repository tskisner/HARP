/*
  High Performance Astrophysical Reconstruction and Processing (HARP)
  (c) 2014-2015, The Regents of the University of California, 
  through Lawrence Berkeley National Laboratory.  See top
  level LICENSE file for details.
*/


#include <harp/static_plugins.hpp>

using namespace std;
using namespace harp;


// parameters

static const char * psf_gh_key_corr = "corr";
static const char * psf_gh_key_path = "path";


harp::psf_gh::psf_gh ( ) : psf () {
  nspec_ = 0;
  nlambda_ = 0;
  rows_ = 0;
  cols_ = 0;
  nglobal_ = 0;
  npix_ = 0;
  pixcorr_ = 0;
  path_ = "";
}


harp::psf_gh::psf_gh ( boost::property_tree::ptree const & props ) : psf ( "gh", props ) {

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


size_t harp::psf_gh::n_spec ( ) const {
  return nspec_;
}


size_t harp::psf_gh::n_lambda ( ) const {
  return nlambda_;
}


size_t harp::psf_gh::img_rows ( ) const {
  return rows_;
}


size_t harp::psf_gh::img_cols ( ) const {
  return cols_;
}


vector_double harp::psf_gh::lambda ( ) const {
  return lambda_;
}


void harp::psf_gh::extent ( size_t spec, size_t lambda, size_t & x_offset, size_t & y_offset, size_t & n_x, size_t & n_y ) const {
  return;
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


BOOST_CLASS_EXPORT(harp::psf_gh)

psf * harp::psf_gh_create ( boost::property_tree::ptree const & props ) {
  return new psf_gh ( props );
}
