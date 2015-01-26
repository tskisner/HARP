/*
  High Performance Astrophysical Reconstruction and Processing (HARP)
  (c) 2014-2015, The Regents of the University of California, 
  through Lawrence Berkeley National Laboratory.  See top
  level LICENSE file for details.
*/


#include <harp/static_plugins.hpp>

#include <boost/random.hpp>

using namespace std;
using namespace harp;


static const char * image_sim_key_psf_type = "psf_type";
static const char * image_sim_key_spec_type = "spec_type";
static const char * image_sim_key_psf = "psf";
static const char * image_sim_key_spec = "spec";


harp::image_sim::image_sim ( ) : image () {
  rows_ = 0;
  cols_ = 0;
}
      

harp::image_sim::image_sim ( boost::property_tree::ptree const & props ) : image ( "sim", props ) {

  psf_type_ = props.get < string > ( image_sim_key_psf_type );
  spec_type_ = props.get < string > ( image_sim_key_spec_type );

  psf_props_ = props.get_child ( image_sim_key_psf );
  spec_props_ = props.get_child ( image_sim_key_spec );

  plugin_registry & reg = plugin_registry::get();

  spec_p child_spec ( reg.create_spec ( spec_type_, spec_props_ ) );

  size_t spec_nspec = child_spec->n_spec();
  size_t spec_nlambda = child_spec->n_lambda();

  psf_p child_psf ( reg.create_psf ( psf_type_, psf_props_ ) );

  size_t psf_nspec = child_psf->n_spec();
  size_t psf_nlambda = child_psf->n_lambda();

  if ( ( spec_nspec != psf_nspec ) || ( spec_nlambda != psf_nlambda ) ) {
    ostringstream o;
    o << "input spec size (" << spec_nspec << " x " << spec_nlambda << ") is not consistent with input PSF spectral size (" << psf_nspec << " x " << psf_nlambda << ")";
    HARP_THROW( o.str().c_str() );
  }

  rows_ = child_psf->img_rows();
  cols_ = child_psf->img_cols();

  size_t npix = rows_ * cols_;

  size_t nglobal = spec_nspec * spec_nlambda;

  // read spectra

  vector_double spec_data ( nglobal );
  vector_double spec_lambda ( spec_nlambda );

  child_spec->values ( spec_data );
  child_spec->lambda ( spec_lambda );

  // check wavelength solution against the one from the PSF

  vector_double psf_lambda = child_psf->lambda();

  for ( size_t i = 0; i < psf_nlambda; ++i ) {
    if ( fabs ( psf_lambda[i] - spec_lambda[i] ) / psf_lambda[i] > 0.001 ) {
      ostringstream o;
      o << "wavelength point " << i << " does not match between PSF (" << psf_lambda[i] << ") and spec (" << spec_lambda[i] << ")";
      HARP_THROW( o.str().c_str() );
    }
  }

  // get design matrix from PSF and project

  matrix_double_sparse AT;

  child_psf->project_transpose ( AT );

  // multiply to get signal data

  signal_.resize ( npix );

  sparse_mv_trans ( AT, spec_data, signal_ );

  // free up memory

  AT.clear();
  AT.resize ( 0, 0, false );

  // construct noise

  noise_.resize ( npix );
  invcov_.resize ( npix );

  typedef boost::ecuyer1988 base_generator_type;
  base_generator_type generator(42u);

  double rms;

  for ( size_t i = 0; i < npix; ++i ) {
    rms = sqrt( 16.0 + signal_[i] );
    invcov_[i] = 1.0 / (rms*rms);

    boost::normal_distribution < double > dist ( 0.0, rms );
    boost::variate_generator < base_generator_type&, boost::normal_distribution < double > > gauss ( generator, dist );

    noise_[i] = gauss();
  }
  
}


harp::image_sim::~image_sim ( ) {
  
}


void harp::image_sim::values ( vector_double & data ) const {

  data = signal_ + noise_;

  return;
}


void harp::image_sim::inv_variance ( vector_double & invvar ) const {

  invvar = invcov_;

  return;
}


BOOST_CLASS_EXPORT(harp::image_sim)

image * harp::image_sim_create ( boost::property_tree::ptree const & props ) {
  return new image_sim ( props );
}

