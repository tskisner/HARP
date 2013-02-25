// @COPYRIGHT@

#include <harp_internal.hpp>

#include <boost/random.hpp>



using namespace std;
using namespace harp;

static const char * format_sandbox_sim = "sandbox_sim";

static const char * sandbox_sim_image_key_psf = "psf";
static const char * sandbox_sim_image_key_spec = "spec";


harp::image_sandbox_sim::image_sandbox_sim ( boost::property_tree::ptree const & props ) : image ( props ) {

  psf_props_ = props.get_child ( sandbox_sim_image_key_psf );
  spec_props_ = props.get_child ( sandbox_sim_image_key_spec );

  spec_p child_spec ( spec::create ( spec_props_ ) );

  size_t spec_nspec = child_spec->nspec();
  size_t spec_nlambda = child_spec->nlambda();

  psf_p child_psf ( psf::create ( psf_props_ ) );

  size_t psf_nspec = child_psf->nspec();
  size_t psf_nlambda = child_psf->nlambda();

  if ( ( spec_nspec != psf_nspec ) || ( spec_nlambda != psf_nlambda ) ) {
    ostringstream o;
    o << "input spec size (" << spec_nspec << " x " << spec_nlambda << ") is not consistent with input PSF spectral size (" << psf_nspec << " x " << psf_nlambda << ")";
    HARP_THROW( o.str().c_str() );
  }

  rows_ = child_psf->pixrows();
  cols_ = child_psf->pixcols();

  size_t npix = rows_ * cols_;

  size_t nglobal = spec_nspec * spec_nlambda;

  // read spectra

  matrix_dist specdata ( nglobal, 1 );
  vector < double > spec_lambda;

  child_spec->read ( specdata, spec_lambda, sky_ );

  // check wavelength solution against the one from the PSF

  vector < double > psf_lambda = child_psf->lambda();

  for ( size_t i = 0; i < psf_nlambda; ++i ) {
    if ( fabs ( psf_lambda[i] - spec_lambda[i] ) / psf_lambda[i] > 0.001 ) {
      ostringstream o;
      o << "wavelength point " << i << " does not match between PSF (" << psf_lambda[i] << ") and spec (" << spec_lambda[i] << ")";
      HARP_THROW( o.str().c_str() );
    }
  }

  // get design matrix from PSF and project

  matrix_sparse design;

  measured_.ResizeTo ( npix, 1 );
  local_matrix_zero ( measured_ );

  child_psf->projection ( 0, psf_nspec - 1, 0, psf_nlambda - 1, design );

  spec_project ( design, specdata, measured_ );

  // add noise to get measured image

  invcov_.ResizeTo ( npix, 1 );
  local_matrix_zero ( invcov_ );

  typedef boost::ecuyer1988 base_generator_type;
  base_generator_type generator(42u);

  double rms;

  for ( size_t i = 0; i < npix; ++i ) {
    rms = sqrt( 16.0 + measured_.Get ( i, 0 ) );
    invcov_.Set ( i, 0, 1.0 / (rms*rms) );

    boost::normal_distribution < double > dist ( 0.0, rms );
    boost::variate_generator < base_generator_type&, boost::normal_distribution < double > > gauss ( generator, dist );

    measured_.Set ( i, 0, measured_.Get(i,0) + gauss() );
  }
  
}


harp::image_sandbox_sim::~image_sandbox_sim ( ) {
  
  cleanup();
  
}


boost::property_tree::ptree harp::image_sandbox_sim::serialize ( ) {
  boost::property_tree::ptree ret;

  ret.put ( "format", image::format() );

  return ret;
}


void harp::image_sandbox_sim::read ( matrix_local & data ) {

  data = measured_;

  return;
}


void harp::image_sandbox_sim::write ( std::string const & path, matrix_local & data ) {

  HARP_THROW( "sandbox_sim image does not support writing" );
  
  return;
}


void harp::image_sandbox_sim::read_noise ( matrix_local & data ) {

  data = invcov_;

  return;
}


void harp::image_sandbox_sim::write_noise ( std::string const & path, matrix_local & data ) {
  
  HARP_THROW( "sandbox_sim image does not support writing" );

  return;
}



