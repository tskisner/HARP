// @COPYRIGHT@

#include <harp_internal.hpp>

#include <boost/random.hpp>


using namespace std;
using namespace harp;


static const char * image_sim_key_psf = "psf";
static const char * image_sim_key_spec = "spec";


harp::image_sim::image_sim ( boost::property_tree::ptree const & props ) : image ( props ) {

  psf_props_ = props.get_child ( image_sim_key_psf );
  spec_props_ = props.get_child ( image_sim_key_spec );

  spec_p child_spec ( spec::create ( spec_props_ ) );

  size_t spec_nspec = child_spec->n_spec();
  size_t spec_nlambda = child_spec->n_lambda();

  psf_p child_psf ( psf::create ( psf_props_ ) );

  size_t psf_nspec = child_psf->n_spec();
  size_t psf_nlambda = child_psf->n_lambda();

  if ( ( spec_nspec != psf_nspec ) || ( spec_nlambda != psf_nlambda ) ) {
    ostringstream o;
    o << "input spec size (" << spec_nspec << " x " << spec_nlambda << ") is not consistent with input PSF spectral size (" << psf_nspec << " x " << psf_nlambda << ")";
    HARP_THROW( o.str().c_str() );
  }

  rows_ = child_psf->img_rows();
  cols_ = child_psf->img_cols();

  cerr << "child psf img " << rows_ << " x " << cols_ << endl;

  size_t npix = rows_ * cols_;

  size_t nglobal = spec_nspec * spec_nlambda;

  // read spectra

  vector_double spec_data ( nglobal );
  vector_double spec_lambda ( spec_nlambda );

  child_spec->read ( spec_data, spec_lambda, sky_ );

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

  matrix_double A;

  child_psf->project ( A );

  vector_double img_data ( npix );

  // multiply to get signal data


  /*
  matrix_sparse design;

  matrix_local signal ( npix, 1 );
  local_matrix_zero ( signal );

  child_psf->projection ( 0, psf_nspec - 1, 0, psf_nlambda - 1, design );

  spec_project ( design, specdata, signal );

  fitsfile * fp;

  */

  // add noise to get measured image

  /*

  matrix_local noise ( npix, 1 );
  local_matrix_zero ( noise );

  measured_.ResizeTo ( npix, 1 );
  local_matrix_zero ( measured_ );

  invcov_.ResizeTo ( npix, 1 );
  local_matrix_zero ( invcov_ );

  typedef boost::ecuyer1988 base_generator_type;
  base_generator_type generator(42u);

  double rms;

  for ( size_t i = 0; i < npix; ++i ) {
    rms = sqrt( 16.0 + signal.Get ( i, 0 ) );
    invcov_.Set ( i, 0, 1.0 / (rms*rms) );

    boost::normal_distribution < double > dist ( 0.0, rms );
    boost::variate_generator < base_generator_type&, boost::normal_distribution < double > > gauss ( generator, dist );

    noise.Set ( i, 0, gauss() );

    measured_.Set ( i, 0, signal.Get(i,0) + noise.Get(i,0) );
  }

  */
  
}


harp::image_sim::~image_sim ( ) {
  
}


void harp::image_sim::read ( vector_double & data, vector_double & invvar, std::vector < bool > & sky ) {

  size_t npix = rows_ * cols_;

  cerr << "read npix = " << rows_ << " x " << cols_ << " = " << npix << endl;

  data.resize ( npix );
  invvar.resize ( npix );

  data = signal_ + noise_;

  invvar = invcov_;

  sky = sky_;

  return;
}


void harp::image_sim::write ( std::string const & path, vector_double & data, vector_double & invvar, std::vector < bool > & sky ) {

  fitsfile * fp;
    
  fits::create ( fp, path );
  
  fits::img_append < double > ( fp, rows_, cols_ );
  
  fits::img_write ( fp, data );

  fits::img_append < double > ( fp, rows_, cols_ );
  
  fits::img_write ( fp, invvar );

  specter_write_sky ( fp, sky );

  fits::img_append < double > ( fp, rows_, cols_ );

  fits::write_key ( fp, "EXTNAME", "signal", "simulated signal" );
  fits::img_write ( fp, signal_ );

  fits::img_append < double > ( fp, rows_, cols_ );

  fits::write_key ( fp, "EXTNAME", "noise", "simulated noise" );
  fits::img_write ( fp, noise_ );

  fits::close ( fp );

  return;
}



