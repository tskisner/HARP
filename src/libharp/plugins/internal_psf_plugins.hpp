/*
  High Performance Astrophysical Reconstruction and Processing (HARP)
  (c) 2014-2015, The Regents of the University of California, 
  through Lawrence Berkeley National Laboratory.  See top
  level LICENSE file for details.
*/

#ifndef HARP_INTERNAL_PSF_PLUGINS_HPP
#define HARP_INTERNAL_PSF_PLUGINS_HPP

namespace harp {
  

  // class for sampling from an elliptical gaussian

  class psf_gauss_resp {

    friend class boost::serialization::access;

    public :

      psf_gauss_resp ( ) :
        x( 0.0 ),
        y( 0.0 ),
        lambda( 0.0 ),
        amp( 0.0 ),
        maj( 0.0 ),
        min( 0.0 ),
        ang( 0.0 ) { }

      ~psf_gauss_resp ( ) { }

      psf_gauss_resp ( psf_gauss_resp const & orig ) :
        x( orig.x ),
        y( orig.y ),
        lambda( orig.lambda ),
        amp( orig.amp ),
        maj( orig.maj ),
        min( orig.min ),
        ang( orig.ang ) { }

      psf_gauss_resp & operator= ( psf_gauss_resp const & rhs );

      void sample ( size_t x_offset, size_t y_offset, matrix_double & patch ) const;

      double x;
      double y;
      double lambda;
      double amp;
      double maj;
      double min;
      double ang;

    private :

      template < class Archive >
      void serialize ( Archive & ar, const unsigned int version ) {
        ar & BOOST_SERIALIZATION_NVP(x);
        ar & BOOST_SERIALIZATION_NVP(y);
        ar & BOOST_SERIALIZATION_NVP(lambda);
        ar & BOOST_SERIALIZATION_NVP(amp);
        ar & BOOST_SERIALIZATION_NVP(maj);
        ar & BOOST_SERIALIZATION_NVP(min);
        ar & BOOST_SERIALIZATION_NVP(ang);
        return;
      }

  };

  BOOST_SERIALIZATION_SHARED_PTR(psf_gauss_resp)


  // elliptical gaussian class read from a file

  class psf_gauss : public psf {

    friend class boost::serialization::access;
    
    public :

      psf_gauss ( );

      psf_gauss ( boost::property_tree::ptree const & props );

      ~psf_gauss ( );

      void write ( std::string const & path );

      // overloaded virtual methods from base class

      size_t n_spec ( ) const;

      size_t n_lambda ( ) const;
      
      size_t img_rows ( ) const;
      
      size_t img_cols ( ) const;
      
      vector_double lambda ( ) const;

      void extent ( size_t spec, size_t lambda, size_t & x_offset, size_t & y_offset, size_t & n_x, size_t & n_y ) const;

      void response ( size_t spec, size_t lambda, size_t & x_offset, size_t & y_offset, matrix_double & patch ) const;

      size_t response_nnz_estimate ( ) const;

    private :

      template < class Archive >
      void serialize ( Archive & ar, const unsigned int version ) {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(psf);
        ar & BOOST_SERIALIZATION_NVP(path_);
        ar & BOOST_SERIALIZATION_NVP(nspec_);
        ar & BOOST_SERIALIZATION_NVP(nlambda_);
        ar & BOOST_SERIALIZATION_NVP(rows_);
        ar & BOOST_SERIALIZATION_NVP(cols_);
        ar & BOOST_SERIALIZATION_NVP(nglobal_);
        ar & BOOST_SERIALIZATION_NVP(npix_);
        ar & BOOST_SERIALIZATION_NVP(pixcorr_);
        ar & BOOST_SERIALIZATION_NVP(lambda_);
        ar & BOOST_SERIALIZATION_NVP(hdus_);
        ar & BOOST_SERIALIZATION_NVP(resp_);
        return;
      }

      int hdu_info ( fitsfile *fp, const char * gauss_psf_hdu );

      std::string path_;
      size_t nspec_;
      size_t nlambda_;
      size_t rows_;
      size_t cols_;
      size_t nglobal_;
      size_t npix_;
      size_t pixcorr_;
      vector_double lambda_;
      std::map < std::string, int > hdus_;
      std::vector < psf_gauss_resp > resp_;
      
  };

  BOOST_SERIALIZATION_SHARED_PTR(psf_gauss)

  psf * psf_gauss_create ( boost::property_tree::ptree const & props );


  // simulated elliptical gaussian PSF

  class psf_gauss_sim : public psf {

    friend class boost::serialization::access;
    
    public :

      psf_gauss_sim ( );

      psf_gauss_sim ( boost::property_tree::ptree const & props );

      ~psf_gauss_sim ( );

      void write ( std::string const & path );

      // overloaded virtual methods from base class

      size_t n_spec ( ) const;

      size_t n_lambda ( ) const;
      
      size_t img_rows ( ) const;
      
      size_t img_cols ( ) const;
      
      vector_double lambda ( ) const;

      void extent ( size_t spec, size_t lambda, size_t & x_offset, size_t & y_offset, size_t & n_x, size_t & n_y ) const;

      void response ( size_t spec, size_t lambda, size_t & x_offset, size_t & y_offset, matrix_double & patch ) const;

      size_t response_nnz_estimate ( ) const;

    private :

      template < class Archive >
      void serialize ( Archive & ar, const unsigned int version ) {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(psf);
        ar & BOOST_SERIALIZATION_NVP(nspec_);
        ar & BOOST_SERIALIZATION_NVP(nlambda_);
        ar & BOOST_SERIALIZATION_NVP(rows_);
        ar & BOOST_SERIALIZATION_NVP(cols_);
        ar & BOOST_SERIALIZATION_NVP(nglobal_);
        ar & BOOST_SERIALIZATION_NVP(npix_);
        ar & BOOST_SERIALIZATION_NVP(pixcorr_);
        ar & BOOST_SERIALIZATION_NVP(lambda_);
        ar & BOOST_SERIALIZATION_NVP(hdus_);
        ar & BOOST_SERIALIZATION_NVP(resp_);
        ar & BOOST_SERIALIZATION_NVP(n_bundle_);
        ar & BOOST_SERIALIZATION_NVP(bundle_size_);
        ar & BOOST_SERIALIZATION_NVP(pix_margin_);
        ar & BOOST_SERIALIZATION_NVP(pix_gap_);
        ar & BOOST_SERIALIZATION_NVP(pix_bundle_);
        ar & BOOST_SERIALIZATION_NVP(pix_offset_);
        ar & BOOST_SERIALIZATION_NVP(psf_fwhm_);
        ar & BOOST_SERIALIZATION_NVP(first_lambda_);
        ar & BOOST_SERIALIZATION_NVP(last_lambda_);
        ar & BOOST_SERIALIZATION_NVP(lambda_spec_type_);
        ar & BOOST_SERIALIZATION_NVP(lambda_spec_props_);
        return;
      }

      void spec2pix ( size_t spec, size_t bin, double & x, double & y );
      
      size_t n_bundle_;
      size_t bundle_size_;
      double pix_margin_;
      double pix_gap_;
      double pix_bundle_;
      double pix_offset_;
      double psf_fwhm_;
      double first_lambda_;
      double last_lambda_;
      boost::property_tree::ptree lambda_spec_props_;
      std::string lambda_spec_type_;

      size_t nspec_;
      size_t nlambda_;
      size_t rows_;
      size_t cols_;
      size_t nglobal_;
      size_t npix_;
      size_t pixcorr_;
      vector_double lambda_;
      std::map < std::string, int > hdus_;
      std::vector < psf_gauss_resp > resp_;
      
  };

  BOOST_SERIALIZATION_SHARED_PTR(psf_gauss_sim)

  psf * psf_gauss_sim_create ( boost::property_tree::ptree const & props );


  // Gauss-Hermite PSF

  class psf_gh : public psf {

    friend class boost::serialization::access;
    
    public :

      psf_gh ( );

      psf_gh ( boost::property_tree::ptree const & props );

      ~psf_gh ( );

      void write ( std::string const & path );

      // Other public methods for updating internal data structures here..


      // overloaded virtual methods from base class

      size_t n_spec ( ) const;

      size_t n_lambda ( ) const;
      
      size_t img_rows ( ) const;
      
      size_t img_cols ( ) const;
      
      vector_double lambda ( ) const;

      void extent ( size_t spec, size_t lambda, size_t & x_offset, size_t & y_offset, size_t & n_x, size_t & n_y ) const;

      void response ( size_t spec, size_t lambda, size_t & x_offset, size_t & y_offset, matrix_double & patch ) const;

      size_t response_nnz_estimate ( ) const;

    private :

      template < class Archive >
      void serialize ( Archive & ar, const unsigned int version ) {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(psf);
        ar & BOOST_SERIALIZATION_NVP(path_);
        ar & BOOST_SERIALIZATION_NVP(nspec_);
        ar & BOOST_SERIALIZATION_NVP(nlambda_);
        ar & BOOST_SERIALIZATION_NVP(rows_);
        ar & BOOST_SERIALIZATION_NVP(cols_);
        ar & BOOST_SERIALIZATION_NVP(nglobal_);
        ar & BOOST_SERIALIZATION_NVP(npix_);
        ar & BOOST_SERIALIZATION_NVP(pixcorr_);
        ar & BOOST_SERIALIZATION_NVP(lambda_);
        return;
      }

      std::string path_;
      size_t nspec_;
      size_t nlambda_;
      size_t rows_;
      size_t cols_;
      size_t nglobal_;
      size_t npix_;
      size_t pixcorr_;
      vector_double lambda_;
      
  };

  BOOST_SERIALIZATION_SHARED_PTR(psf_gh)

  psf * psf_gh_create ( boost::property_tree::ptree const & props );


  // Specter PSF formats

  /*
  class psf_pyspecter : public psf {

    friend class boost::serialization::access;
    
    public :

      psf_pyspecter ( );

      psf_pyspecter ( boost::property_tree::ptree const & props );

      ~psf_pyspecter ( );

      // overloaded virtual methods from base class

      size_t n_spec ( ) const;

      size_t n_lambda ( ) const;
      
      size_t img_rows ( ) const;
      
      size_t img_cols ( ) const;
      
      vector_double lambda ( ) const;

      size_t response_nnz_estimate ( ) const;

      void extent ( size_t spec, size_t lambda, size_t & x_offset, size_t & y_offset, size_t & n_x, size_t & n_y ) const;

      void response ( size_t spec, size_t lambda, size_t & x_offset, size_t & y_offset, matrix_double & patch ) const;

      void extent_multi ( std::map < size_t, std::set < size_t > > const & speclambda, std::vector < size_t > & x_offset, std::vector < size_t > & y_offset, std::vector < size_t > & n_x, std::vector < size_t > & n_y ) const;

      void project ( std::map < size_t, std::set < size_t > > const & speclambda, matrix_double & A ) const;

      void project_transpose ( std::map < size_t, std::set < size_t > > const & speclambda, matrix_double_sparse & AT ) const;

    private :

      template < class Archive >
      void serialize ( Archive & ar, const unsigned int version ) {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(psf);
        ar & BOOST_SERIALIZATION_NVP(path_);
        ar & BOOST_SERIALIZATION_NVP(lambda_min_);
        ar & BOOST_SERIALIZATION_NVP(lambda_max_);
        ar & BOOST_SERIALIZATION_NVP(wavemin_);
        ar & BOOST_SERIALIZATION_NVP(wavemax_);
        ar & BOOST_SERIALIZATION_NVP(wavebin_);
        ar & BOOST_SERIALIZATION_NVP(nspec_);
        ar & BOOST_SERIALIZATION_NVP(nlambda_);
        ar & BOOST_SERIALIZATION_NVP(lambda_);
        ar & BOOST_SERIALIZATION_NVP(nglobal_);
        ar & BOOST_SERIALIZATION_NVP(imgrows_);
        ar & BOOST_SERIALIZATION_NVP(imgcols_);
        ar & BOOST_SERIALIZATION_NVP(npix_);
        ar & BOOST_SERIALIZATION_NVP(nnz_);
        ar & BOOST_SERIALIZATION_NVP(pickled_);
        return;
      }

      std::string path_;
      // maximum range of lambda available
      double lambda_min_;
      double lambda_max_;
      // user-requested wavelengths
      double wavemin_;
      double wavemax_;
      double wavebin_;
      // wavelengths used
      size_t nspec_;
      size_t nlambda_;
      vector_double lambda_;
      size_t nglobal_;
      // image dims
      size_t imgrows_;
      size_t imgcols_;
      size_t npix_;
      // estimate
      size_t nnz_;
      // python class
      std::string pickled_;
      
  };

  BOOST_SERIALIZATION_SHARED_PTR(psf_pyspecter)

  psf * psf_pyspecter_create ( boost::property_tree::ptree const & props );
  */

}

#endif
