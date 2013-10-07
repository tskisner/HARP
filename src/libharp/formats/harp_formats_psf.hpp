// @COPYRIGHT@

namespace harp {
  

  // class for sampling from an elliptical gaussian

  class psf_gauss_resp {

    friend class boost::serialization::access;

    public :

      psf_gauss_resp ( ) { }

      ~psf_gauss_resp ( ) { }

      psf_gauss_resp ( psf_gauss_resp const & orig ) :
        x( orig.x ),
        y( orig.y ),
        lambda( orig.lambda ),
        amp( orig.amp ),
        maj( orig.maj ),
        min( orig.min ),
        ang( orig.ang ) { }

      psf_gauss_resp & operator= ( psf_gauss_resp const & rhs ) {
        if ( &rhs != this ) {
          x = rhs.x;
          y = rhs.y;
          lambda = rhs.lambda;
          amp = rhs.amp;
          maj = rhs.maj;
          min = rhs.min;
          ang = rhs.ang;
        }
        return *this;
      }

      void sample ( size_t x_offset, size_t y_offset, matrix_double & patch );

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

  //BOOST_CLASS_EXPORT_GUID(psf_gauss_resp, "psf_gauss_resp")


  // elliptical gaussian class read from a file

  class psf_gauss : public psf {

    friend class boost::serialization::access;
    
    public :

      psf_gauss ( ) : psf () {
        nspec_ = 0;
        nlambda_ = 0;
        rows_ = 0;
        cols_ = 0;
        nglobal_ = 0;
        npix_ = 0;
        pixcorr_ = 0;
        path_ = "";
      }

      psf_gauss ( boost::property_tree::ptree const & props );

      ~psf_gauss ( );

      size_t n_spec ( ) { return nspec_; }

      size_t n_lambda ( ) { return nlambda_; }
      
      size_t img_rows ( ) { return rows_; }
      
      size_t img_cols ( ) { return cols_; }
      
      vector_double lambda ( ) { return lambda_; }

      void write ( std::string const & path );

      void response ( size_t spec, size_t lambda, size_t & x_offset, size_t & y_offset, matrix_double & patch );

      // This provides a public interface to the gaussian parameters for each flux bin

      psf_gauss_resp & parameters ( size_t spec, size_t lambda ) {
        size_t bin = spec * nlambda_ + lambda;
        return resp_[ bin ];
      }

      psf_gauss_resp & parameters ( size_t bin ) {
        return resp_[ bin ];
      }

    private :

      template < class Archive >
      void serialize ( Archive & ar, const unsigned int version ) const {
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

  //BOOST_CLASS_EXPORT_GUID(psf_gauss, "psf_gauss")


  // simulated elliptical gaussian PSF

  class psf_gauss_sim : public psf {

    friend class boost::serialization::access;
    
    public :

      psf_gauss_sim ( ) : psf () {
        nspec_ = 0;
        nlambda_ = 0;
        rows_ = 0;
        cols_ = 0;
        nglobal_ = 0;
        npix_ = 0;
        pixcorr_ = 0;
        n_bundle_ = 0;
        bundle_size_ = 0;
        pix_margin_ = 0.0;
        pix_gap_ = 0.0;
        pix_bundle_ = 0.0;
        pix_offset_ = 0.0;
        psf_fwhm_ = 0.0;
        first_lambda_ = 0.0;
        last_lambda_ = 0.0;
      }

      psf_gauss_sim ( boost::property_tree::ptree const & props );

      ~psf_gauss_sim ( );

      size_t n_spec ( ) { return nspec_; }

      size_t n_lambda ( ) { return nlambda_; }
      
      size_t img_rows ( ) { return rows_; }
      
      size_t img_cols ( ) { return cols_; }
      
      vector_double lambda ( ) { return lambda_; }

      void write ( std::string const & path );

      void response ( size_t spec, size_t lambda, size_t & x_offset, size_t & y_offset, matrix_double & patch );

      // This provides a public interface to the gaussian parameters for each flux bin

      psf_gauss_resp & parameters ( size_t spec, size_t lambda ) {
        size_t bin = spec * nlambda_ + lambda;
        return resp_[ bin ];
      }

      psf_gauss_resp & parameters ( size_t bin ) {
        return resp_[ bin ];
      }

    private :

      template < class Archive >
      void serialize ( Archive & ar, const unsigned int version ) const {
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

  //BOOST_CLASS_EXPORT_GUID(psf_gauss_sim, "psf_gauss_sim")


}

