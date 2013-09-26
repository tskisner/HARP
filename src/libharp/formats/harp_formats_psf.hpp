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
          ar & x;
          ar & y;
          ar & lambda;
          ar & amp;
          ar & maj;
          ar & min;
          ar & ang;
      }

  };
  
  // common function to project a



  // elliptical gaussian class read from a file

  class psf_gauss : public psf {

    friend class boost::serialization::access;
    
    public :
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
      void save ( Archive & ar, const unsigned int version ) const {
          ar << boost::serialization::base_object < psf > (*this);
          ar << path_;
          ar << nspec_;
          ar << nlambda_;
          ar << rows_;
          ar << cols_;
          ar << nglobal_;
          ar << npix_;
          ar << pixcorr_;
          ar << lambda_;
          ar << hdus_;
          ar << resp_;
      }
      template < class Archive >
      void load ( Archive & ar, const unsigned int version ) {
          ar >> boost::serialization::base_object < psf > (*this);
          ar >> path_;
          ar >> nspec_;
          ar >> nlambda_;
          ar >> rows_;
          ar >> cols_;
          ar >> nglobal_;
          ar >> npix_;
          ar >> pixcorr_;
          ar >> lambda_;
          ar >> hdus_;
          ar >> resp_;
      }
      BOOST_SERIALIZATION_SPLIT_MEMBER()

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
  

  // simulated elliptical gaussian PSF

  class psf_gauss_sim : public psf {

    friend class boost::serialization::access;
    
    public :
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
      void save ( Archive & ar, const unsigned int version ) const {
          ar << boost::serialization::base_object < psf > (*this);
          ar << nspec_;
          ar << nlambda_;
          ar << rows_;
          ar << cols_;
          ar << nglobal_;
          ar << npix_;
          ar << pixcorr_;
          ar << lambda_;
          ar << hdus_;
          ar << resp_;
          ar << n_bundle_;
          ar << bundle_size_;
          ar << pix_margin_;
          ar << pix_gap_;
          ar << pix_bundle_;
          ar << pix_offset_;
          ar << psf_fwhm_;
          ar << first_lambda_;
          ar << last_lambda_;
          ar << lambda_spec_props_;
      }
      template < class Archive >
      void load ( Archive & ar, const unsigned int version ) {
          ar >> boost::serialization::base_object < psf > (*this);
          ar >> nspec_;
          ar >> nlambda_;
          ar >> rows_;
          ar >> cols_;
          ar >> nglobal_;
          ar >> npix_;
          ar >> pixcorr_;
          ar >> lambda_;
          ar >> hdus_;
          ar >> resp_;
          ar >> n_bundle_;
          ar >> bundle_size_;
          ar >> pix_margin_;
          ar >> pix_gap_;
          ar >> pix_bundle_;
          ar >> pix_offset_;
          ar >> psf_fwhm_;
          ar >> first_lambda_;
          ar >> last_lambda_;
          ar >> lambda_spec_props_;
      }
      BOOST_SERIALIZATION_SPLIT_MEMBER()

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
  

}

