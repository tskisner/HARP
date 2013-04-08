// @COPYRIGHT@

namespace harp {
  
  class psf_gauss_resp {

    public :

      psf_gauss_resp ( ) { }

      ~psf_gauss_resp ( ) {
        x.Empty();
        y.Empty();
        lambda.Empty();
        amp.Empty();
        maj.Empty();
        min.Empty();
        ang.Empty();
      }

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

      matrix_local x;
      matrix_local y;
      matrix_local lambda;
      matrix_local amp;
      matrix_local maj;
      matrix_local min;
      matrix_local ang;
  };

  typedef boost::shared_ptr < psf_gauss_resp > psf_gauss_resp_p;
  

  class psf_gauss : public psf {
    
    public :
      psf_gauss ( boost::property_tree::ptree const & props );
      ~psf_gauss ( );

      boost::property_tree::ptree serialize ( );

      size_t nspec ( ) { return nspec_; }
      size_t nlambda ( ) { return nlambda_; }
      size_t pixrows ( ) { return rows_; }
      size_t pixcols ( ) { return cols_; }
      std::vector < double > lambda ( ) { return lambda_; }
      void projection ( size_t first_spec, size_t last_spec, size_t first_lambda, size_t last_lambda, matrix_sparse & data );
      
    private :

      int hdu_info ( fitsfile *fp, const char * gauss_psf_hdu );
    
      void cache ( size_t first_spec, size_t last_spec );

      void extent ( size_t firstspec, size_t lastspec, size_t firstbin, size_t lastbin, size_t & firstcol, size_t & firstrow, size_t & lastcol, size_t & lastrow );
      
      void gauss_sample ( matrix_local & vals, matrix_local & xrel, matrix_local & yrel, double amp, double maj, double min, double ang );

      void fake_spec2pix ( size_t spec, size_t bin, size_t & row, size_t & col );
      
      bool dofake_;
      bool cached_;
      size_t cached_first_;
      size_t cached_last_;
      size_t fake_n_bundle_;
      size_t fake_bundle_size_;
      size_t fake_pix_margin_;
      size_t fake_pix_gap_;
      size_t fake_pix_bundle_;
      double fake_psf_fwhm_;
      double fake_first_lambda_;
      double fake_last_lambda_;
      boost::property_tree::ptree fake_spec_props_;

      std::string path_;
      size_t nspec_;
      size_t nlambda_;
      size_t rows_;
      size_t cols_;
      size_t nglobal_;
      size_t npix_;
      size_t pixcorr_;
      std::vector < double > lambda_;
      std::vector < bool > sky_;
      std::map < std::string, int > hdus_;
      std::map < size_t, psf_gauss_resp_p > resp_;
      
  };
  
  
}

