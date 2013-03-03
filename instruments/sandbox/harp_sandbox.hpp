// @COPYRIGHT@

namespace harp {
  
  class image_sandbox : public image {
    
    public :
      image_sandbox ( boost::property_tree::ptree const & props );
      ~image_sandbox ( );

      boost::property_tree::ptree serialize ( );

      size_t rows ( ) { return rows_; }
      size_t cols ( ) { return cols_; }

      std::vector < bool > sky ( ) { return sky_; }
      
      void read ( matrix_local & data );
      void write ( std::string const & path, matrix_local & data );
      
      void read_noise ( matrix_local & data );
      void write_noise ( std::string const & path, matrix_local & data );
    
    private :
    
      size_t rows_;
      size_t cols_;
      std::string path_;
      int sighdu_;
      int nsehdu_;
      int skyhdu_;
      std::vector < bool > sky_;
    
  };


  class image_sandbox_sim : public image {
    
    public :
      image_sandbox_sim ( boost::property_tree::ptree const & props );
      ~image_sandbox_sim ( );

      boost::property_tree::ptree serialize ( );

      size_t rows ( ) { return rows_; }
      size_t cols ( ) { return cols_; }

      std::vector < bool > sky ( ) { return sky_; }
      
      void read ( matrix_local & data );
      void write ( std::string const & path, matrix_local & data );
      
      void read_noise ( matrix_local & data );
      void write_noise ( std::string const & path, matrix_local & data );
    
    private :
    
      boost::property_tree::ptree spec_props_;
      boost::property_tree::ptree psf_props_;
      size_t rows_;
      size_t cols_;
      matrix_local measured_;
      matrix_local invcov_;
      std::vector < bool > sky_;
      bool debug_;
    
  };
  
  
  class spec_sandbox : public spec {
    
    public :
      spec_sandbox ( boost::property_tree::ptree const & props );
      ~spec_sandbox ( );

      boost::property_tree::ptree serialize ( );
      
      size_t nspec ( ) { return nspec_; }
      size_t nlambda ( ) { return nlambda_; }

      void read ( matrix_dist & data, std::vector < double > & lambda, std::vector < bool > & sky );
      void write ( std::string const & path, matrix_dist & data, std::vector < double > const & lambda, std::vector < bool > const & sky );
    
    private :
    
      size_t size_;
      size_t nspec_;
      size_t nlambda_;
      double background_;
      double atmpeak_;
      double objpeak_;
      size_t atmspace_;
      size_t objspace_;
      size_t skymod_;
      double first_lambda_;
      double last_lambda_;
    
  };
  
  
  typedef struct {
    matrix_local x;
    matrix_local y;
    matrix_local lambda;
    matrix_local amp;
    matrix_local maj;
    matrix_local min;
    matrix_local ang;
  } psf_sandbox_resp;
  
  class psf_sandbox : public psf {
    
    public :
      psf_sandbox ( boost::property_tree::ptree const & props );
      ~psf_sandbox ( );

      boost::property_tree::ptree serialize ( );

      size_t nspec ( ) { return nspec_; }
      size_t nlambda ( ) { return nlambda_; }
      size_t pixrows ( ) { return rows_; }
      size_t pixcols ( ) { return cols_; }
      virtual std::vector < double > lambda ( ) { return lambda_; }
      void projection ( size_t first_spec, size_t last_spec, size_t first_lambda, size_t last_lambda, matrix_sparse & data );
      
    private :

      int hdu_info ( fitsfile *fp, const char * sandbox_psf_hdu );
    
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

      std::string path_;
      size_t nspec_;
      size_t nlambda_;
      size_t rows_;
      size_t cols_;
      size_t nglobal_;
      size_t npix_;
      size_t pixcorr_;
      std::vector < double > lambda_;
      std::map < std::string, int > hdus_;
      std::map < size_t, psf_sandbox_resp > resp_;
      
  };
  
  
}