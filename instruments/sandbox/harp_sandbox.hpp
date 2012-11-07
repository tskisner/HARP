// @COPYRIGHT@

namespace harp {
  
  class image_sandbox : public image {
    
    public :
      image_sandbox ( boost::property_tree::ptree const & props );
      ~image_sandbox ( );

      boost::property_tree::ptree serialize ( );

      size_t rows ( ) { return rows_; }
      size_t cols ( ) { return cols_; }
      
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
    
  };
  
  
  class spec_sandbox : public spec {
    
    public :
      spec_sandbox ( boost::property_tree::ptree const & props );
      ~spec_sandbox ( );

      boost::property_tree::ptree serialize ( );
      
      size_t nspectrum ( ) { return nspec_; }
      size_t spectrum_size ( size_t spectrum ) { return specsize_; }
      void read ( matrix_dist & data );
      void write ( std::string const & path, matrix_dist & data );
    
    private :
    
      size_t size_;
      size_t nspec_;
      size_t specsize_;
      std::string path_;
      int hdu_;
    
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
      size_t specsize ( ) { return specsize_; }
      size_t pixrows ( ) { return rows_; }
      size_t pixcols ( ) { return cols_; }
      void projection ( size_t firstbin, size_t lastbin, matrix_sparse & data );

      void fake_spec ( matrix_dist & data );
      
    private :

      int hdu_info ( fitsfile *fp, const char * sandbox_psf_hdu );
    
      void cache ( size_t firstlocal, size_t nlocal );

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
      double fake_background_;
      double fake_peak_amp_;
      size_t fake_peak_space_;
      size_t fake_peak_obj_;
      double fake_psf_fwhm_;

      std::string path_;
      size_t nspec_;
      size_t specsize_;
      size_t rows_;
      size_t cols_;
      size_t nglobal_;
      size_t npix_;
      size_t pixcorr_;
      std::map < std::string, int > hdus_;
      std::map < size_t, psf_sandbox_resp > resp_;
      
  };
  
  
}