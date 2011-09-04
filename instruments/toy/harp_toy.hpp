// @COPYRIGHT@

namespace harp {
  
  class image_toy : public image {
    
    public :
      image_toy ( std::map < std::string, std::string > const & params );
      ~image_toy ( );
      
      size_t rows ( ) { return rows_; }
      size_t cols ( ) { return cols_; }
      void rowscols ( size_t & rows, size_t & cols ) { rows = rows_; cols = cols_; return; }
      
      void read ( size_t startrow, size_t startcol, mat_denserow & data );
      void write ( std::string const & path, size_t startrow, size_t startcol, mat_denserow & data );
      void read ( vec_dense & data );
      void write ( std::string const & path, vec_dense & data );
      
      void read_noise ( size_t startrow, size_t startcol, mat_denserow & data );
      void write_noise ( std::string const & path, size_t startrow, size_t startcol, mat_denserow & data );
      void read_noise ( vec_dense & data );
      void write_noise ( std::string const & path, vec_dense & data );
    
    private :
    
      size_t rows_;
      size_t cols_;
      std::string path_;
      int sighdu_;
      int nsehdu_;
    
  };
  
  
  class spec_toy : public spec {
    
    public :
      spec_toy ( std::map < std::string, std::string > const & params );
      ~spec_toy ( );
      
      size_t size ( ) { return size_; }
      size_t nspectrum ( ) { return nspec_; }
      size_t spectrum_size ( size_t spectrum ) { return specsize_; }
      void read ( vec_dense & data );
      void write ( std::string const & path, vec_dense & data );
      void read_spectrum ( size_t spectrum, vec_dense & data );
      void write_spectrum ( size_t spectrum, std::string const & path, vec_dense & data );
    
    private :
    
      size_t size_;
      size_t nspec_;
      size_t specsize_;
      std::string path_;
      int hdu_;
    
  };
  
  
  typedef struct {
    vec_dense x;
    vec_dense y;
    vec_dense lambda;
    vec_dense amp;
    vec_dense maj;
    vec_dense min;
    vec_dense ang;
  } psf_toy_resp;
  
  class psf_toy : public psf {
    
    public :
      psf_toy ( std::map < std::string, std::string > const & params );
      ~psf_toy ( );
      
      size_t nspec ( ) { return nspec_; }
      size_t specsize ( size_t specnum ) { return nreduced_; }
      void lambda ( size_t specnum, vec_dense & data );
      void extent ( size_t firstspec, size_t lastspec, size_t firstbin, size_t lastbin, size_t & firstX, size_t & firstY, size_t & lastX, size_t & lastY );
      void projection ( std::string profcalc, std::string profremap, size_t firstspec, size_t lastspec, size_t firstbin, size_t lastbin, size_t firstX, size_t lastX, size_t firstY, size_t lastY, mat_comprow & data );
      
            
    private :
    
      void cache_spec ( size_t first, size_t last );
      
      size_t valid_range ( size_t const & firstX, size_t const & lastX, size_t const & firstY, size_t const & lastY, size_t & startX, size_t & stopX, size_t & startY, size_t & stopY, size_t & spec, size_t & bin );
      
      void gauss_sample ( vec_dense & vals, vec_dense & xrel, vec_dense & yrel, double amp, double maj, double min, double ang );
      
      void gauss_sample_alt ( vec_dense & vals, vec_dense & xrel, vec_dense & yrel, double amp, double maj, double min, double ang );
    
      std::string path_;
      size_t binning_;
      size_t nspec_;
      size_t nbins_;
      size_t nreduced_;
      size_t pixcorr_;
      std::map < std::string, int > hdus_;
      std::map < size_t, psf_toy_resp > resp_;
      
  };
  
  
}