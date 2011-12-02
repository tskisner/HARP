// @COPYRIGHT@


namespace harp {
  
  class image_boss : public image {
    
    public :
      image_boss ( std::map < std::string, std::string > const & params );
      ~image_boss ( );
      
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
  
  
  class spec_boss : public spec {
    
    public :
      spec_boss ( std::map < std::string, std::string > const & params );
      ~spec_boss ( );
      
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
  } psf_boss_resp;
  
  typedef struct {
    size_t igroup;
    double x0;
    double xscale;
    double y0;
    double yscale;
  } psf_boss_pixscale;
  
  class psf_boss : public psf {
    
    public :
      psf_boss ( std::map < std::string, std::string > const & params );
      ~psf_boss ( );
      
      size_t nspec ( ) { return nspec_; }
      size_t specsize ( size_t specnum ) { return nbins_; }
      void lambda ( size_t specnum, vec_dense & data );
      void extent ( size_t firstspec, size_t lastspec, size_t firstbin, size_t lastbin, size_t & firstX, size_t & firstY, size_t & lastX, size_t & lastY );
      void projection ( std::string profcalc, std::string profremap, size_t firstspec, size_t lastspec, size_t firstbin, size_t lastbin, size_t firstX, size_t lastX, size_t firstY, size_t lastY, mat_sparserow & data );
      
            
    private :
    
      void cache_spec ( size_t first, size_t last );
      
      void gauss_sample ( vec_dense & vals, vec_dense & xrel, vec_dense & yrel, double amp, double maj, double min, double ang );
      
      void shift_kernel ( vec_dense & vals, double delta, vec_dense & data );
      
      void sshift ( double dx, double dy, size_t xoff, size_t yoff, size_t xsize, size_t ysize );
      
      size_t valid_range ( size_t const & firstX, size_t const & lastX, size_t const & firstY, size_t const & lastY, size_t & startX, size_t & stopX, size_t & startY, size_t & stopY, size_t & spec, size_t & bin );
    
      std::string path_;
      size_t nspec_;
      size_t nbins_;
      size_t xpix_;
      size_t ypix_;
      size_t xpixfft_;
      size_t ypixfft_;
      size_t xpixcorr_;
      size_t ypixcorr_;
      size_t xpixwidth_;
      size_t ypixwidth_;
      size_t ncoeff_;
      size_t ngroup_;
      std::map < std::string, int > hdus_;
      std::map < size_t, psf_boss_resp > resp_;
      std::map < size_t, psf_boss_pixscale > xyscale_;
      vec_dense xexp_;
      vec_dense yexp_;
      std::map < size_t, std::map < size_t, vec_dense > > psfimage_;
      std::string type_;
      fftw_plan fsincplan_;
      fftw_plan fpixplan_;
      fftw_plan rpixplan_;
      double * infft_;
      fftw_complex * outfft_;
      double * insinc_;
      fftw_complex * outsinc_;
      
  };
  
  
}