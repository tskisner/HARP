// @COPYRIGHT@

namespace harp {
  
  class image_boss : public image {
    
    public :
      image_boss ( std::map < std::string, std::string > const & params );
      ~image_boss ( );
      
      size_t rows ( ) { return rows_; }
      size_t cols ( ) { return cols_; }
      void rowscols ( size_t & rows, size_t & cols ) { rows = rows_; cols = cols_; return; }
      void read ( size_t startrow, size_t startcol, harp::dense_rowmat_view & data );
      void write ( std::string const & path, size_t startrow, size_t startcol, harp::dense_rowmat_view & data );
    
    private :
    
      size_t rows_;
      size_t cols_;
      std::string path_;
      int hdu_;
    
  };
  
  
  class spectrum_boss : public spectrum {
    
    public :
      spectrum_boss ( std::map < std::string, std::string > const & params );
      ~spectrum_boss ( );
      
      size_t size ( ) { return size_; }
      void read ( harp::data_vec_view & data );
      void write ( std::string const & path, harp::data_vec_view & data );
    
    private :
    
      size_t size_;
      size_t pos_;
      std::string path_;
      int hdu_;
    
  };
  
  
  typedef struct {
    data_vec x;
    data_vec y;
  } psf_boss_resp;
  
  class psf_boss : public psf {
    
    public :
      psf_boss ( std::map < std::string, std::string > const & params );
      ~psf_boss ( );
      
      size_t nspec ( ) { return nspec_; }
      size_t specsize ( size_t specnum ) { return nbins_; }
      void lambda ( size_t specnum, data_vec & data );
      void extent ( size_t firstspec, size_t lastspec, size_t firstbin, size_t lastbin, size_t & firstX, size_t & firstY, size_t & lastX, size_t & lastY );
      void projection ( size_t firstspec, size_t lastspec, size_t firstbin, size_t lastbin, size_t firstX, size_t lastX, size_t firstY, size_t lastY, comp_rowmat & data );
      
            
    private :
    
      void cache_spec ( size_t first, size_t last );
      
      size_t valid_range ( size_t const & firstX, size_t const & lastX, size_t const & firstY, size_t const & lastY, size_t & startX, size_t & stopX, size_t & startY, size_t & stopY, size_t & spec, size_t & bin );
    
      std::string path_;
      size_t nspec_;
      size_t nbins_;
      size_t pixcorr_;
      std::map < std::string, int > hdus_;
      std::map < size_t, psf_boss_resp > resp_;
      
  };
  
  
}