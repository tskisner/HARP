// @COPYRIGHT@

namespace harp {
  
  static const char * format_toy = "toy";
  
  class image_toy : public image {
    
    public :
      image_toy ( std::map < std::string, std::string > const & params );
      ~image_toy ( );
      
      size_t rows ( ) { return rows_; }
      size_t cols ( ) { return cols_; }
      void rowscols ( size_t & rows, size_t & cols ) { rows = rows_; cols = cols_; return; }
      void read ( size_t startrow, size_t startcol, harp::dense_mat_view & data );
      void write ( std::string const & path, size_t startrow, size_t startcol, harp::dense_mat_view & data ) { }
    
    private :
    
      size_t rows_;
      size_t cols_;
      std::string path_;
      int hdu_;
    
  };
  
  
  class spectra_toy : public spectra {
    
    public :
      spectra_toy ( std::map < std::string, std::string > const & params );
      ~spectra_toy ( );
      
      size_t size ( ) { return size_; }
      void read ( harp::data_vec & data );
      void write ( std::string const & path, harp::data_vec & data ) { }
    
    private :
    
      size_t size_;
      size_t pos_;
      std::string path_;
      int hdu_;
    
  };
  
  
}