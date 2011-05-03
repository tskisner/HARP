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
      std::string name_;
    
  };
  
  
}