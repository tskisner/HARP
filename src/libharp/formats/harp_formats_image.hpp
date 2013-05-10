// @COPYRIGHT@

namespace harp {
  
  class image_fits : public image {
    
    public :
      image_fits ( boost::property_tree::ptree const & props );
      ~image_fits ( );

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


  class image_sim : public image {
    
    public :
      image_sim ( boost::property_tree::ptree const & props );
      ~image_sim ( );

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
    
  };
  
  
}

