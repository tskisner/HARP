// @COPYRIGHT@

namespace harp {
  
  class image_fits : public image {

    friend class boost::serialization::access;
    
    public :
      image_fits ( boost::property_tree::ptree const & props );
      ~image_fits ( );
      
      size_t n_rows ( ) { return rows_; }
      
      size_t n_cols ( ) { return cols_; }
      
      void read ( vector_double & data, vector_double & invvar, std::vector < bool > & sky );

      void write ( std::string const & path, vector_double & data, vector_double & invvar, std::vector < bool > & sky );
    
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

    friend class boost::serialization::access;
    
    public :
      image_sim ( boost::property_tree::ptree const & props );
      ~image_sim ( );

      size_t n_rows ( ) { return rows_; }
      
      size_t n_cols ( ) { return cols_; }
      
      void read ( vector_double & data, vector_double & invvar, std::vector < bool > & sky );

      void write ( std::string const & path, vector_double & data, vector_double & invvar, std::vector < bool > & sky );
    
    private :
    
      boost::property_tree::ptree spec_props_;
      boost::property_tree::ptree psf_props_;
      size_t rows_;
      size_t cols_;
      vector_double measured_;
      vector_double invcov_;
      std::vector < bool > sky_;
    
  };
  
  
}

