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

      template < class Archive >
      void save ( Archive & ar, const unsigned int version ) const {
          ar << boost::serialization::base_object < image > (*this);
          ar << rows_;
          ar << cols_;
          ar << path_;
          ar << sighdu_;
          ar << nsehdu_;
          ar << skyhdu_;
          ar << sky_;
      }
      template < class Archive >
      void load ( Archive & ar, const unsigned int version ) {
          ar >> boost::serialization::base_object < image > (*this);
          ar >> rows_;
          ar >> cols_;
          ar >> path_;
          ar >> sighdu_;
          ar >> nsehdu_;
          ar >> skyhdu_;
          ar >> sky_;
      }
      BOOST_SERIALIZATION_SPLIT_MEMBER()
    
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

      template < class Archive >
      void save ( Archive & ar, const unsigned int version ) const {
          ar << boost::serialization::base_object < image > (*this);
          ar << spec_props_;
          ar << psf_props_;
          ar << rows_;
          ar << cols_;
          ar << signal_;
          ar << noise_;
          ar << invcov_;
          ar << sky_;
      }
      template < class Archive >
      void load ( Archive & ar, const unsigned int version ) {
          ar >> boost::serialization::base_object < image > (*this);
          ar >> spec_props_;
          ar >> psf_props_;
          ar >> rows_;
          ar >> cols_;
          ar >> signal_;
          ar >> noise_;
          ar >> invcov_;
          ar >> sky_;
      }
      BOOST_SERIALIZATION_SPLIT_MEMBER()
    
      boost::property_tree::ptree spec_props_;
      boost::property_tree::ptree psf_props_;
      size_t rows_;
      size_t cols_;
      vector_double signal_;
      vector_double noise_;
      vector_double invcov_;
      std::vector < bool > sky_;
    
  };
  
  
}

