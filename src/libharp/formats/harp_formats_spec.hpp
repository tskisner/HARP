// @COPYRIGHT@

namespace harp {
  
  class spec_sim : public spec {

    friend class boost::serialization::access;
    
    public :
      spec_sim ( boost::property_tree::ptree const & props );
      ~spec_sim ( );
      
      size_t n_spec ( ) { return nspec_; }

      size_t n_lambda ( ) { return nlambda_; }

      void read ( vector_double & data, vector_double & lambda, std::vector < bool > & sky );
      
      void write ( std::string const & path, vector_double & data, vector_double const & lambda, std::vector < bool > const & sky );

      void sky_truth ( vector_double & data );
    
    private :

      template < class Archive >
      void save ( Archive & ar, const unsigned int version ) const {
          ar << boost::serialization::base_object < spec > (*this);
          ar << size_;
          ar << nspec_;
          ar << nlambda_;
          ar << background_;
          ar << atmpeak_;
          ar << objpeak_;
          ar << atmspace_;
          ar << objspace_;
          ar << skymod_;
          ar << first_lambda_;
          ar << last_lambda_;
      }
      template < class Archive >
      void load ( Archive & ar, const unsigned int version ) {
          ar >> boost::serialization::base_object < spec > (*this);
          ar >> size_;
          ar >> nspec_;
          ar >> nlambda_;
          ar >> background_;
          ar >> atmpeak_;
          ar >> objpeak_;
          ar >> atmspace_;
          ar >> objspace_;
          ar >> skymod_;
          ar >> first_lambda_;
          ar >> last_lambda_;
      }
      BOOST_SERIALIZATION_SPLIT_MEMBER()
    
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


  void specter_read_sky ( fitsfile * fp, std::vector < bool > & sky );


  void specter_write_sky ( fitsfile * fp, std::vector < bool > const & sky );
  

  class spec_specter : public spec {

    friend class boost::serialization::access;
    
    public :
      spec_specter ( boost::property_tree::ptree const & props );
      ~spec_specter ( );
      
      size_t n_spec ( ) { return nspec_; }

      size_t n_lambda ( ) { return nlambda_; }

      void read ( vector_double & data, vector_double & lambda, std::vector < bool > & sky );
      
      void write ( std::string const & path, vector_double & data, vector_double const & lambda, std::vector < bool > const & sky );
    
    private :

      template < class Archive >
      void save ( Archive & ar, const unsigned int version ) const {
          ar << boost::serialization::base_object < spec > (*this);
          ar << size_;
          ar << nspec_;
          ar << nlambda_;
          ar << nglobal_;
          ar << path_;
          ar << objonly_;
          ar << spechdu_;
          ar << lambdahdu_;
          ar << targethdu_;
      }
      template < class Archive >
      void load ( Archive & ar, const unsigned int version ) {
          ar >> boost::serialization::base_object < spec > (*this);
          ar >> size_;
          ar >> nspec_;
          ar >> nlambda_;
          ar >> nglobal_;
          ar >> path_;
          ar >> objonly_;
          ar >> spechdu_;
          ar >> lambdahdu_;
          ar >> targethdu_;
      }
      BOOST_SERIALIZATION_SPLIT_MEMBER()
    
      size_t size_;
      size_t nspec_;
      size_t nlambda_;
      size_t nglobal_;
      std::string path_;
      bool objonly_;
      int spechdu_;
      int lambdahdu_;
      int targethdu_;
    
  };
  
  
}

