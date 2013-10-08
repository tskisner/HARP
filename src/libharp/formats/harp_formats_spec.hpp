// @COPYRIGHT@

namespace harp {
  
  class spec_sim : public spec {

    friend class boost::serialization::access;
    
    public :

      spec_sim ( ) : spec () {
        size_ = 0;
        nspec_ = 0;
        nlambda_ = 0;
        background_ = 0.0;
        atmpeak_ = 0.0;
        objpeak_ = 0.0;
        atmspace_ = 0;
        objspace_ = 0;
        skymod_ = 0;
        first_lambda_ = 0.0;
        last_lambda_ = 0.0;
      }

      spec_sim ( boost::property_tree::ptree const & props );
      
      ~spec_sim ( );
      
      size_t n_spec ( ) { return nspec_; }

      size_t n_lambda ( ) { return nlambda_; }

      void read ( vector_double & data, vector_double & lambda, std::vector < bool > & sky );
      
      void write ( std::string const & path, vector_double & data, vector_double const & lambda, std::vector < bool > const & sky );

      void sky_truth ( vector_double & data );
    
    private :

      template < class Archive >
      void serialize ( Archive & ar, const unsigned int version ) {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(spec);
        ar & BOOST_SERIALIZATION_NVP(size_);
        ar & BOOST_SERIALIZATION_NVP(nspec_);
        ar & BOOST_SERIALIZATION_NVP(nlambda_);
        ar & BOOST_SERIALIZATION_NVP(background_);
        ar & BOOST_SERIALIZATION_NVP(atmpeak_);
        ar & BOOST_SERIALIZATION_NVP(objpeak_);
        ar & BOOST_SERIALIZATION_NVP(atmspace_);
        ar & BOOST_SERIALIZATION_NVP(objspace_);
        ar & BOOST_SERIALIZATION_NVP(skymod_);
        ar & BOOST_SERIALIZATION_NVP(first_lambda_);
        ar & BOOST_SERIALIZATION_NVP(last_lambda_);
        return;
      }
    
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

  BOOST_SERIALIZATION_SHARED_PTR(spec_sim)



  void specter_read_sky ( fitsfile * fp, std::vector < bool > & sky );


  void specter_write_sky ( fitsfile * fp, std::vector < bool > const & sky );
  

  class spec_specter : public spec {

    friend class boost::serialization::access;
    
    public :

      spec_specter ( ) : spec () {
        nspec_ = 0;
        nlambda_ = 0;
        nglobal_ = 0;
        path_ = "";
        objonly_ = false;
        spechdu_ = -1;
        lambdahdu_ = -1;
        targethdu_ = -1;
      }

      spec_specter ( boost::property_tree::ptree const & props );
      
      ~spec_specter ( );
      
      size_t n_spec ( ) { return nspec_; }

      size_t n_lambda ( ) { return nlambda_; }

      void read ( vector_double & data, vector_double & lambda, std::vector < bool > & sky );
      
      void write ( std::string const & path, vector_double & data, vector_double const & lambda, std::vector < bool > const & sky );
    
    private :

      template < class Archive >
      void serialize ( Archive & ar, const unsigned int version ) {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(spec);
        ar & BOOST_SERIALIZATION_NVP(nspec_);
        ar & BOOST_SERIALIZATION_NVP(nlambda_);
        ar & BOOST_SERIALIZATION_NVP(nglobal_);
        ar & BOOST_SERIALIZATION_NVP(path_);
        ar & BOOST_SERIALIZATION_NVP(objonly_);
        ar & BOOST_SERIALIZATION_NVP(spechdu_);
        ar & BOOST_SERIALIZATION_NVP(lambdahdu_);
        ar & BOOST_SERIALIZATION_NVP(targethdu_);
        return;
      }
    
      size_t nspec_;
      size_t nlambda_;
      size_t nglobal_;
      std::string path_;
      bool objonly_;
      int spechdu_;
      int lambdahdu_;
      int targethdu_;
    
  };

  BOOST_SERIALIZATION_SHARED_PTR(spec_specter)

  
}

