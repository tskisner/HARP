// @COPYRIGHT@

namespace harp {
  
  class spec_sim : public spec {
    
    public :
      spec_sim ( boost::property_tree::ptree const & props );
      ~spec_sim ( );

      boost::property_tree::ptree serialize ( );
      
      size_t nspec ( ) { return nspec_; }
      size_t nlambda ( ) { return nlambda_; }

      void read ( matrix_dist & data, std::vector < double > & lambda, std::vector < bool > & sky );
      
      void write ( std::string const & path, matrix_dist & data, std::vector < double > const & lambda, std::vector < bool > const & sky );

      void sky_truth ( matrix_dist & data );
    
    private :
    
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
  

  class spec_specter : public spec {
    
    public :
      spec_specter ( boost::property_tree::ptree const & props );
      ~spec_specter ( );

      boost::property_tree::ptree serialize ( );
      
      size_t nspec ( ) { return nspec_; }
      size_t nlambda ( ) { return nlambda_; }

      void read ( matrix_dist & data, std::vector < double > & lambda, std::vector < bool > & sky );

      void write ( std::string const & path, matrix_dist & data, std::vector < double > const & lambda, std::vector < bool > const & sky );
    
    private :
    
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

