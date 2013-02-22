// @COPYRIGHT@

namespace harp {
  
  
  class spec_boss_specter : public spec {
    
    public :
      spec_boss_specter ( boost::property_tree::ptree const & props );
      ~spec_boss_specter ( );

      boost::property_tree::ptree serialize ( );
      
      size_t nspectrum ( ) { return nspec_; }
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