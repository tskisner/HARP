// @COPYRIGHT@

#ifndef HARP_PSF_HPP
#define HARP_PSF_HPP


namespace harp {

  class psf : public boost::enable_shared_from_this < psf > {
    
    public :
      psf ( std::string const & format, boost::property_tree::ptree const & props );
      virtual ~psf ( ) { }
      void cleanup ( );

      virtual boost::property_tree::ptree serialize ( ) { return boost::property_tree::ptree(); }
      
      virtual size_t nspec ( ) { return 0; }
      
      virtual size_t specsize ( ) { return 0; }

      virtual size_t pixrows ( ) { return 0; }

      virtual size_t pixcols ( ) { return 0; }
      
      virtual void projection ( size_t firstbin, size_t lastbin, matrix_sparse & AT ) { return; }
      
      std::string format ( );
      
      static psf * create ( std::string const & format, boost::property_tree::ptree const & props );
      
      psf * clone ( );
      
      template < class T >
      boost::shared_ptr < T > shared_ref ( ) {
        return boost::shared_polymorphic_downcast < T, psf > ( shared_from_this ( ) );
      }
      
      template < class T >
      boost::weak_ptr < T > weak_ref ( ) {
        boost::shared_ptr < T > temp = boost::shared_polymorphic_downcast < T, psf > ( shared_from_this() );
        return boost::weak_ptr < T > ( temp );
      }
      
    private :
    
      std::string format_;
      boost::property_tree::ptree props_;
      
      
  };
  
  typedef boost::shared_ptr < harp::psf > psf_p;
  typedef boost::weak_ptr < harp::psf > psf_wp;

}


#endif

