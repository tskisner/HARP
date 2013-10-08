// @COPYRIGHT@

#ifndef HARP_PSF_HPP
#define HARP_PSF_HPP


namespace harp {

  class psf : public boost::enable_shared_from_this < psf > {
    
    public :
      psf ( boost::property_tree::ptree const & props );
      virtual ~psf ( ) { }
      void cleanup ( );

      virtual boost::property_tree::ptree serialize ( ) { 
        HARP_THROW( "fell through to virtual method" );
        return boost::property_tree::ptree();
      }
      
      virtual size_t nspec ( ) {
        HARP_THROW( "fell through to virtual method" );
        return 0;
      }
      
      virtual size_t nlambda ( ) {
        HARP_THROW( "fell through to virtual method" );
        return 0;
      }

      virtual size_t pixrows ( ) {
        HARP_THROW( "fell through to virtual method" );
        return 0;
      }

      virtual size_t pixcols ( ) {
        HARP_THROW( "fell through to virtual method" );
        return 0;
      }

      virtual std::vector < double > lambda ( ) {
        HARP_THROW( "fell through to virtual method" );
        return std::vector < double > ();
      }
      
      virtual void projection ( size_t first_spec, size_t last_spec, size_t first_lambda, size_t last_lambda, matrix_sparse & AT ) {
        HARP_THROW( "fell through to virtual method" );
        return;
      }
      
      std::string format ( );
      
      static psf * create ( boost::property_tree::ptree const & props );
      
      psf * clone ( );
      
      template < class T >
      boost::shared_ptr < T > shared_ref ( ) {
        return boost::dynamic_pointer_cast < T, psf > ( shared_from_this ( ) );
      }
      
      template < class T >
      boost::weak_ptr < T > weak_ref ( ) {
        boost::shared_ptr < T > temp = boost::dynamic_pointer_cast < T, psf > ( shared_from_this() );
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

