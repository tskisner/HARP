// @COPYRIGHT@

#ifndef HARP_IMAGE_HPP
#define HARP_IMAGE_HPP



namespace harp {

  class image : public boost::enable_shared_from_this < image > {
    
    public :
      image ( boost::property_tree::ptree const & props );
      virtual ~image ( ) { }
      void cleanup ( );

      virtual boost::property_tree::ptree serialize ( ) {
        HARP_THROW( "fell through to virtual method" );
        return boost::property_tree::ptree();
      }
      
      virtual size_t rows ( ) {
        HARP_THROW( "fell through to virtual method" );
        return 0;
      }
      
      virtual size_t cols ( ) {
        HARP_THROW( "fell through to virtual method" );
        return 0;
      }

      virtual std::vector < bool > sky ( ) {
        HARP_THROW( "fell through to virtual method" );
        return std::vector < bool > ();
      }
      
      virtual void read ( matrix_local & data ) {
        HARP_THROW( "fell through to virtual method" );
        return;
      }

      virtual void write ( std::string const & path, matrix_local & data ) {
        HARP_THROW( "fell through to virtual method" );
        return;
      }
      
      virtual void read_noise ( matrix_local & data ) {
        HARP_THROW( "fell through to virtual method" );
        return;
      }

      virtual void write_noise ( std::string const & path, matrix_local & data ) {
        HARP_THROW( "fell through to virtual method" );
        return;
      }

      std::string format ( );
      
      static image * create ( boost::property_tree::ptree const & props );
      
      image * clone ( );
      
      template < class T >
      boost::shared_ptr < T > shared_ref ( ) {
        return boost::shared_polymorphic_downcast < T, image > ( shared_from_this ( ) );
      }
      
      template < class T >
      boost::weak_ptr < T > weak_ref ( ) {
        boost::shared_ptr < T > temp = boost::shared_polymorphic_downcast < T, image > ( shared_from_this() );
        return boost::weak_ptr < T > ( temp );
      }
      
    private :
    
      std::string format_;
      boost::property_tree::ptree props_;
      
      
  };
  
  typedef boost::shared_ptr < harp::image > image_p;
  typedef boost::weak_ptr < harp::image > image_wp;

}


#endif

