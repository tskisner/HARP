// @COPYRIGHT@

#ifndef HARP_IMAGE_HPP
#define HARP_IMAGE_HPP



namespace harp {

  class image : public boost::enable_shared_from_this < image > {
    
    public :
      image ( std::string const & format, std::map < std::string, std::string > const & params );
      virtual ~image ( ) { }
      void cleanup ( );
      
      virtual size_t rows ( ) { return 0; }
      virtual size_t cols ( ) { return 0; }
      virtual void rowscols ( size_t & rows, size_t & cols ) { rows = 0; cols = 0; return; }
      virtual void read ( size_t startrow, size_t startcol, harp::dense_rowmat_view & data ) { return; }
      virtual void write ( std::string const & path, size_t startrow, size_t startcol, harp::dense_rowmat_view & data ) { return; }

      std::string format ( );
      
      static image * create ( std::string const & format, std::map < std::string, std::string > const & params );
      
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
      std::map < std::string, std::string > params_;
      
      
  };
  
  typedef boost::shared_ptr < harp::image > image_p;
  typedef boost::weak_ptr < harp::image > image_wp;

}


#endif

