// @COPYRIGHT@

#ifndef HARP_SPECTRUM_HPP
#define HARP_SPECTRUM_HPP


namespace harp {

  class spectrum : public boost::enable_shared_from_this < spectrum > {
    
    public :
      spectrum ( std::string const & format, std::map < std::string, std::string > const & params );
      virtual ~spectrum ( ) { }
      void cleanup ( );
      
      virtual size_t size ( ) { return 0; }
      virtual void read ( harp::data_vec_view & data ) { return; }
      virtual void write ( std::string const & path, harp::data_vec_view & data ) { return; }

      std::string format ( );
      
      static spectrum * create ( std::string const & format, std::map < std::string, std::string > const & params );
      
      spectrum * clone ( );
      
      template < class T >
      boost::shared_ptr < T > shared_ref ( ) {
        return boost::shared_polymorphic_downcast < T, spectrum > ( shared_from_this ( ) );
      }
      
      template < class T >
      boost::weak_ptr < T > weak_ref ( ) {
        boost::shared_ptr < T > temp = boost::shared_polymorphic_downcast < T, spectrum > ( shared_from_this() );
        return boost::weak_ptr < T > ( temp );
      }
      
    private :
    
      std::string format_;
      std::map < std::string, std::string > params_;
      
      
  };
  
  typedef boost::shared_ptr < harp::spectrum > spectrum_p;
  typedef boost::weak_ptr < harp::spectrum > spectrum_wp;

}


#endif

