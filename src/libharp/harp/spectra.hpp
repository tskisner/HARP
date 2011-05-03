// @COPYRIGHT@

#ifndef HARP_SPECTRA_HPP
#define HARP_SPECTRA_HPP


namespace harp {

  class spectra : public boost::enable_shared_from_this < spectra > {
    
    public :
      spectra ( std::string const & format, std::map < std::string, std::string > const & params );
      virtual ~spectra ( ) { }
      void cleanup ( );
      
      virtual size_t size ( ) { return 0; }
      virtual void read ( size_t startrow, size_t startcol, harp::data_vec & data ) { return; }
      virtual void write ( std::string const & path, harp::data_vec & data ) { return; }

      std::string format ( );
      
      static spectra * create ( std::string const & format, std::map < std::string, std::string > const & params );
      
      spectra * clone ( );
      
      template < class T >
      boost::shared_ptr < T > shared_ref ( ) {
        return boost::shared_polymorphic_downcast < T, spectra > ( shared_from_this ( ) );
      }
      
      template < class T >
      boost::weak_ptr < T > weak_ref ( ) {
        boost::shared_ptr < T > temp = boost::shared_polymorphic_downcast < T, spectra > ( shared_from_this() );
        return boost::weak_ptr < T > ( temp );
      }
      
    private :
    
      std::string format_;
      std::map < std::string, std::string > params_;
      
      
  };
  
  typedef boost::shared_ptr < harp::spectra > spectra_p;
  typedef boost::weak_ptr < harp::spectra > spectra_wp;

}


#endif

