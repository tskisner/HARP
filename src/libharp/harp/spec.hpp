// @COPYRIGHT@

#ifndef HARP_SPEC_HPP
#define HARP_SPEC_HPP


namespace harp {

  class spec : public boost::enable_shared_from_this < spec > {
    
    public :
      spec ( std::string const & format, boost::property_tree::ptree const & props );
      virtual ~spec ( ) { }
      void cleanup ( );

      virtual boost::property_tree::ptree serialize ( ) { return boost::property_tree::ptree(); }

      virtual size_t nspectrum ( ) { return 0; }
      virtual size_t spectrum_size ( ) { return 0; }
      
      virtual void read ( matrix_dist & data ) { return; }
      virtual void write ( std::string const & path, matrix_dist & data ) { return; }

      std::string format ( );
      
      static spec * create ( std::string const & format, boost::property_tree::ptree const & props );
      
      spec * clone ( );
      
      template < class T >
      boost::shared_ptr < T > shared_ref ( ) {
        return boost::shared_polymorphic_downcast < T, spec > ( shared_from_this ( ) );
      }
      
      template < class T >
      boost::weak_ptr < T > weak_ref ( ) {
        boost::shared_ptr < T > temp = boost::shared_polymorphic_downcast < T, spec > ( shared_from_this() );
        return boost::weak_ptr < T > ( temp );
      }
      
    private :
    
      std::string format_;
      boost::property_tree::ptree props_;
      
      
  };
  
  typedef boost::shared_ptr < harp::spec > spec_p;
  typedef boost::weak_ptr < harp::spec > spec_wp;

}


#endif

