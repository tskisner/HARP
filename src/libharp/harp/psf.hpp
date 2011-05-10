// @COPYRIGHT@

#ifndef HARP_PSF_HPP
#define HARP_PSF_HPP


namespace harp {

  class psf : public boost::enable_shared_from_this < psf > {
    
    public :
      psf ( std::string const & format, std::map < std::string, std::string > const & params );
      virtual ~psf ( ) { }
      void cleanup ( );
      
      virtual size_t npix ( ) { return 0; }
      virtual size_t nspec ( ) { return 0; }
      virtual void read ( size_t firstpix, size_t firstspec, sparse_mat_view & data ) { return; }
      virtual void write ( size_t firstpix, size_t firstspec, sparse_mat_view & data ) { return; }

      std::string format ( );
      
      static psf * create ( std::string const & format, std::map < std::string, std::string > const & params );
      
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
      std::map < std::string, std::string > params_;
      
      
  };
  
  typedef boost::shared_ptr < harp::psf > psf_p;
  typedef boost::weak_ptr < harp::psf > psf_wp;

}


#endif

