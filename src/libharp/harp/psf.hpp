// @COPYRIGHT@

#ifndef HARP_PSF_HPP
#define HARP_PSF_HPP


namespace harp {

  class psf : public boost::enable_shared_from_this < psf > {

    friend class boost::serialization::access;
    
    public :
      psf ( boost::property_tree::ptree const & props );

      virtual ~psf ( ) { }

      virtual size_t n_spec ( ) {
        HARP_THROW( "fell through to virtual method" );
        return 0;
      }
      
      virtual size_t n_lambda ( ) {
        HARP_THROW( "fell through to virtual method" );
        return 0;
      }

      virtual size_t img_rows ( ) {
        HARP_THROW( "fell through to virtual method" );
        return 0;
      }

      virtual size_t img_cols ( ) {
        HARP_THROW( "fell through to virtual method" );
        return 0;
      }

      virtual vector_double lambda ( ) {
        HARP_THROW( "fell through to virtual method" );
        return vector_double();
      }

      virtual void response ( size_t spec, size_t lambda, size_t & x_offset, size_t & y_offset, matrix_double & patch ) {
        HARP_THROW( "fell through to virtual method" );
        return;
      }

      virtual void write ( std::string const & path ) {
        HARP_THROW( "fell through to virtual method" );
        return;
      }
      
      virtual void project ( size_t first_spec, size_t last_spec, size_t first_lambda, size_t last_lambda, matrix_double & A ) {
        

        return;
      }

      virtual void project_transpose ( size_t first_spec, size_t last_spec, size_t first_lambda, size_t last_lambda, matrix_double_sparse & AT ) {

        
        
        return;
      }
      
      std::string format ( );
      
      static psf * create ( boost::property_tree::ptree const & props );
      
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

      template < class Archive >
      void save ( Archive & ar, const unsigned int version ) const {
          ar << format_;
          ar << props_;
      }
      template < class Archive >
      void load ( Archive & ar, const unsigned int version ) {
          ar >> format_;
          ar >> props_;
      }
      BOOST_SERIALIZATION_SPLIT_MEMBER()
    
      std::string format_;
      boost::property_tree::ptree props_;
      
      
  };
  
  typedef boost::shared_ptr < harp::psf > psf_p;
  typedef boost::weak_ptr < harp::psf > psf_wp;

}


#endif

