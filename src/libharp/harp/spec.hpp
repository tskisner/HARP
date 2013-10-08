// @COPYRIGHT@

#ifndef HARP_SPEC_HPP
#define HARP_SPEC_HPP


namespace harp {

  class spec : public boost::enable_shared_from_this < spec > {

    friend class boost::serialization::access;
    
    public :

      spec ( ) {
        format_ = "";
      }

      spec ( boost::property_tree::ptree const & props );
      
      virtual ~spec ( ) { }

      virtual size_t n_spec ( ) {
        HARP_THROW( "fell through to virtual method" );
        return 0;
      }

      virtual size_t n_lambda ( ) {
        HARP_THROW( "fell through to virtual method" );
        return 0;
      }

      virtual void read ( vector_double & data, vector_double & lambda, std::vector < bool > & sky ) {
        HARP_THROW( "fell through to virtual method" );
        return;
      }

      virtual void write ( std::string const & path, vector_double & data, vector_double const & lambda, std::vector < bool > const & sky ) {
        HARP_THROW( "fell through to virtual method" );
        return;
      }

      void read ( matrix_double & data, vector_double & lambda, std::vector < bool > & sky ) {

        size_t nspec = n_spec();
        size_t nlambda = n_lambda();

        size_t nelem = nspec * nlambda;

        data.resize ( nspec, nlambda );

        vector_double tempdata ( nelem );

        read ( tempdata, lambda, sky );

        for ( size_t i = 0; i < nspec; ++i ) {
          for ( size_t j = 0; j < nlambda; ++j ) {
            data( i, j ) = tempdata[ i * nlambda + j ];
          }
        }

        return;
      }

      void write ( std::string const & path, matrix_double & data, vector_double & lambda, std::vector < bool > & sky ) {

        size_t nspec = n_spec();
        size_t nlambda = n_lambda();

        size_t nelem = nspec * nlambda;

        if ( ( nlambda != data.size2() ) || ( nspec != data.size1() ) ) {
          HARP_THROW( "data size does not match spec dimensions" );
        }

        vector_double tempdata ( nelem );

        for ( size_t i = 0; i < nspec; ++i ) {
          for ( size_t j = 0; j < nlambda; ++j ) {
            tempdata[ i * nlambda + j ] = data( i, j );
          }
        }

        write ( path, tempdata, lambda, sky );

        return;
      }

      std::string format ( );
      
      static spec * create ( boost::property_tree::ptree const & props );
      
      spec * clone ( );
      
      template < class T >
      boost::shared_ptr < T > shared_ref ( ) {
        return boost::dynamic_pointer_cast < T, spec > ( shared_from_this ( ) );
      }
      
      template < class T >
      boost::weak_ptr < T > weak_ref ( ) {
        boost::shared_ptr < T > temp = boost::dynamic_pointer_cast < T, spec > ( shared_from_this() );
        return boost::weak_ptr < T > ( temp );
      }
      
    private :

      template < class Archive >
      void serialize ( Archive & ar, const unsigned int version ) {
        ar & BOOST_SERIALIZATION_NVP(format_);
        ar & BOOST_SERIALIZATION_NVP(props_);
        return;
      }
    
      std::string format_;
      boost::property_tree::ptree props_;
      
      
  };

  //BOOST_CLASS_EXPORT_GUID(spec, "spec")

  BOOST_SERIALIZATION_ASSUME_ABSTRACT(spec)

  BOOST_SERIALIZATION_SHARED_PTR(spec)
  
  typedef boost::shared_ptr < harp::spec > spec_p;
  typedef boost::weak_ptr < harp::spec > spec_wp;

}


#endif

