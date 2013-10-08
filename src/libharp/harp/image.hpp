// @COPYRIGHT@

#ifndef HARP_IMAGE_HPP
#define HARP_IMAGE_HPP



namespace harp {

  class image : public boost::enable_shared_from_this < image > {

    friend class boost::serialization::access;
    
    public :

      image ( ) {
        format_ = "";
      }

      image ( boost::property_tree::ptree const & props );
      
      virtual ~image ( ) { }
      
      virtual size_t n_rows ( ) {
        HARP_THROW( "fell through to virtual method" );
        return 0;
      }
      
      virtual size_t n_cols ( ) {
        HARP_THROW( "fell through to virtual method" );
        return 0;
      }
      
      virtual void read ( vector_double & data, vector_double & invvar, std::vector < bool > & sky ) {
        HARP_THROW( "fell through to virtual method" );
        return;
      }

      virtual void write ( std::string const & path, vector_double & data, vector_double & invvar, std::vector < bool > & sky ) {
        HARP_THROW( "fell through to virtual method" );
        return;
      }
      
      void read ( matrix_double & data, matrix_double & invvar, std::vector < bool > & sky ) {

        size_t imgrows = n_rows();
        size_t imgcols = n_cols();

        size_t nelem = imgrows * imgcols;

        data.resize ( imgcols, imgrows );
        invvar.resize ( imgcols, imgrows );

        vector_double tempdata ( nelem );
        vector_double tempvar ( nelem );

        read ( tempdata, tempvar, sky );

        for ( size_t i = 0; i < imgcols; ++i ) {
          for ( size_t j = 0; j < imgrows; ++j ) {
            data( i, j ) = tempdata[ i * imgrows + j ];
            invvar( i, j ) = tempvar[ i * imgrows + j ];
          }
        }

        return;
      }

      void write ( std::string const & path, matrix_double & data, matrix_double & invvar, std::vector < bool > & sky ) {

        size_t imgrows = n_rows();
        size_t imgcols = n_cols();

        size_t nelem = imgrows * imgcols;

        if ( ( imgrows != data.size2() ) || ( imgcols != data.size1() ) ) {
          HARP_THROW( "data size does not match image dimensions" );
        }

        if ( ( imgrows != invvar.size2() ) || ( imgcols != invvar.size1() ) ) {
          HARP_THROW( "inverse variance size does not match image dimensions" );
        }

        vector_double tempdata ( nelem );
        vector_double tempvar ( nelem );

        for ( size_t i = 0; i < imgcols; ++i ) {
          for ( size_t j = 0; j < imgrows; ++j ) {
            tempdata[ i * imgrows + j ] = data( i, j );
            tempvar[ i * imgrows + j ] = invvar( i, j );
          }
        }

        write ( path, tempdata, tempvar, sky );

        return;
      }

      std::string format ( );
      
      static image * create ( boost::property_tree::ptree const & props );
      
      image * clone ( );
      
      template < class T >
      boost::shared_ptr < T > shared_ref ( ) {
        return boost::dynamic_pointer_cast < T, image > ( shared_from_this ( ) );
      }
      
      template < class T >
      boost::weak_ptr < T > weak_ref ( ) {
        boost::shared_ptr < T > temp = boost::dynamic_pointer_cast < T, image > ( shared_from_this() );
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

  BOOST_SERIALIZATION_ASSUME_ABSTRACT(image)

  BOOST_SERIALIZATION_SHARED_PTR(image)
  
  typedef boost::shared_ptr < harp::image > image_p;
  typedef boost::weak_ptr < harp::image > image_wp;

}


#endif

