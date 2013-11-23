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
      
      virtual size_t n_rows ( ) const {
        HARP_THROW( "fell through to virtual method" );
        return 0;
      }
      
      virtual size_t n_cols ( ) const {
        HARP_THROW( "fell through to virtual method" );
        return 0;
      }
      
      virtual void values ( vector_double & data ) const {
        HARP_THROW( "fell through to virtual method" );
        return;
      }

      virtual void inv_variance ( vector_double & invvar ) const {
        HARP_THROW( "fell through to virtual method" );
        return;
      }

      virtual boost::property_tree::ptree metadata ( ) const {
        HARP_THROW( "fell through to virtual method" );
        return boost::property_tree::ptree();
      }
      
      void values ( matrix_double & data ) const {

        size_t imgrows = n_rows();
        size_t imgcols = n_cols();

        size_t nelem = imgrows * imgcols;

        data.resize ( imgrows, imgcols );

        vector_double tempdata ( nelem );

        values ( tempdata );

        for ( size_t i = 0; i < imgcols; ++i ) {
          for ( size_t j = 0; j < imgrows; ++j ) {
            data( j, i ) = tempdata[ i * imgrows + j ];
          }
        }

        return;
      }

      void inv_variance ( matrix_double & invvar ) const {

        size_t imgrows = n_rows();
        size_t imgcols = n_cols();

        size_t nelem = imgrows * imgcols;

        invvar.resize ( imgrows, imgcols );

        vector_double tempvar ( nelem );

        inv_variance ( tempvar );

        for ( size_t i = 0; i < imgcols; ++i ) {
          for ( size_t j = 0; j < imgrows; ++j ) {
            invvar( j, i ) = tempvar[ i * imgrows + j ];
          }
        }

        return;
      }

      std::string format ( ) const;
      
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

