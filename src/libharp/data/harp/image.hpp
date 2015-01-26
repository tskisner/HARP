/*
  High Performance Astrophysical Reconstruction and Processing (HARP)
  (c) 2014-2015, The Regents of the University of California, 
  through Lawrence Berkeley National Laboratory.  See top
  level LICENSE file for details.
*/

#ifndef HARP_IMAGE_HPP
#define HARP_IMAGE_HPP



namespace harp {

  class image : public boost::enable_shared_from_this < image > {

    friend class boost::serialization::access;
    
    public :

      image ( );

      image ( std::string const & type, boost::property_tree::ptree const & props );
      
      virtual ~image ( );
      
      virtual size_t n_rows ( ) const;
      
      virtual size_t n_cols ( ) const;
      
      virtual void values ( vector_double & data ) const;

      virtual void inv_variance ( vector_double & invvar ) const;

      virtual void mask ( vector_mask & msk ) const;
      
      void values ( matrix_double & data ) const;

      void inv_variance ( matrix_double & invvar ) const;

      void mask ( matrix_mask & msk ) const;


      std::string type ( ) const;
      
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
        ar & BOOST_SERIALIZATION_NVP(type_);
        ar & BOOST_SERIALIZATION_NVP(props_);
        return;
      }
    
      std::string type_;
      boost::property_tree::ptree props_;
      
      
  };

  BOOST_SERIALIZATION_ASSUME_ABSTRACT(image)

  BOOST_SERIALIZATION_SHARED_PTR(image)
  
  typedef boost::shared_ptr < harp::image > image_p;
  typedef boost::weak_ptr < harp::image > image_wp;

}


#endif

