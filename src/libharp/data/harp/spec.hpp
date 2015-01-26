/*
  High Performance Astrophysical Reconstruction and Processing (HARP)
  (c) 2014-2015, The Regents of the University of California, 
  through Lawrence Berkeley National Laboratory.  See top
  level LICENSE file for details.
*/

#ifndef HARP_SPEC_HPP
#define HARP_SPEC_HPP


namespace harp {

  class spec : public boost::enable_shared_from_this < spec > {

    friend class boost::serialization::access;
    
    public :

      spec ( );
      
      spec ( std::string const & type, boost::property_tree::ptree const & props );
      
      virtual ~spec ( );

      virtual size_t n_spec ( ) const;

      virtual size_t n_lambda ( ) const;

      virtual void values ( vector_double & data ) const;

      virtual void inv_variance ( vector_double & data ) const;

      virtual void lambda ( vector_double & lambda_vals ) const;

      void values ( matrix_double & data ) const;

      void inv_variance ( matrix_double & data ) const;

      std::string type ( ) const;
      
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
        ar & BOOST_SERIALIZATION_NVP(type_);
        ar & BOOST_SERIALIZATION_NVP(props_);
        return;
      }
    
      std::string type_;
      boost::property_tree::ptree props_;
      
      
  };

  BOOST_SERIALIZATION_ASSUME_ABSTRACT(spec)

  BOOST_SERIALIZATION_SHARED_PTR(spec)
  
  typedef boost::shared_ptr < harp::spec > spec_p;
  typedef boost::weak_ptr < harp::spec > spec_wp;

}


#endif

