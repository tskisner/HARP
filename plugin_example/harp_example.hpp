// @COPYRIGHT@

namespace harp {
  
  class spec_example : public spec {

    friend class boost::serialization::access;
    
    public :

      spec_example ( ) : spec () {
      }

      spec_example ( boost::property_tree::ptree const & props );
      
      ~spec_example ( );

      // overloaded virtual methods from base class

      boost::property_tree::ptree metadata ( ) const;

      size_t n_spec ( ) const { return nspec_; }

      size_t n_lambda ( ) const { return nlambda_; }

      void values ( vector_double & data ) const;

      void lambda ( vector_double & lambda_vals ) const;

      void targets ( std::vector < target > & target_list ) const;
    
    private :

      template < class Archive >
      void serialize ( Archive & ar, const unsigned int version ) {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(spec);
        ar & BOOST_SERIALIZATION_NVP(nspec_);
        ar & BOOST_SERIALIZATION_NVP(nlambda_);
        return;
      }
    
      size_t nspec_;
      size_t nlambda_;
    
  };

  BOOST_SERIALIZATION_SHARED_PTR(spec_example)



  
}

