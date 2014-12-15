/*
  High Performance Astrophysical Reconstruction and Processing (HARP)
  (c) 2014-2015, The Regents of the University of California, 
  through Lawrence Berkeley National Laboratory.  See top
  level LICENSE file for details.
*/

#ifndef HARP_INTERNAL_SPEC_PLUGINS_HPP
#define HARP_INTERNAL_SPEC_PLUGINS_HPP

namespace harp {
  
  class spec_sim : public spec {

    friend class boost::serialization::access;
    
    public :

      spec_sim ( ) : spec () {
        size_ = 0;
        nspec_ = 0;
        nlambda_ = 0;
        background_ = 0.0;
        atmpeak_ = 0.0;
        objpeak_ = 0.0;
        atmspace_ = 0;
        objspace_ = 0;
        skymod_ = 0;
        first_lambda_ = 0.0;
        last_lambda_ = 0.0;
      }

      spec_sim ( boost::property_tree::ptree const & props );
      
      ~spec_sim ( );

      void sky_truth ( vector_double & data ) const;
      
      // overloaded virtual methods from base class

      size_t n_spec ( ) const { return nspec_; }

      size_t n_lambda ( ) const { return nlambda_; }

      void values ( vector_double & data ) const;

      void inv_variance ( vector_double & data ) const;

      void lambda ( vector_double & lambda_vals ) const;
    
    private :

      template < class Archive >
      void serialize ( Archive & ar, const unsigned int version ) {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(spec);
        ar & BOOST_SERIALIZATION_NVP(size_);
        ar & BOOST_SERIALIZATION_NVP(nspec_);
        ar & BOOST_SERIALIZATION_NVP(nlambda_);
        ar & BOOST_SERIALIZATION_NVP(background_);
        ar & BOOST_SERIALIZATION_NVP(atmpeak_);
        ar & BOOST_SERIALIZATION_NVP(objpeak_);
        ar & BOOST_SERIALIZATION_NVP(atmspace_);
        ar & BOOST_SERIALIZATION_NVP(objspace_);
        ar & BOOST_SERIALIZATION_NVP(skymod_);
        ar & BOOST_SERIALIZATION_NVP(first_lambda_);
        ar & BOOST_SERIALIZATION_NVP(last_lambda_);
        return;
      }
    
      size_t size_;
      size_t nspec_;
      size_t nlambda_;
      double background_;
      double atmpeak_;
      double objpeak_;
      size_t atmspace_;
      size_t objspace_;
      size_t skymod_;
      double first_lambda_;
      double last_lambda_;
    
  };

  BOOST_SERIALIZATION_SHARED_PTR(spec_sim)

  spec * spec_sim_create ( boost::property_tree::ptree const & props );
  

  class spec_fits : public spec {

    friend class boost::serialization::access;
    
    public :

      spec_fits ( ) : spec () {
        nspec_ = 0;
        nlambda_ = 0;
        nglobal_ = 0;
        path_ = "";
        spechdu_ = -1;
        invvarhdu_ = -1;
        lambdahdu_ = -1;
      }

      spec_fits ( boost::property_tree::ptree const & props );
      
      ~spec_fits ( );

      static void write ( std::string const & path, vector_double const & data, vector_double const & invvar, vector_double const & lambda );

      // overloaded virtual methods from base class
      
      size_t n_spec ( ) const { return nspec_; }

      size_t n_lambda ( ) const { return nlambda_; }

      void values ( vector_double & data ) const;

      void inv_variance ( vector_double & data ) const;

      void lambda ( vector_double & lambda_vals ) const;
    
    private :

      template < class Archive >
      void serialize ( Archive & ar, const unsigned int version ) {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(spec);
        ar & BOOST_SERIALIZATION_NVP(nspec_);
        ar & BOOST_SERIALIZATION_NVP(nlambda_);
        ar & BOOST_SERIALIZATION_NVP(nglobal_);
        ar & BOOST_SERIALIZATION_NVP(path_);
        ar & BOOST_SERIALIZATION_NVP(spechdu_);
        ar & BOOST_SERIALIZATION_NVP(invvarhdu_);
        ar & BOOST_SERIALIZATION_NVP(lambdahdu_);
        return;
      }
    
      size_t nspec_;
      size_t nlambda_;
      size_t nglobal_;
      std::string path_;
      int spechdu_;
      int invvarhdu_;
      int lambdahdu_;
    
  };

  BOOST_SERIALIZATION_SHARED_PTR(spec_fits)

  spec * spec_fits_create ( boost::property_tree::ptree const & props );

  
}

#endif