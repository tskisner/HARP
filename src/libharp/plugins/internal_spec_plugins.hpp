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

      spec_sim ( );

      spec_sim ( boost::property_tree::ptree const & props );
      
      ~spec_sim ( );

      void sky_truth ( vector_double & data ) const;
      
      // overloaded virtual methods from base class

      size_t n_spec ( ) const;

      size_t n_lambda ( ) const;

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

      spec_fits ( );

      spec_fits ( boost::property_tree::ptree const & props );
      
      ~spec_fits ( );

      static void write ( std::string const & path, vector_double const & data, vector_double const & invvar, vector_double const & lambda );

      // overloaded virtual methods from base class
      
      size_t n_spec ( ) const;

      size_t n_lambda ( ) const;

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

  
  
  class spec_desi : public spec {

    friend class boost::serialization::access;
    
    public :

      spec_desi ( );

      spec_desi ( boost::property_tree::ptree const & props );
      
      ~spec_desi ( );

      static void write ( std::string const & path, boost::property_tree::ptree const & meta, vector_double const & data, vector_double const & invvar, vector_double const & lambda );

      // overloaded virtual methods from base class
      
      size_t n_spec ( ) const;

      size_t n_lambda ( ) const;

      void values ( vector_double & data ) const;

      void inv_variance ( vector_double & data ) const;

      void lambda ( vector_double & lambda_vals ) const;

      // desi-specific

      boost::property_tree::ptree meta () const;

      float crval;
      float cdelt;
      std::string airorvac;
      int loglam;
      std::string simfile;
      std::string camera;
      std::string vharp;
      float exptime;
      float rdnoise;
      std::string flavor;
      std::string in_psf;
      std::string in_img;

    private :

      template < class Archive >
      void serialize ( Archive & ar, const unsigned int version ) {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(spec);
        ar & BOOST_SERIALIZATION_NVP(nspec_);
        ar & BOOST_SERIALIZATION_NVP(nlambda_);
        ar & BOOST_SERIALIZATION_NVP(nglobal_);
        ar & BOOST_SERIALIZATION_NVP(path_);
        ar & BOOST_SERIALIZATION_NVP(meta_);
        ar & BOOST_SERIALIZATION_NVP(spechdu_);
        ar & BOOST_SERIALIZATION_NVP(invvarhdu_);
        ar & BOOST_SERIALIZATION_NVP(lambdahdu_);
        ar & BOOST_SERIALIZATION_NVP(reshdu_);
        ar & BOOST_SERIALIZATION_NVP(crval);
        ar & BOOST_SERIALIZATION_NVP(cdelt);
        ar & BOOST_SERIALIZATION_NVP(airorvac);
        ar & BOOST_SERIALIZATION_NVP(loglam);
        ar & BOOST_SERIALIZATION_NVP(simfile);
        ar & BOOST_SERIALIZATION_NVP(camera);
        ar & BOOST_SERIALIZATION_NVP(vharp);
        ar & BOOST_SERIALIZATION_NVP(exptime);
        ar & BOOST_SERIALIZATION_NVP(rdnoise);
        ar & BOOST_SERIALIZATION_NVP(flavor);
        ar & BOOST_SERIALIZATION_NVP(in_psf);
        ar & BOOST_SERIALIZATION_NVP(in_img);
        return;
      }
    
      boost::property_tree::ptree meta_;
      size_t nspec_;
      size_t nlambda_;
      size_t nglobal_;
      std::string path_;
      int spechdu_;
      int invvarhdu_;
      int lambdahdu_;
      int reshdu_;
    
  };

  BOOST_SERIALIZATION_SHARED_PTR(spec_desi)

  spec * spec_desi_create ( boost::property_tree::ptree const & props );


  class spec_desisim : public spec {

    friend class boost::serialization::access;
    
    public :

      spec_desisim ( );

      spec_desisim ( boost::property_tree::ptree const & props );
      
      ~spec_desisim ( );

      // overloaded virtual methods from base class
      
      size_t n_spec ( ) const;

      size_t n_lambda ( ) const;

      void values ( vector_double & data ) const;

      void inv_variance ( vector_double & data ) const;

      void lambda ( vector_double & lambda_vals ) const;

      // desi-specific

      boost::property_tree::ptree meta () const;

      void sky ( vector_double & data ) const;

      float crval;
      float cdelt;
      std::string airorvac;
      int loglam;

    private :

      template < class Archive >
      void serialize ( Archive & ar, const unsigned int version ) {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(spec);
        ar & BOOST_SERIALIZATION_NVP(nspec_);
        ar & BOOST_SERIALIZATION_NVP(nlambda_);
        ar & BOOST_SERIALIZATION_NVP(nglobal_);
        ar & BOOST_SERIALIZATION_NVP(path_);
        ar & BOOST_SERIALIZATION_NVP(meta_);
        ar & BOOST_SERIALIZATION_NVP(objhdu_);
        ar & BOOST_SERIALIZATION_NVP(skyhdu_);
        ar & BOOST_SERIALIZATION_NVP(crval);
        ar & BOOST_SERIALIZATION_NVP(cdelt);
        ar & BOOST_SERIALIZATION_NVP(airorvac);
        ar & BOOST_SERIALIZATION_NVP(loglam);
        return;
      }
    
      boost::property_tree::ptree meta_;
      size_t nspec_;
      size_t nlambda_;
      size_t nglobal_;
      std::string path_;
      int objhdu_;
      int skyhdu_;
    
  };

  BOOST_SERIALIZATION_SHARED_PTR(spec_desisim)

  spec * spec_desisim_create ( boost::property_tree::ptree const & props );

}

#endif
