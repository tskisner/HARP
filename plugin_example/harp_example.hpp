// @COPYRIGHT@

// To get the usual type definitions from HARP...
#ifdef USE_MPI
#  include <harp_mpi.hpp>
#else
#  include <harp.hpp>
#endif


#include <iostream>


// For this example, we define 3 "fake" classes (one spec, one image, and one psf) and
// define the functions which allow registration of these 3 classes as external plugins.


namespace harp {

  // example spec class which does nothing.  This is just to show how to compile
  // an external plugin.
  
  class spec_example : public spec {

    // this is needed for serialization
    friend class boost::serialization::access;
    
    public :

      // default constructor is needed by serialization
      spec_example ( ) : spec () { }

      // this is the "real" constructor which takes a property tree
      spec_example ( boost::property_tree::ptree const & props );
      
      ~spec_example ( ) { }

      // overloaded virtual methods from base class

      boost::property_tree::ptree metadata ( ) const { return boost::property_tree::ptree(); }

      size_t n_spec ( ) const { return nspec_; }

      size_t n_lambda ( ) const { return nlambda_; }

      void values ( vector_double & data ) const { return; }

      void lambda ( vector_double & lambda_vals ) const { return; }

      void targets ( std::vector < obs_target > & target_list ) const { return; }
    
    private :

      // you must define a serialize function which declares the member variables to
      // dump / restore.

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

  // This macro allows serialization of a shared_ptr to the object, and must be defined.

  BOOST_SERIALIZATION_SHARED_PTR(spec_example)

  // This function should allocate and return a new instance of the plugin.  It will be called
  // by the plugin registry when creating instances of your plugin.  The name of the function must
  // be of the form:
  //
  // <type> * harp::<type>_<name>_create ( boost::property_tree::ptree const & props )

  spec * spec_example_create ( boost::property_tree::ptree const & props );


  // Here are some more example instances of the psf and image classes


  class psf_example : public psf {

    friend class boost::serialization::access;
    
    public :

      psf_example ( ) : psf () { }

      psf_example ( boost::property_tree::ptree const & props );

      ~psf_example ( ) { }

      // overloaded virtual methods from base class

      boost::property_tree::ptree metadata ( ) const { return boost::property_tree::ptree(); }

      size_t n_spec ( ) const { return nspec_; }

      size_t n_lambda ( ) const { return nlambda_; }
      
      size_t img_rows ( ) const { return rows_; }
      
      size_t img_cols ( ) const { return cols_; }
      
      vector_double lambda ( ) const { return vector_double(); }

      void response ( size_t spec, size_t lambda, size_t & x_offset, size_t & y_offset, matrix_double & patch ) const { return; }

      size_t response_nnz_estimate ( ) const { return 0; }

    private :

      template < class Archive >
      void serialize ( Archive & ar, const unsigned int version ) {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(psf);
        ar & BOOST_SERIALIZATION_NVP(nspec_);
        ar & BOOST_SERIALIZATION_NVP(nlambda_);
        ar & BOOST_SERIALIZATION_NVP(rows_);
        ar & BOOST_SERIALIZATION_NVP(cols_);
        return;
      }

      size_t nspec_;
      size_t nlambda_;
      size_t rows_;
      size_t cols_;
      
  };

  BOOST_SERIALIZATION_SHARED_PTR(psf_example)

  psf * psf_example_create ( boost::property_tree::ptree const & props );


  class image_example : public image {

    friend class boost::serialization::access;
    
    public :

      image_example ( ) : image () { }
      
      image_example ( boost::property_tree::ptree const & props );
      
      ~image_example ( ) { }

      // overloaded virtual methods from base class
      
      size_t n_rows ( ) const { return rows_; }
      
      size_t n_cols ( ) const { return cols_; }
      
      void values ( vector_double & data ) const { return; }

      void inv_variance ( vector_double & invvar ) const { return; }

      boost::property_tree::ptree metadata ( ) const { return boost::property_tree::ptree(); }

    private :

      template < class Archive >
      void serialize ( Archive & ar, const unsigned int version ) {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(image);
        ar & BOOST_SERIALIZATION_NVP(rows_);
        ar & BOOST_SERIALIZATION_NVP(cols_);
        return;
      }
    
      size_t rows_;
      size_t cols_;
    
  };

  BOOST_SERIALIZATION_SHARED_PTR(image_example)

  image * image_example_create ( boost::property_tree::ptree const & props );

}


// This global function in the shared object must be unique and outside the harp namespace.
// this is the function actually called by the plugin registry.

extern "C" {
  void initialize ( void * registry );
}


