// @COPYRIGHT@

#ifndef HARP_INTERNAL_IMAGE_PLUGINS_HPP
#define HARP_INTERNAL_IMAGE_PLUGINS_HPP

namespace harp {
  
  class image_fits : public image {

    friend class boost::serialization::access;
    
    public :

      image_fits ( ) : image () {
        rows_ = 0;
        cols_ = 0;
        path_ = "";
        sighdu_ = -1;
        nsehdu_ = -1;
        skyhdu_ = -1;
      }
      
      image_fits ( boost::property_tree::ptree const & props );
      
      ~image_fits ( );

      void write ( std::string const & path, vector_double & data, vector_double & invvar );

      void write ( std::string const & path, matrix_double & data, matrix_double & invvar );

      // overloaded virtual methods from base class
      
      size_t n_rows ( ) const { return rows_; }
      
      size_t n_cols ( ) const { return cols_; }
      
      void values ( vector_double & data ) const;

      void inv_variance ( vector_double & invvar ) const;

    private :

      template < class Archive >
      void serialize ( Archive & ar, const unsigned int version ) {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(image);
        ar & BOOST_SERIALIZATION_NVP(rows_);
        ar & BOOST_SERIALIZATION_NVP(cols_);
        ar & BOOST_SERIALIZATION_NVP(path_);
        ar & BOOST_SERIALIZATION_NVP(sighdu_);
        ar & BOOST_SERIALIZATION_NVP(nsehdu_);
        ar & BOOST_SERIALIZATION_NVP(skyhdu_);
        ar & BOOST_SERIALIZATION_NVP(sky_);
        return;
      }
    
      size_t rows_;
      size_t cols_;
      std::string path_;
      int sighdu_;
      int nsehdu_;
      int skyhdu_;
      std::vector < bool > sky_;
    
  };

  BOOST_SERIALIZATION_SHARED_PTR(image_fits)

  image * image_fits_create ( boost::property_tree::ptree const & props );


  class image_sim : public image {

    friend class boost::serialization::access;
    
    public :

      image_sim ( ) : image () {
        rows_ = 0;
        cols_ = 0;
      }

      image_sim ( boost::property_tree::ptree const & props );
      
      ~image_sim ( );

      // overloaded virtual methods from base class

      size_t n_rows ( ) const { return rows_; }
      
      size_t n_cols ( ) const { return cols_; }

      void values ( vector_double & data ) const;

      void inv_variance ( vector_double & invvar ) const;

    private :

      template < class Archive >
      void serialize ( Archive & ar, const unsigned int version ) {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(image);
        ar & BOOST_SERIALIZATION_NVP(spec_type_);
        ar & BOOST_SERIALIZATION_NVP(spec_props_);
        ar & BOOST_SERIALIZATION_NVP(psf_type_);
        ar & BOOST_SERIALIZATION_NVP(psf_props_);
        ar & BOOST_SERIALIZATION_NVP(rows_);
        ar & BOOST_SERIALIZATION_NVP(cols_);
        ar & BOOST_SERIALIZATION_NVP(signal_);
        ar & BOOST_SERIALIZATION_NVP(noise_);
        ar & BOOST_SERIALIZATION_NVP(invcov_);
        ar & BOOST_SERIALIZATION_NVP(sky_);
        return;
      }
    
      std::string spec_type_;
      std::string psf_type_;
      boost::property_tree::ptree spec_props_;
      boost::property_tree::ptree psf_props_;
      size_t rows_;
      size_t cols_;
      vector_double signal_;
      vector_double noise_;
      vector_double invcov_;
      std::vector < bool > sky_;
    
  };

  BOOST_SERIALIZATION_SHARED_PTR(image_sim)

  image * image_sim_create ( boost::property_tree::ptree const & props );
  
  
}

#endif
