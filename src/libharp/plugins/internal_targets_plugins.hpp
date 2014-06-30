// @COPYRIGHT@

#ifndef HARP_INTERNAL_TARGETS_PLUGINS_HPP
#define HARP_INTERNAL_TARGETS_PLUGINS_HPP

namespace harp {
    

  class targets_fits : public targets {

    friend class boost::serialization::access;
    
    public :

      targets_fits ( ) : targets () {
        nobjects_ = 0;
        objects_.clear();
        hdu_ = 1;
        path_ = "";
      }

      targets_fits ( boost::property_tree::ptree const & props );
      
      ~targets_fits ( );

      static void write ( std::string const & path, int hdu, std::vector < object_p > const & objects );

      // overloaded virtual methods from base class

      virtual size_t n_objects ( ) const { return nobjects_; }

      virtual std::vector < object_p > objects ( ) const { return objects_; }
    
    private :

      template < class Archive >
      void serialize ( Archive & ar, const unsigned int version ) {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(targets);
        ar & BOOST_SERIALIZATION_NVP(nobjects_);
        ar & BOOST_SERIALIZATION_NVP(path_);
        ar & BOOST_SERIALIZATION_NVP(hdu_);
        ar & BOOST_SERIALIZATION_NVP(objects_);
        return;
      }
    
      size_t nobjects_;
      std::string path_;
      int hdu_;
      std::vector < object_p > objects_;
    
  };

  BOOST_SERIALIZATION_SHARED_PTR(targets_fits)

  targets * targets_fits_create ( boost::property_tree::ptree const & props );



  class targets_sim : public targets {

    friend class boost::serialization::access;
    
    public :

      targets_sim ( ) : targets () {
        nobjects_ = 0;
        objects_.clear();
        skymod_ = 0;
      }

      targets_sim ( boost::property_tree::ptree const & props );
      
      ~targets_sim ( );

      // overloaded virtual methods from base class

      virtual size_t n_objects ( ) const { return nobjects_; }

      virtual std::vector < object_p > objects ( ) const { return objects_; }
    
    private :

      template < class Archive >
      void serialize ( Archive & ar, const unsigned int version ) {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(targets);
        ar & BOOST_SERIALIZATION_NVP(nobjects_);
        ar & BOOST_SERIALIZATION_NVP(skymod_);
        ar & BOOST_SERIALIZATION_NVP(objects_);
        return;
      }
    
      size_t nobjects_;
      size_t skymod_;
      std::vector < object_p > objects_;
    
  };

  BOOST_SERIALIZATION_SHARED_PTR(targets_sim)

  targets * targets_sim_create ( boost::property_tree::ptree const & props );

  
}

#endif
