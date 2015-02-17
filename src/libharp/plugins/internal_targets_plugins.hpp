/*
  High Performance Astrophysical Reconstruction and Processing (HARP)
  (c) 2014-2015, The Regents of the University of California, 
  through Lawrence Berkeley National Laboratory.  See top
  level LICENSE file for details.
*/

#ifndef HARP_INTERNAL_TARGETS_PLUGINS_HPP
#define HARP_INTERNAL_TARGETS_PLUGINS_HPP

namespace harp {
    

  class targets_fits : public targets {

    friend class boost::serialization::access;
    
    public :

      targets_fits ( );

      targets_fits ( boost::property_tree::ptree const & props );
      
      ~targets_fits ( );

      static void write ( std::string const & path, int hdu, std::vector < object_p > const & objects );

      // overloaded virtual methods from base class

      size_t n_objects ( ) const;

      std::vector < object_p > objects ( ) const;
    
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

      targets_sim ( );

      targets_sim ( boost::property_tree::ptree const & props );
      
      ~targets_sim ( );

      // overloaded virtual methods from base class

      size_t n_objects ( ) const;

      std::vector < object_p > objects ( ) const;
    
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


  class object_desi : public object {

    friend class boost::serialization::access;
    
    public :

      object_desi ( );

      object_desi ( object_type type,
        std::string in_targetcat,
        int64_t in_targetid,
        int64_t in_targetmask,
        float * in_mag,
        std::string in_filter,
        int64_t in_spectroid,
        int64_t in_positioner,
        int32_t in_fiber,
        float in_lambdaref,
        double in_ra_target,
        double in_dec_target,
        double in_ra_obs,
        double in_dec_obs,
        double in_x_target,
        double in_y_target,
        double in_x_fvcobs,
        double in_y_fvcobs,
        float in_x_fvcerr,
        float in_y_fvcerr
      );

      std::string targetcat;
      int64_t targetid;
      int64_t targetmask;
      float mag[5];
      std::string filter;
      int64_t spectroid;
      int64_t positioner;
      int32_t fiber;
      float lambdaref;
      double ra_target;
      double dec_target;
      double ra_obs;
      double dec_obs;
      double x_target;
      double y_target;
      double x_fvcobs;
      double y_fvcobs;
      float x_fvcerr;
      float y_fvcerr;

    private :

      template < class Archive >
      void serialize ( Archive & ar, const unsigned int version ) {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(object);
        ar & BOOST_SERIALIZATION_NVP(targetcat);
        ar & BOOST_SERIALIZATION_NVP(targetid);
        ar & BOOST_SERIALIZATION_NVP(targetmask);
        ar & BOOST_SERIALIZATION_NVP(mag);
        ar & BOOST_SERIALIZATION_NVP(filter);
        ar & BOOST_SERIALIZATION_NVP(spectroid);
        ar & BOOST_SERIALIZATION_NVP(positioner);
        ar & BOOST_SERIALIZATION_NVP(fiber);
        ar & BOOST_SERIALIZATION_NVP(lambdaref);
        ar & BOOST_SERIALIZATION_NVP(ra_target);
        ar & BOOST_SERIALIZATION_NVP(dec_target);
        ar & BOOST_SERIALIZATION_NVP(ra_obs);
        ar & BOOST_SERIALIZATION_NVP(dec_obs);
        ar & BOOST_SERIALIZATION_NVP(x_target);
        ar & BOOST_SERIALIZATION_NVP(y_target);
        ar & BOOST_SERIALIZATION_NVP(x_fvcobs);
        ar & BOOST_SERIALIZATION_NVP(y_fvcobs);
        ar & BOOST_SERIALIZATION_NVP(x_fvcerr);
        ar & BOOST_SERIALIZATION_NVP(y_fvcerr);
        return;
      }

  };

  BOOST_SERIALIZATION_SHARED_PTR(object_desi)
  
  typedef boost::shared_ptr < harp::object_desi > object_desi_p;


  class targets_desi : public targets {

    friend class boost::serialization::access;
    
    public :

      targets_desi ( );

      targets_desi ( boost::property_tree::ptree const & props );
      
      ~targets_desi ( );

      static void write ( std::string const & path, boost::property_tree::ptree const & meta, std::vector < object_desi_p > const & objects );

      // overloaded virtual methods from base class

      size_t n_objects ( ) const;

      std::vector < object_p > objects ( ) const;

      // desi-specific

      boost::property_tree::ptree meta () const;

      std::vector < object_desi_p > desi_objects ( ) const;

      static std::vector < std::string > colnames ();

      static object_type desi2harp ( std::string const & name );

      static std::string harp2desi ( object_type type );

      int tileid;
      float telera;
      float teledec;
      int expid;
      std::string night;
      std::string vdmodel;
      std::string voptics;
      std::string vfibvcam;
      float hexpdrot;
      std::string dateobs;
    
    private :

      template < class Archive >
      void serialize ( Archive & ar, const unsigned int version ) {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(targets);
        ar & BOOST_SERIALIZATION_NVP(nobjects_);
        ar & BOOST_SERIALIZATION_NVP(path_);
        ar & BOOST_SERIALIZATION_NVP(hdu_);
        ar & BOOST_SERIALIZATION_NVP(objects_);
        ar & BOOST_SERIALIZATION_NVP(meta_);
        ar & BOOST_SERIALIZATION_NVP(tileid);
        ar & BOOST_SERIALIZATION_NVP(telera);
        ar & BOOST_SERIALIZATION_NVP(teledec);
        ar & BOOST_SERIALIZATION_NVP(expid);
        ar & BOOST_SERIALIZATION_NVP(night);
        ar & BOOST_SERIALIZATION_NVP(vdmodel);
        ar & BOOST_SERIALIZATION_NVP(voptics);
        ar & BOOST_SERIALIZATION_NVP(vfibvcam);
        ar & BOOST_SERIALIZATION_NVP(hexpdrot);
        ar & BOOST_SERIALIZATION_NVP(dateobs);
        return;
      }
    
      boost::property_tree::ptree meta_;
      size_t nobjects_;
      std::string path_;
      int hdu_;
      std::vector < object_desi_p > objects_;
    
  };

  BOOST_SERIALIZATION_SHARED_PTR(targets_desi)

  targets * targets_desi_create ( boost::property_tree::ptree const & props );


  
}

#endif
