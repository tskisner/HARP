/*
  High Performance Astrophysical Reconstruction and Processing (HARP)
  (c) 2014-2015, The Regents of the University of California, 
  through Lawrence Berkeley National Laboratory.  See top
  level LICENSE file for details.
*/

#ifndef HARP_TARGETS_HPP
#define HARP_TARGETS_HPP


namespace harp {

  typedef enum {
    OBJECT_CALIB,
    OBJECT_SKY,
    OBJECT_LRG,
    OBJECT_ELG,
    OBJECT_QSO,
    OBJECT_STAR,
    OBJECT_UNKNOWN
  } object_type;

  object_type object_str2type ( std::string const & in );

  std::string object_type2str ( object_type const & in );


  // FIXME: add information here about sky location, etc.

  class object : public boost::enable_shared_from_this < object > {

    friend class boost::serialization::access;
    
    public :

      object ( );

      object ( object_type type, std::string name );

      virtual ~object ( );

      // default copy / assignment operators are fine for now.  if 
      // class becomes more complicated, define them here.

      // accessor methods

      object_type type ( ) const;

      std::string name ( ) const;

      void set_type ( object_type t );

      void set_name ( std::string const & n );

      bool is_sky ( ) const;
      
    private :

      template < class Archive >
      void serialize ( Archive & ar, const unsigned int version ) {
        ar & BOOST_SERIALIZATION_NVP(name_);
        ar & BOOST_SERIALIZATION_NVP(type_);
        return;
      }
    
      std::string name_;
      object_type type_;
      
  };

  BOOST_SERIALIZATION_ASSUME_ABSTRACT(object)

  BOOST_SERIALIZATION_SHARED_PTR(object)
  
  typedef boost::shared_ptr < harp::object > object_p;


  class targets : public boost::enable_shared_from_this < targets > {

    friend class boost::serialization::access;
    
    public :

      targets ( );

      targets ( std::string const & type, boost::property_tree::ptree const & props );
      
      virtual ~targets ( );

      virtual size_t n_objects ( ) const;

      virtual std::vector < object_p > objects ( ) const;
      
      std::string type ( ) const;
      
      template < class T >
      boost::shared_ptr < T > shared_ref ( ) {
        return boost::dynamic_pointer_cast < T, targets > ( shared_from_this ( ) );
      }
      
      template < class T >
      boost::weak_ptr < T > weak_ref ( ) {
        boost::shared_ptr < T > temp = boost::dynamic_pointer_cast < T, targets > ( shared_from_this() );
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

  BOOST_SERIALIZATION_ASSUME_ABSTRACT(targets)

  BOOST_SERIALIZATION_SHARED_PTR(targets)
  
  typedef boost::shared_ptr < harp::targets > targets_p;
  typedef boost::weak_ptr < harp::targets > targets_wp;


}

#endif

