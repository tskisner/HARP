// @COPYRIGHT@

#ifndef HARP_TARGET_HPP
#define HARP_TARGET_HPP


namespace harp {

  typedef enum {
    TARGET_SKY,
    TARGET_LRG,
    TARGET_ELG,
    TARGET_QSO,
    TARGET_STAR,
    TARGET_UNKNOWN
  } target_type;


  // FIXME: add information here about sky location, etc.

  class target {

    friend class boost::serialization::access;
    
    public :

      target ( ) {
        type_ = TARGET_UNKNOWN;
        name_ = "";
      }

      target ( target_type type, std::string name ) {
        type_ = type;
        name_ = name;
      }

      ~target ( ) { }

      // default copy / assignment operators are fine for now.  if 
      // class becomes more complicated, define them here.

      // accessor methods

      target_type type ( ) const {
        return type_;
      }

      std::string name ( ) const {
        return name_;
      }

      bool is_sky ( ) const {
        return ( type_ == TARGET_SKY );
      }
      
    private :

      template < class Archive >
      void serialize ( Archive & ar, const unsigned int version ) {
        ar & BOOST_SERIALIZATION_NVP(name_);
        ar & BOOST_SERIALIZATION_NVP(type_);
        return;
      }
    
      std::string name_;
      target_type type_;
      
  };


  BOOST_SERIALIZATION_SHARED_PTR(target)
  
  typedef boost::shared_ptr < harp::target > target_p;


}

#endif

