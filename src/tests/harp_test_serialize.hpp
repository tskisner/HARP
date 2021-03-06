/*
  High Performance Astrophysical Reconstruction and Processing (HARP)
  (c) 2014-2015, The Regents of the University of California, 
  through Lawrence Berkeley National Laboratory.  See top
  level LICENSE file for details.
*/

#ifndef HARP_TEST_SERIALIZE_HPP
#define HARP_TEST_SERIALIZE_HPP


namespace harp {

  class dwrapper {

    friend class boost::serialization::access;

    public :
      double val;

    private :

      template < class Archive >
      void serialize ( Archive & ar, const unsigned int version ) {
        ar & BOOST_SERIALIZATION_NVP(val);
        return;
      }
  };

  BOOST_SERIALIZATION_SHARED_PTR(dwrapper)

	class data_base : public boost::enable_shared_from_this < data_base > {

	  friend class boost::serialization::access;
	  
	  public :

	    data_base ( ) {
	      format_ = "";
	    }

	    data_base ( boost::property_tree::ptree const & props ) {
        format_ = props.get < std::string > ( "format", "" );
        props_ = props;
        dshr_.reset ( new dwrapper );
        (*dshr_).val = 543.2;
      }
	    
	    virtual ~data_base ( ) { }

	    virtual size_t get_sizet ( ) {
	      return 0;
	    }

	    virtual size_t get_double ( ) {
	      return 0.0;
	    }

	    std::string format ( ) { return format_; }
	    
	    static data_base * create ( boost::property_tree::ptree const & props );
	    
	    template < class T >
	    boost::shared_ptr < T > shared_ref ( ) {
	      return boost::dynamic_pointer_cast < T, data_base > ( shared_from_this ( ) );
	    }
	    
	    template < class T >
	    boost::weak_ptr < T > weak_ref ( ) {
	      boost::shared_ptr < T > temp = boost::dynamic_pointer_cast < T, data_base > ( shared_from_this() );
	      return boost::weak_ptr < T > ( temp );
	    }
	    
	  private :

	    template < class Archive >
	    void serialize ( Archive & ar, const unsigned int version ) {
	      ar & BOOST_SERIALIZATION_NVP(format_);
	      ar & BOOST_SERIALIZATION_NVP(props_);
        ar & BOOST_SERIALIZATION_NVP(dshr_);
	      return;
	    }
	  
	    std::string format_;
	    boost::property_tree::ptree props_;
	    boost::shared_ptr < dwrapper > dshr_;
	    
	};

  BOOST_SERIALIZATION_SHARED_PTR(data_base)

	BOOST_SERIALIZATION_ASSUME_ABSTRACT(data_base)

	typedef boost::shared_ptr < data_base > data_base_p;
	typedef boost::weak_ptr < data_base > data_base_wp;



	class data_one : public data_base {

	  friend class boost::serialization::access;
	    
	  public :

	    data_one ( ) : data_base () {
	      sizet_val_ = 0;
	      double_val_ = 0.0;
	    }

	    data_one ( boost::property_tree::ptree const & props ) : data_base ( props ) {
	      sizet_val_ = props.get < size_t > ( "sizet", 0 );
	      double_val_ = props.get < double > ( "double", 0.0 );
        doneshr_.reset ( new dwrapper );
        (*doneshr_).val = 789.0;
	    }
	    
	    ~data_one ( ) { }

	    size_t get_sizet ( ) {
	      return sizet_val_;
	    }

	    size_t get_double ( ) {
	      return double_val_;
	    }
	  
	  private :

	    template < class Archive >
	    void serialize ( Archive & ar, const unsigned int version ) {
	      ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(data_base);
	      ar & BOOST_SERIALIZATION_NVP(sizet_val_);
	      ar & BOOST_SERIALIZATION_NVP(double_val_);
        ar & BOOST_SERIALIZATION_NVP(doneshr_);
	      return;
	    }
	  
	    size_t sizet_val_;
	    double double_val_;
	    boost::shared_ptr < dwrapper > doneshr_;

	};

  BOOST_SERIALIZATION_SHARED_PTR(data_one)


}

  
#endif