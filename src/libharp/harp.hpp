// @COPYRIGHT@

#ifndef HARP_HPP
#define HARP_HPP

#include <string>
#include <memory>
#include <vector>
#include <map>
#include <set>
#include <sstream>


#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/storage.hpp>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree_serialization.hpp>

#include <boost/smart_ptr/shared_ptr.hpp>
#include <boost/smart_ptr/weak_ptr.hpp>
#include <boost/enable_shared_from_this.hpp>

#include <boost/serialization/string.hpp>
#include <boost/serialization/version.hpp>
#include <boost/serialization/split_member.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/nvp.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/list.hpp>
#include <boost/serialization/export.hpp>
#include <boost/serialization/assume_abstract.hpp>


#include <boost/foreach.hpp>


namespace harp {

  static const int STRLEN = 256;
  static const int BIGSTRLEN = 1024;

  typedef enum {
    EIG_NONE,
    EIG_SQRT,
    EIG_INVSQRT
  } eigen_op;

  // Exception handling

  class exception : public std::exception {
    
    public:
      exception ( char const * msg, char const * file, int line );
      ~exception ( ) throw ();
      char const * what() const throw();
    
    private:  
      char msg_[BIGSTRLEN];
      
  };  

  typedef void (*HARP_EXCEPTION_HANDLER) ( harp::exception & e );

  #define HARP_THROW(msg) \
  throw harp::exception ( msg, __FILE__, __LINE__ )

  #define HARP_TRY \
  try {

  #define HARP_CATCH \
  } catch ( harp::exception & e ) { \
    std::cerr << e.what() << std::endl; \
    throw; \
  }

  #define HARP_CATCH_CUSTOM(handler) \
  } catch ( harp::exception & e ) { \
    (*handler) ( e ); \
  }

  // utilities

  std::string ptree_quote ( const std::string & s );
  void ptree_print ( const boost::property_tree::ptree & pt, int level );
  void ptree_print ( const boost::property_tree::ptree & pt );
  
}

#include <harp/matrix.hpp>

#include <harp/fits.hpp>

#include <harp/image.hpp>
#include <harp/spec.hpp>
#include <harp/psf.hpp>

#include <harp/extraction.hpp>

#include <harp/sky.hpp>

#endif



