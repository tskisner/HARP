// @COPYRIGHT@

#ifndef HARP_MATH_HPP
#define HARP_MATH_HPP

#include <string>
#include <memory>
#include <vector>
#include <map>
#include <set>
#include <sstream>

#include <boost/smart_ptr/shared_ptr.hpp>
#include <boost/smart_ptr/weak_ptr.hpp>
#include <boost/enable_shared_from_this.hpp>

#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/portable_oarchive.hpp>
#include <boost/archive/portable_iarchive.hpp>

#include <boost/serialization/string.hpp>
#include <boost/serialization/version.hpp>
#include <boost/serialization/split_member.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/nvp.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/export.hpp>
#include <boost/serialization/assume_abstract.hpp>

#include <boost/serialization/list.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/map.hpp>


namespace harp {

  static const int STRLEN = 256;
  static const int BIGSTRLEN = 1024;

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

  // if we have been built with python support, enable functions to
  // repackage python exceptions into thrown C++ exceptions.

  std::string parse_python_exception();

  // helper wrapper to gettimeofday

  double wtime ( );

  // 1D distribution across OpenMP thread

  void omp_dist_1D ( int n, int & rank, int & nthreads, int & myn, int & offset );

}

#include <harp/linalg.hpp>


#endif
