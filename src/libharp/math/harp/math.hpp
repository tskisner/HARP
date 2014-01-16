// @COPYRIGHT@

#ifndef HARP_MATH_HPP
#define HARP_MATH_HPP

#include <string>
#include <memory>
#include <vector>
#include <map>
#include <set>
#include <sstream>


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

  // helper wrapper to gettimeofday

  double wtime ( );

  // 1D distribution across OpenMP thread

  void omp_dist_1D ( int n, int & rank, int & nthreads, int & myn, int & offset );

}

#include <harp/linalg.hpp>


#endif
