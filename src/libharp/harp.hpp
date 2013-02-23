// @COPYRIGHT@

#ifndef HARP_HPP
#define HARP_HPP

#include <string>
#include <memory>
#include <vector>
#include <map>
#include <sstream>

#include <boost/smart_ptr/shared_ptr.hpp>
#include <boost/smart_ptr/weak_ptr.hpp>
#include <boost/enable_shared_from_this.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include <clique.hpp>


namespace harp {

  static const int STRLEN = 256;
  static const int BIGSTRLEN = 1024;

  typedef enum {
    EIG_NONE,
    EIG_SQRT,
    EIG_INVSQRT
  } eigen_op;

  // typedefs for linear algebra

  typedef elem::Matrix < double > matrix_local;

  typedef elem::DistMatrix < double, elem::MC, elem::MR > matrix_dist;

  typedef cliq::DistSparseMatrix < double > matrix_sparse;

  typedef cliq::DistVector < double > vector_dist;

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

  void mpi_check ( MPI_Comm comm, int status );

  // utilities

  void dist_matrix_zero ( matrix_dist & mat );

  void local_matrix_zero ( matrix_local & mat );
  
}


#include <harp/fits.hpp>
#include <harp/matrix.hpp>

#include <harp/image.hpp>
#include <harp/spec.hpp>
#include <harp/psf.hpp>

#include <harp/extract.hpp>
#include <harp/metadata.hpp>

#endif



