// @COPYRIGHT@

#include <harp_math_internal.hpp>

#include <cstdio>

extern "C" {
  #include <sys/time.h>
}

#ifdef _OPENMP
#  include <omp.h>
#endif

#include <cstdio>


using namespace std;
using namespace harp;


double harp::wtime ( ) {
  struct timeval now;
  int ret = gettimeofday ( &now, NULL );

  if ( ret != 0 ) {
    HARP_THROW( "gettimeofday() failed" );
  }

  double nowtime = static_cast < double > ( now.tv_sec ) + 1.0e-6 * static_cast < double > ( now.tv_usec );

  return nowtime;
}


void harp::omp_dist_1D ( int n, int & rank, int & nthreads, int & myn, int & offset ) {
#ifdef _OPENMP  
  rank = omp_get_thread_num();
  nthreads = omp_get_num_threads();
#else
  rank = 0;
  nthreads = 1;
#endif
  myn = (int) ( n / nthreads );
  
  int leftover = n % nthreads;
  
  if ( rank < leftover ) {
    ++myn;
    offset = myn * rank;
  } else {
    offset = ( (myn + 1) * leftover ) + ( myn * (rank - leftover) );
  }
  
  return;
}


harp::exception::exception ( char const * msg, char const * file, int line ) {
  int ret;
  ret = snprintf ( msg_, BIGSTRLEN, "Exception at line %d of file %s:  %s", line, file, msg );
}


harp::exception::~exception ( ) throw() { }


const char* harp::exception::what() const throw() { 
  return msg_;
}




