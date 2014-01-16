// @COPYRIGHT@

#include <harp_math_internal.hpp>

extern "C" {
  #include <sys/time.h>
}


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

