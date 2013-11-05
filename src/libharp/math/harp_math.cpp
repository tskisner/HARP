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

