// @COPYRIGHT@

#include <harp_math_internal.hpp>


using namespace std;
using namespace harp;


void harp::check_column_major ( boost::numeric::ublas::column_major_tag ) {
  return;
}


void harp::check_column_major ( boost::numeric::ublas::row_major_tag ) {
  HARP_THROW( "Matrices must be column-major, in order to maintain compatibility with LAPACK" );
  return;
}

