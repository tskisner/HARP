// @COPYRIGHT@

#include <harp.hpp>


using namespace std;
using namespace harp;


void harp::fits_check ( int status ) {
  if ( status ) {
    char msg[ FLEN_ERRMSG ];
    fits_get_errstatus ( status, msg );

    std::ostringstream o;
    o << "cfitsio library error: " << msg;
    MOAT_THROW( o.str().c_str() );
  }
  return;
}