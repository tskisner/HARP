// @COPYRIGHT@

#include <harp_internal.hpp>


using namespace std;
using namespace harp;


harp::exception::exception ( char const * msg, char const * file, int line ) {
  int ret;
  ret = snprintf ( msg_, BIGSTRLEN, "Exception at line %d of file %s:  %s", line, file, msg );
}


harp::exception::~exception ( ) throw() { }


const char* harp::exception::what() const throw() { 
  return msg_;
}


// MPI return value checking

void harp::mpi_check ( MPI_Comm comm, int status ) {
  if ( status != MPI_SUCCESS ) {
    int myp;
    MPI_Comm_rank ( comm, &myp );
    
    char msg[ MPI_MAX_ERROR_STRING ];
    
    int len;
    int ret = MPI_Error_string( status, msg, &len );
    
    std::ostringstream o;
    o << "MPI error (process " << myp << "): " << msg;
    HARP_THROW( o.str().c_str() );
  }
  return;
}


// utilities


void harp::dist_matrix_zero ( matrix_dist & mat ) {
  size_t n = mat.AllocatedMemory();
  double * data = mat.LocalBuffer();
  for ( size_t i = 0; i < n; ++i ) {
    data[i] = 0.0;
  }
  return;
}


void harp::local_matrix_zero ( matrix_local & mat ) {
  size_t n = mat.MemorySize();
  double * data = mat.Buffer();
  for ( size_t i = 0; i < n; ++i ) {
    data[i] = 0.0;
  }
  return;
}


