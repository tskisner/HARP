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

std::string harp::ptree_quote ( const std::string & s ) {
  return "\"" + s + "\"";
}

void harp::ptree_print ( const boost::property_tree::ptree & pt, int level ) {
  const std::string sep ( 2 * level, ' ' );
  BOOST_FOREACH ( const boost::property_tree::ptree::value_type & v, pt ) {
    std::cerr << sep << ptree_quote ( v.first ) << " : " << ptree_quote ( v.second.data() ) << "\n";
    ptree_print ( v.second, level + 1 );
  }
  return;
}

void harp::ptree_print ( const boost::property_tree::ptree & pt ) {
  ptree_print ( pt, 0 );
  return;
}


