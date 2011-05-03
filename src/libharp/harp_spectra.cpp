// @COPYRIGHT@

#include <harp_internal.hpp>


using namespace std;
using namespace harp;


harp::spectra::spectra ( string const & format, std::map < std::string, std::string > const & params ) {
  format_ = format;
  params_ = params;
}


void harp::spectra::cleanup ( ) {
  // nothing for now
  return;
}


string harp::spectra::format ( ) {
  return format_;
}


spectra * harp::spectra::clone ( ) {
  return create ( format_, params_ );
}


spectra * harp::spectra::create ( std::string const & format, std::map < std::string, std::string > const & params ) {
  
  #include "harp_spectra_formats.cpp"
  
  std::ostringstream o;
  o << "Cannot create spectra of unknown format (" << format << ")";
  MOAT_THROW( o.str().c_str() );

  return NULL;
  
}

