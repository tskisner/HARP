// @COPYRIGHT@

#include <harp_internal.hpp>


using namespace std;
using namespace harp;

static const char * format_boss = "boss";

static const char * boss_psf_key_path = "path";

static const char * boss_psf_key_name = "PSFPARAM";

static const char * boss_psf_hdu_x = "X";
static const char * boss_psf_hdu_y = "Y";

static const char * boss_psf_key_corr = "corr";



harp::psf_boss::psf_boss ( std::map < std::string, std::string > const & params ) : psf ( format_boss, params ) {
  
}


harp::psf_boss::~psf_boss ( ) {
  
}


void harp::psf_boss::lambda ( size_t specnum, data_vec & data ) {
  
  return;
}


void harp::psf_boss::extent ( size_t firstspec, size_t lastspec, size_t firstbin, size_t lastbin, size_t & firstX, size_t & firstY, size_t & lastX, size_t & lastY ) {
  
  return;
}


void harp::psf_boss::projection ( size_t firstspec, size_t lastspec, size_t firstbin, size_t lastbin, size_t firstX, size_t lastX, size_t firstY, size_t lastY, comp_rowmat & data ) {
  
  return;
}


