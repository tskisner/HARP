// @COPYRIGHT@

#include <harp_internal.hpp>


using namespace std;
using namespace harp;

static const char * format_boss_specter = "boss_specter";

static const char * boss_specter_key_path = "path";
static const char * boss_specter_key_hdu = "hdu";
static const char * boss_specter_key_nspec = "nspec";
static const char * boss_specter_key_specsize = "specsize";


harp::spec_boss_specter::spec_boss_specter ( boost::property_tree::ptree const & props ) : spec ( format_boss_specter, props ) {
  
  
}


harp::spec_boss_specter::~spec_boss_specter ( ) {
  
  cleanup();
  
}


boost::property_tree::ptree harp::spec_boss_specter::serialize ( ) {
  boost::property_tree::ptree ret;


  return ret;
}


void harp::spec_boss_specter::read ( matrix_dist & data, std::vector < double > & lambda, std::vector < bool > & sky ) {

  
  return;
}


void harp::spec_boss_specter::write ( std::string const & path, matrix_dist & data, std::vector < double > const & lambda, std::vector < bool > const & sky ) {

  
  return;
}



