// @COPYRIGHT@

#include <harp.hpp>
#include <harp_example.hpp>

using namespace std;
using namespace harp;


harp::spec_example::spec_example ( boost::property_tree::ptree const & props ) : spec ( props ) {
  
  nspec_ = 1;

  nlambda_ = 1;
  
}


harp::spec_example::~spec_example ( ) {
  
}


boost::property_tree::ptree harp::spec_example::metadata ( ) const {

  return boost::property_tree::ptree();
}


void harp::spec_example::values ( vector_double & data ) const {
  
  return;
}


void harp::spec_example::lambda ( vector_double & lambda_vals ) const {
  
  return;
}


void harp::spec_example::targets ( std::vector < target > & target_list ) const {

  return;
}





