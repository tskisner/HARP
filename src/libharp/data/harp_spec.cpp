/*
  High Performance Astrophysical Reconstruction and Processing (HARP)
  (c) 2014-2015, The Regents of the University of California, 
  through Lawrence Berkeley National Laboratory.  See top
  level LICENSE file for details.
*/

#include <harp_data_internal.hpp>


using namespace std;
using namespace harp;


harp::spec::spec ( std::string const & type, boost::property_tree::ptree const & props ) {
  props_ = props;
  type_ = type;
}


string harp::spec::type ( ) const {
  return type_;
}


void harp::spec::values ( matrix_double & data ) const {

  size_t nspec = n_spec();
  size_t nlambda = n_lambda();

  size_t nelem = nspec * nlambda;

  data.resize ( nlambda, nspec );

  vector_double tempdata ( nelem );

  values ( tempdata );

  for ( size_t i = 0; i < nspec; ++i ) {
    for ( size_t j = 0; j < nlambda; ++j ) {
      data( j, i ) = tempdata[ i * nlambda + j ];
    }
  }

  return;
}

void harp::spec::inv_variance ( matrix_double & data ) const {

  size_t nspec = n_spec();
  size_t nlambda = n_lambda();

  size_t nelem = nspec * nlambda;

  data.resize ( nlambda, nspec );

  vector_double tempdata ( nelem );

  inv_variance ( tempdata );

  for ( size_t i = 0; i < nspec; ++i ) {
    for ( size_t j = 0; j < nlambda; ++j ) {
      data( j, i ) = tempdata[ i * nlambda + j ];
    }
  }

  return;
}

