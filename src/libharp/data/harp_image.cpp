/*
  High Performance Astrophysical Reconstruction and Processing (HARP)
  (c) 2014-2015, The Regents of the University of California, 
  through Lawrence Berkeley National Laboratory.  See top
  level LICENSE file for details.
*/

#include <harp_data_internal.hpp>


using namespace std;
using namespace harp;


harp::image::image ( std::string const & type, boost::property_tree::ptree const & props ) {
  props_ = props;
  type_ = type;
}


size_t harp::image::n_rows ( ) const {
  HARP_THROW( "fell through to virtual method" );
  return 0;
}


size_t harp::image::n_cols ( ) const {
  HARP_THROW( "fell through to virtual method" );
  return 0;
}


void harp::image::values ( vector_double & data ) const {
  HARP_THROW( "fell through to virtual method" );
  return;
}


void harp::image::inv_variance ( vector_double & invvar ) const {
  HARP_THROW( "fell through to virtual method" );
  return;
}


void harp::image::mask ( vector_mask & msk ) const {
  HARP_THROW( "fell through to virtual method" );
  return;
}


string harp::image::type ( ) const {
  return type_;
}


void harp::image::values ( matrix_double & data ) const {

  size_t imgrows = n_rows();
  size_t imgcols = n_cols();

  size_t nelem = imgrows * imgcols;

  data.resize ( imgrows, imgcols );

  vector_double tempdata ( nelem );

  values ( tempdata );

  for ( size_t i = 0; i < imgcols; ++i ) {
    for ( size_t j = 0; j < imgrows; ++j ) {
      data( j, i ) = tempdata[ i * imgrows + j ];
    }
  }

  return;
}


void harp::image::inv_variance ( matrix_double & invvar ) const {

  size_t imgrows = n_rows();
  size_t imgcols = n_cols();

  size_t nelem = imgrows * imgcols;

  invvar.resize ( imgrows, imgcols );

  vector_double tempvar ( nelem );

  inv_variance ( tempvar );

  for ( size_t i = 0; i < imgcols; ++i ) {
    for ( size_t j = 0; j < imgrows; ++j ) {
      invvar( j, i ) = tempvar[ i * imgrows + j ];
    }
  }

  return;
}


void harp::image::mask ( matrix_mask & msk ) const {

  size_t imgrows = n_rows();
  size_t imgcols = n_cols();

  size_t nelem = imgrows * imgcols;

  msk.resize ( imgrows, imgcols );

  vector_mask tempmask ( nelem );

  mask ( tempmask );

  for ( size_t i = 0; i < imgcols; ++i ) {
    for ( size_t j = 0; j < imgrows; ++j ) {
      msk( j, i ) = tempmask[ i * imgrows + j ];
    }
  }

  return;
}

