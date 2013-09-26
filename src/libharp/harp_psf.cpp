// @COPYRIGHT@

#include <harp_internal.hpp>


using namespace std;
using namespace harp;


harp::psf::psf ( boost::property_tree::ptree const & props ) {
  props_ = props;
  format_ = props.get < string > ( "format" );
}


string harp::psf::format ( ) {
  return format_;
}


psf * harp::psf::clone ( ) {
  return create ( props_ );
}


psf * harp::psf::create ( boost::property_tree::ptree const & props ) {
  
  string format = props.get < string > ( "format" );

  /*
  if ( format == "gauss" ) {
    return static_cast < psf * > ( new psf_gauss ( props ) );
  }
  */
  
  std::ostringstream o;
  o << "Cannot create psf of unknown format (" << format << ")";
  HARP_THROW( o.str().c_str() );

  return NULL;
  
}


size_t harp::psf::total_bins ( std::map < size_t, std::set < size_t > > const & speclambda ) {

  size_t total = 0;

  for ( std::map < size_t, std::set < size_t > > :: const_iterator itspec = speclambda.begin(); itspec != speclambda.end(); ++itspec ) {
    for ( std::set < size_t > :: const_iterator itlambda = itspec->second.begin(); itlambda != itspec->second.end(); ++itlambda ) {
      ++total;
    }
  }

  return total;
}


void harp::psf::project ( std::map < size_t, std::set < size_t > > const & speclambda, matrix_double & A ) {

  size_t total = total_bins ( speclambda );

  size_t pixrows = img_rows();
  size_t pixcols = img_cols();
  size_t pixtotal = pixrows * pixcols;

  // resize output to correct dimensions

  A.resize ( total, pixtotal );
  A.clear();

  // iterate over spectral bins and populate the matrix elements

  matrix_double patch;
  size_t xoff;
  size_t yoff;

  size_t col = 0;
  size_t row;

  for ( std::map < size_t, std::set < size_t > > :: const_iterator itspec = speclambda.begin(); itspec != speclambda.end(); ++itspec ) {

    for ( std::set < size_t > :: const_iterator itlambda = itspec->second.begin(); itlambda != itspec->second.end(); ++itlambda ) {

      response ( itspec->first, (*itlambda), xoff, yoff, patch );

      for ( size_t patch_col = 0; patch_col < patch.size1(); ++patch_col ) {

        if ( xoff + patch_col < pixcols ) {
          // this column is within the image dimensions

          for ( size_t patch_row = 0; patch_row < patch.size2(); ++patch_row ) {

            if ( yoff + patch_row < pixrows ) {
              // this row is within the image dimensions

              row = ( xoff + patch_col ) * pixrows + yoff + patch_row;

              A ( col, row ) = patch ( patch_col, patch_row );

            }

          }

        }

      }

      ++col;

    }

  }

  return;
}


void harp::psf::project_transpose ( std::map < size_t, std::set < size_t > > const & speclambda, matrix_double_sparse & AT ) {

  size_t total = total_bins ( speclambda );

  size_t pixrows = img_rows();
  size_t pixcols = img_cols();
  size_t pixtotal = pixrows * pixcols;

  // resize output to correct dimensions

  AT.resize ( total, pixtotal );
  AT.clear();

  // iterate over spectral bins and populate the matrix elements

  matrix_double patch;
  size_t xoff;
  size_t yoff;

  size_t row = 0;
  size_t col;

  for ( std::map < size_t, std::set < size_t > > :: const_iterator itspec = speclambda.begin(); itspec != speclambda.end(); ++itspec ) {

    for ( std::set < size_t > :: const_iterator itlambda = itspec->second.begin(); itlambda != itspec->second.end(); ++itlambda ) {

      response ( itspec->first, (*itlambda), xoff, yoff, patch );

      for ( size_t patch_col = 0; patch_col < patch.size1(); ++patch_col ) {

        if ( xoff + patch_col < pixcols ) {
          // this column is within the image dimensions

          for ( size_t patch_row = 0; patch_row < patch.size2(); ++patch_row ) {

            if ( yoff + patch_row < pixrows ) {
              // this row is within the image dimensions

              col = ( xoff + patch_col ) * pixrows + yoff + patch_row;

              AT ( row, col ) = patch ( patch_col, patch_row );

            }

          }

        }

      }

      ++row;

    }

  }
  
  return;
}


void harp::psf::project ( matrix_double & A ) {

  std::map < size_t, std::set < size_t > > speclambda;

  for ( size_t spec = 0; spec < n_spec(); ++spec ) {
    for ( size_t lambda = 0; lambda < n_lambda(); ++lambda ) {
      speclambda[ spec ].insert ( lambda );
    }
  }

  project ( speclambda, A );

  return;
}


void harp::psf::project_transpose ( matrix_double_sparse & AT ) {

  std::map < size_t, std::set < size_t > > speclambda;

  for ( size_t spec = 0; spec < n_spec(); ++spec ) {
    for ( size_t lambda = 0; lambda < n_lambda(); ++lambda ) {
      speclambda[ spec ].insert ( lambda );
    }
  }

  project_transpose ( speclambda, AT );
  
  return;
}


