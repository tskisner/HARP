// @COPYRIGHT@

#include <harp_data_internal.hpp>


using namespace std;
using namespace harp;


harp::psf::psf ( std::string const & name, boost::property_tree::ptree const & props ) {
  props_ = props;
  name_ = name;
}


string harp::psf::name ( ) const {
  return name_;
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


void harp::psf::project ( std::map < size_t, std::set < size_t > > const & speclambda, matrix_double & A ) const {

  size_t total = total_bins ( speclambda );

  size_t pixrows = img_rows();
  size_t pixcols = img_cols();
  size_t pixtotal = pixrows * pixcols;

  // resize output to correct dimensions

  A.resize ( pixtotal, total, false );
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

      for ( size_t patch_col = 0; patch_col < patch.size2(); ++patch_col ) {

        if ( xoff + patch_col < pixcols ) {
          // this column is within the image dimensions

          for ( size_t patch_row = 0; patch_row < patch.size1(); ++patch_row ) {

            if ( yoff + patch_row < pixrows ) {
              // this row is within the image dimensions

              row = ( xoff + patch_col ) * pixrows + yoff + patch_row;

              A ( row, col ) = patch ( patch_row, patch_col );

            }

          }

        }

      }

      ++col;

    }

  }

  return;
}


void harp::psf::project_transpose ( std::map < size_t, std::set < size_t > > const & speclambda, matrix_double_sparse & AT ) const {

  size_t total = total_bins ( speclambda );

  size_t pixrows = img_rows();
  size_t pixcols = img_cols();
  size_t pixtotal = pixrows * pixcols;

  // resize output to correct dimensions

  AT.resize ( total, pixtotal, false );
  AT.clear();

  // we use the estimated number of non-zeros per row times the number 
  // of rows to reserve space in the sparse matrix.

  size_t nnz_estimate = response_nnz_estimate();
  nnz_estimate *= total;

  AT.reserve ( nnz_estimate, false );

  // iterate over spectral bins and populate the matrix elements

  matrix_double patch;
  size_t xoff;
  size_t yoff;

  size_t row = 0;
  size_t col;

  for ( std::map < size_t, std::set < size_t > > :: const_iterator itspec = speclambda.begin(); itspec != speclambda.end(); ++itspec ) {

    for ( std::set < size_t > :: const_iterator itlambda = itspec->second.begin(); itlambda != itspec->second.end(); ++itlambda ) {

      response ( itspec->first, (*itlambda), xoff, yoff, patch );

      for ( size_t patch_col = 0; patch_col < patch.size2(); ++patch_col ) {

        if ( xoff + patch_col < pixcols ) {
          // this column is within the image dimensions

          for ( size_t patch_row = 0; patch_row < patch.size1(); ++patch_row ) {

            if ( yoff + patch_row < pixrows ) {
              // this row is within the image dimensions

              col = ( xoff + patch_col ) * pixrows + yoff + patch_row;

              AT ( row, col ) = patch ( patch_row, patch_col );

            }

          }

        }

      }

      ++row;

    }

  }
  
  return;
}


void harp::psf::project ( matrix_double & A ) const {

  std::map < size_t, std::set < size_t > > speclambda;

  for ( size_t spec = 0; spec < n_spec(); ++spec ) {
    for ( size_t lambda = 0; lambda < n_lambda(); ++lambda ) {
      speclambda[ spec ].insert ( lambda );
    }
  }

  project ( speclambda, A );

  return;
}


void harp::psf::project_transpose ( matrix_double_sparse & AT ) const {

  std::map < size_t, std::set < size_t > > speclambda;

  for ( size_t spec = 0; spec < n_spec(); ++spec ) {
    for ( size_t lambda = 0; lambda < n_lambda(); ++lambda ) {
      speclambda[ spec ].insert ( lambda );
    }
  }

  project_transpose ( speclambda, AT );
  
  return;
}


