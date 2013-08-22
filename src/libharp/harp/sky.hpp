// @COPYRIGHT@

#ifndef HARP_EXTRACT_HPP
#define HARP_EXTRACT_HPP


namespace harp {

  void sky_design ( matrix_sparse const & AT, std::vector < bool > const & sky, matrix_sparse & skyAT, bool skysub = true );

  void sky_design2 ( matrix_sparse const & AT, std::vector < bool > const & sky, matrix_sparse & skyAT );

  void sky_subtract ( size_t nobj, matrix_dist const & Rf_orig, matrix_dist const & err_orig, matrix_dist & Rf, matrix_dist & err );

}

#endif

