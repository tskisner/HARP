// @COPYRIGHT@

#ifndef HARP_EXTRACT_HPP
#define HARP_EXTRACT_HPP


namespace harp {

  void sub_spec ( matrix_dist const & in, size_t total_nspec, size_t first_spec, size_t nspec, size_t first_lambda, size_t nlambda, matrix_dist & out );

  void accum_spec ( matrix_dist & full, size_t total_nspec, size_t first_spec, size_t nspec, size_t first_lambda, size_t nlambda, matrix_dist const & chunk );

  void spec_project ( matrix_sparse const & m, matrix_dist const & in, matrix_local & out );

  void noise_weighted_spec ( matrix_sparse const & psf, matrix_local const & invnoise, matrix_local const & img, matrix_dist & z );

  void inverse_covariance ( matrix_sparse const & psf, matrix_local const & invnoise, matrix_dist & invcov );

  void resolution ( matrix_dist & D, matrix_dist & W, matrix_dist & S, matrix_dist & R );

  void extract ( matrix_dist & D, matrix_dist & W, matrix_dist & S, matrix_dist & z, matrix_dist & Rf, matrix_dist & f );

  void sky_design ( matrix_sparse const & AT, std::vector < bool > const & sky, matrix_sparse & skyAT, bool skysub = true );

  void sky_subtract ( size_t nobj, matrix_dist const & Rf_orig, matrix_dist const & err_orig, matrix_dist & Rf, matrix_dist & err );

}

#endif

