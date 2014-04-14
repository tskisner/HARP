// @COPYRIGHT@

#ifndef HARP_MPI_EXTRACT_HPP
#define HARP_MPI_EXTRACT_HPP


namespace harp {


  void mpi_sub_spec ( spec_slice_region const & full_region, spec_slice_region const & sub_region, mpi_matrix const & full_data, bool use_good_sub, mpi_matrix & sub_data );

  void mpi_accum_spec ( spec_slice_region const & sub_region, spec_slice_region const & full_region, mpi_matrix const & sub_data, bool use_good_sub, mpi_matrix & full_data );


  /*

  void sub_spec ( matrix_dist const & in, size_t total_nspec, size_t first_spec, size_t nspec, size_t first_lambda, size_t nlambda, matrix_dist & out );

  void accum_spec ( matrix_dist & full, size_t total_nspec, size_t first_spec, size_t nspec, size_t first_lambda, size_t nlambda, matrix_dist const & chunk );

  void noise_weighted_spec ( matrix_sparse const & psf, matrix_local const & invnoise, matrix_local const & img, matrix_dist & z );

  void inverse_covariance ( matrix_sparse const & psf, matrix_local const & invnoise, matrix_dist & invcov );

  void resolution ( matrix_dist & D, matrix_dist & W, matrix_dist & S, matrix_dist & R );

  void extract ( matrix_dist & D, matrix_dist & W, matrix_dist & S, matrix_dist & z, matrix_dist & Rf, matrix_dist & f );

  */

}

#endif

