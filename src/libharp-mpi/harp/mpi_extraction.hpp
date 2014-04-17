// @COPYRIGHT@

#ifndef HARP_MPI_EXTRACT_HPP
#define HARP_MPI_EXTRACT_HPP


namespace harp {


  void mpi_sub_spec ( spec_slice_region const & full_region, spec_slice_region const & sub_region, mpi_matrix const & full_data, bool use_good_sub, mpi_matrix & sub_data );

  void mpi_accum_spec ( spec_slice_region const & sub_region, spec_slice_region const & full_region, mpi_matrix const & sub_data, bool use_good_sub, mpi_matrix & full_data );

  void mpi_noise_weighted_spec ( mpi_matrix_sparse const & AT, elem_matrix_local const & invnoise, vector_mask const & mask, elem_matrix_local const & img, mpi_matrix & z );

  void mpi_inverse_covariance ( mpi_matrix_sparse const & AT, elem_matrix_local const & invnoise, vector_mask const & mask, mpi_matrix & invcov );


  /*
  void mpi_resolution ( mpi_matrix & D, mpi_matrix & W, mpi_matrix & S, mpi_matrix & R );

  void mpi_extract ( mpi_matrix & D, mpi_matrix & W, mpi_matrix & S, mpi_matrix & z, mpi_matrix & Rf, mpi_matrix & f );

  void mpi_extract_slices ( mpi_spec_slice_p slice, mpi_psf_p design, elem_matrix_local const & img, elem_matrix_local const & img_inv_var, mpi_matrix const & truth, mpi_matrix & Rf, mpi_matrix & f, mpi_matrix & err, mpi_matrix & Rtruth, std::map < std::string, double > & profile, bool region_threads, bool lambda_mask, std::string const & status_prefix = "" );

  */

}

#endif

