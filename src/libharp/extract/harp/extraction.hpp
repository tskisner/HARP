/*
  High Performance Astrophysical Reconstruction and Processing (HARP)
  (c) 2014-2015, The Regents of the University of California, 
  through Lawrence Berkeley National Laboratory.  See top
  level LICENSE file for details.
*/

#ifndef HARP_EXTRACTION_HPP
#define HARP_EXTRACTION_HPP

#include <ostream>


namespace harp {

  void sub_spec ( spec_slice_region const & full_region, spec_slice_region const & sub_region, vector_double const & full_data, bool use_good_sub, vector_double & sub_data );

  void accum_spec ( spec_slice_region const & sub_region, spec_slice_region const & full_region, vector_double const & sub_data, bool use_good_sub, vector_double & full_data );

  void noise_weighted_spec ( matrix_double_sparse const & AT, vector_double const & invnoise, vector_mask const & mask, vector_double const & img, vector_double & z );

  void inverse_covariance ( matrix_double_sparse const & AT, vector_double const & invnoise, vector_mask const & mask, matrix_double & invC );

  void resolution ( vector_double const & D, matrix_double const & W, vector_double & S, matrix_double & R );

  void accum_diag_resolution ( spec_slice_region const & region, matrix_double const & R, matrix_double_sparse & diag );

  void extract ( vector_double const & D, matrix_double const & W, vector_double const & S, vector_double const & z, vector_double & Rf, vector_double & f );

  void extract_slices ( spec_slice_p slice, psf_p design, vector_double const & img, vector_double const & img_inv_var, vector_double const & truth, vector_double & Rf, vector_double & f, vector_double & err, vector_double & Rtruth, std::map < std::string, double > & profile, bool region_threads, bool lambda_mask, std::string const & status_prefix = "" );

}

#endif

