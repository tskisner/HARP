// @COPYRIGHT@

#ifndef HARP_EXTRACT_HPP
#define HARP_EXTRACT_HPP


namespace harp {

  void spec_project ( matrix_sparse const & m, matrix_dist const & in, matrix_local & out );

  void noise_weighted_spec ( matrix_sparse const & psf, matrix_local const & invnoise, matrix_local const & img, matrix_dist & z );

  void inverse_covariance ( matrix_sparse const & psf, matrix_local const & invnoise, matrix_dist & invcov );

  void resolution ( matrix_dist & D, matrix_dist & W, matrix_dist & S, matrix_dist & R );

  void extract ( matrix_dist & D, matrix_dist & W, matrix_dist & S, matrix_dist & z, matrix_dist & f );


}

#endif

