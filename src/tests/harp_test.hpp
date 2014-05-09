// @COPYRIGHT@

#ifndef HARP_TEST_HPP
#define HARP_TEST_HPP

#include <harp_internal.hpp>

#include <harp_test_serialize.hpp>



namespace harp {

  void test_serialize ( std::string const & datadir );

  void test_linalg ( std::string const & datadir );

  void test_spec_simspecter ( std::string const & datadir );

  void test_psf_gauss_sim ( std::string const & datadir );

  void test_psf_gauss ( std::string const & datadir );

  void test_image_simfits ( std::string const & datadir );

  void test_specslice ( std::string const & datadir );

  void test_extract ( std::string const & datadir );

}

  
#endif