// @COPYRIGHT@

#ifndef HARP_TEST_HPP
#define HARP_TEST_HPP

#include <harp_internal.hpp>

#include <harp_test_serialize.hpp>

// boost serialization headers for dump / load tests
#include <fstream>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

namespace harp {

  void test_serialize ( std::string const & datadir );

  void test_spec_specter ( std::string const & datadir );

  void test_spec_sim ( std::string const & datadir );

  //void test_image_fits ( std::string const & datadir );

  //void test_image_sim ( std::string const & datadir );

  //void test_psf_gauss ( std::string const & datadir );

  //void test_psf_gauss_sim ( std::string const & datadir );

}

  
#endif