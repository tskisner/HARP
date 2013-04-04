// @COPYRIGHT@

#ifndef HARP_TEST_HPP
#define HARP_TEST_HPP

#include <harp_internal.hpp>

namespace harp {

  void test_elemental ( std::string const & datadir );

  void test_invcov ( std::string const & datadir );

  //void test_tinyKLT ( std::string const & datadir );

  void test_spec_specter ( std::string const & datadir );

  void test_psf_gauss ( std::string const & datadir );
  
}

  
#endif