/*
  High Performance Astrophysical Reconstruction and Processing (HARP)
  (c) 2014-2015, The Regents of the University of California, 
  through Lawrence Berkeley National Laboratory.  See top
  level LICENSE file for details.
*/

#ifndef HARP_METADATA_HPP
#define HARP_METADATA_HPP

#include <list>


namespace harp {

  spec_p load_spec ( boost::property_tree::ptree const & tree );

  image_p load_image ( boost::property_tree::ptree const & tree );

  psf_p load_psf ( boost::property_tree::ptree const & tree );

  targets_p load_targets ( boost::property_tree::ptree const & tree );


  class group {

    public :
      
      group ( boost::property_tree::ptree const & tree );

      void load ( boost::property_tree::ptree const & tree );

      psf_p psf ( ) { return handle_; }
      std::list < image_p > images ( ) { return imgs_; }
      
    private :
      psf_p handle_;
      std::list < image_p > imgs_;

  };

  typedef boost::shared_ptr < group > group_p;

  std::list < group_p > load_groups ( std::string const & path );

}

#endif

