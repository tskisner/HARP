// @COPYRIGHT@

#ifndef HARP_METADATA_HPP
#define HARP_METADATA_HPP


#include <boost/property_tree/json_parser.hpp>

namespace harp { namespace meta {

  class psfgroup {

    public :
      psfgroup ( );
      psfgroup ( boost::ptree & tree );
      ~psfgroup ( ) { }

      void load ( boost::ptree & tree );
      boost::ptree dump ( );

      psf_p psf ( ) { return handle_; }
      std::list < image_p > images ( ) { return imgs_; }
      
    private :
      psf_p handle_;
      std::list < image_p > imgs_;

  };

  typedef boost::shared_ptr < psfgroup > psfgroup_p;

  std::list < psfgroup_p > job_load_extraction ( std::string const & path );

  std::list < psfgroup_p > job_load_extraction ( std::istream & strm );

  void job_dump_extraction ( std::list < psfgroup_p > & psfs, std::string const & path );

  void job_dump_extraction ( std::list < psfgroup_p > & psfs, std::ostream & strm );  




} }

#endif

