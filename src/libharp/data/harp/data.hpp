// @COPYRIGHT@

#ifndef HARP_DATA_HPP
#define HARP_DATA_HPP

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree_serialization.hpp>

#include <boost/smart_ptr/shared_ptr.hpp>
#include <boost/smart_ptr/weak_ptr.hpp>
#include <boost/enable_shared_from_this.hpp>

#include <boost/serialization/string.hpp>
#include <boost/serialization/version.hpp>
#include <boost/serialization/split_member.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/nvp.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/export.hpp>
#include <boost/serialization/assume_abstract.hpp>

#include <boost/serialization/list.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/map.hpp>

#include <boost/foreach.hpp>

#include <harp/math.hpp>

namespace harp {

  // utilities

  std::string ptree_quote ( const std::string & s );
  void ptree_print ( const boost::property_tree::ptree & pt, int level );
  void ptree_print ( const boost::property_tree::ptree & pt );

}

#include <harp/fits.hpp>

#include <harp/target.hpp>
#include <harp/image.hpp>
#include <harp/spec.hpp>
#include <harp/psf.hpp>

namespace harp {

  // These are convenience typedefs for the function signatures needed during plugin registration

  typedef harp::image * (*image_factory)( boost::property_tree::ptree const & );
  typedef harp::spec * (*spec_factory)( boost::property_tree::ptree const & );
  typedef harp::psf * (*psf_factory)( boost::property_tree::ptree const & );

  // plugin registry

  class plugin_registry {
    
    public :

      ~plugin_registry ( );

      static plugin_registry & get ( ) {
        static plugin_registry instance;
        return instance;
      }

      void register_image ( std::string const & name, image_factory create );
      void register_spec ( std::string const & name, spec_factory create );
      void register_psf ( std::string const & name, psf_factory create );

      image * create_image ( std::string const & name, boost::property_tree::ptree const & props );
      psf * create_psf ( std::string const & name, boost::property_tree::ptree const & props );
      spec * create_spec ( std::string const & name, boost::property_tree::ptree const & props );
      
    private :

      plugin_registry ( );
      
      // do not implement these two, to prevent copying
      plugin_registry ( plugin_registry const & orig );
      void operator= ( plugin_registry const & rhs );

      void find_dlls ( std::string const & dirpath, std::vector < std::string > & files );

      std::string path_;
      std::vector < std::string > files_;
      std::map < std::string, image_factory > image_plugins_;
      std::map < std::string, spec_factory > spec_plugins_;
      std::map < std::string, psf_factory > psf_plugins_;
      std::map < std::string, void * > handles_;
      
  };

}


#endif
