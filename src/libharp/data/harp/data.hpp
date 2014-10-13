// @COPYRIGHT@

#ifndef HARP_DATA_HPP
#define HARP_DATA_HPP

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree_serialization.hpp>

#include <boost/foreach.hpp>

#include <harp/math.hpp>

namespace harp {

  // utilities

  std::string ptree_quote ( const std::string & s );
  void ptree_print ( const boost::property_tree::ptree & pt, int level );
  void ptree_print ( const boost::property_tree::ptree & pt );

  std::string const & source_version ( );

}

#include <harp/fits.hpp>

#include <harp/targets.hpp>
#include <harp/image.hpp>
#include <harp/spec.hpp>
#include <harp/psf.hpp>

namespace harp {

  // These are convenience typedefs for the function signatures needed during plugin registration

  typedef harp::image * (*image_factory)( boost::property_tree::ptree const & );
  typedef harp::spec * (*spec_factory)( boost::property_tree::ptree const & );
  typedef harp::psf * (*psf_factory)( boost::property_tree::ptree const & );
  typedef harp::targets * (*targets_factory)( boost::property_tree::ptree const & );

  // plugin registry

  class plugin_registry {
    
    public :

      ~plugin_registry ( );

      static plugin_registry & get ( bool debug = false ) {
        static plugin_registry instance ( debug );
        return instance;
      }

      void register_image ( std::string const & type, image_factory create, std::string const & version );
      void register_spec ( std::string const & type, spec_factory create, std::string const & version );
      void register_psf ( std::string const & type, psf_factory create, std::string const & version );
      void register_targets ( std::string const & type, targets_factory create, std::string const & version );

      image * create_image ( std::string const & type, boost::property_tree::ptree const & props );
      psf * create_psf ( std::string const & type, boost::property_tree::ptree const & props );
      spec * create_spec ( std::string const & type, boost::property_tree::ptree const & props );
      targets * create_targets ( std::string const & type, boost::property_tree::ptree const & props );
      
    private :

      plugin_registry ( bool debug );
      
      // do not implement these two, to prevent copying
      plugin_registry ( plugin_registry const & orig );
      void operator= ( plugin_registry const & rhs );

      void find_dlls ( std::string const & dirpath, std::vector < std::string > & files );

      std::string path_;
      std::vector < std::string > files_;
      std::map < std::string, image_factory > image_plugins_;
      std::map < std::string, spec_factory > spec_plugins_;
      std::map < std::string, psf_factory > psf_plugins_;
      std::map < std::string, targets_factory > targets_plugins_;
      std::map < std::string, void * > handles_;
      
  };

}

#include <harp/metadata.hpp>


#endif
