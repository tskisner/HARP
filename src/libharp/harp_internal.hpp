// @COPYRIGHT@

#ifndef HARP_INTERNAL_HPP
#define HARP_INTERNAL_HPP

#include <config.h>
#include <harp_math_internal.hpp>
#include <harp.hpp>


namespace harp {

  class plugin_registry {
    
    public :

      ~plugin_registry ( );

      static plugin_registry & get ( ) {
        static plugin_registry instance;
        return instance;
      }

      
    private :

      plugin_registry ( );
      
      // do not implement these two, to prevent copying
      plugin_registry ( plugin_registry const & orig );
      void operator= ( plugin_registry const & rhs );

      void find_dlls ( std::string const & dirpath, std::vector < std::string > & files );

      std::string path_;
      std::vector < std::string > files_;

      
  };
  
  
}

  
#endif