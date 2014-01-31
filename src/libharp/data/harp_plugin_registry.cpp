// @COPYRIGHT@

#include <harp_data_internal.hpp>

#include <boost/regex.hpp>

extern "C" {
  #include <sys/stat.h>
  #include <dirent.h>
}

#ifdef USE_PLUGINS
extern "C" {
  #include <dlfcn.h>
}
#endif

#include "../plugin/harp/static_plugins.hpp"


using namespace std;
using namespace harp;


// Include the generated file that contains the git revision
#include "git-version.cpp"

// List of directories to search for plugins.  This should be colon-separated, like
// the standard PATH variable.  In each directory, we will look for files with
// names like "harp_plugin_*.so".
static const char plugin_path[] = "HARP_PLUGIN_PATH";


// plugin registry

void harp::plugin_registry::find_dlls ( string const & dirpath, vector < string > & files ) {

  DIR * dip;
  struct dirent * dit;
    
  if ( ( dip = opendir ( dirpath.c_str() ) ) == NULL ) {
    ostringstream o;
    o << "cannot open plugin directory \"" << dirpath << "\"";
    HARP_THROW( o.str().c_str() );
  }

  string shared_ext = ".so";
#ifdef LT_MODULE_EXT
  cerr << "DBG: using module extension " << LT_MODULE_EXT << endl;
  shared_ext = LT_MODULE_EXT;
#endif

  boost::cmatch what;
  boost::regex expr ( "harp_plugin_*.so" );

  while ( ( dit = readdir ( dip ) ) != NULL ) {
    // if this entry is a file with the right name convention, add to list
    if ( regex_match ( dit->d_name, what, expr ) ) {
      string path = dirpath + "/" + dit->d_name;
      files.push_back ( path );
    }
  }

  if ( closedir ( dip ) == -1 ) {
    ostringstream o;
    o << "error closing plugin directory \"" << dirpath << "\"";
    HARP_THROW( o.str().c_str() );
  }

  return;
}


void harp::plugin_registry::register_image ( std::string const & name, image_factory create ) {
  if ( image_plugins_.count ( name ) > 0 ) {
    ostringstream o;
    o << "image plugin \"" << name << "\" is already registered";
    HARP_THROW( o.str().c_str() );
  }

  image_plugins_[ name ] = create;

  cerr << "DBG: registered image plugin \"" << name << "\"" << endl;

  return;
}


void harp::plugin_registry::register_spec ( std::string const & name, spec_factory create ) {
  if ( spec_plugins_.count ( name ) > 0 ) {
    ostringstream o;
    o << "spec plugin \"" << name << "\" is already registered";
    HARP_THROW( o.str().c_str() );
  }

  spec_plugins_[ name ] = create;

  cerr << "DBG: registered spec plugin \"" << name << "\"" << endl;

  return;
}


void harp::plugin_registry::register_psf ( std::string const & name, psf_factory create ) {
  if ( psf_plugins_.count ( name ) > 0 ) {
    ostringstream o;
    o << "psf plugin \"" << name << "\" is already registered";
    HARP_THROW( o.str().c_str() );
  }

  psf_plugins_[ name ] = create;

  cerr << "DBG: registered psf plugin \"" << name << "\"" << endl;

  return;
}


harp::plugin_registry::plugin_registry ( ) {

  // register internal, static plugins

  #include "../plugin/plugin_register.cpp"

  // parse the plugin path and get the list of files to use

  path_ = "";
  files_.clear();

#ifdef USE_PLUGINS

  char * envval = getenv ( plugin_path );
  path_ = envval;

  size_t offset = 0;
  size_t len;
  size_t split = 0;

  vector < string > dirs;

  while ( split != string::npos ) {

    split = path_.find ( ':', offset );

    if ( split != string::npos ) {
      len = split - offset;
      dirs.push_back ( path_.substr ( offset, len ) );
      offset = split + 1;
    }

  }

  dirs.push_back ( path_.substr ( offset, string::npos ) );

  for ( vector < string > :: iterator it = dirs.begin(); it != dirs.end(); ++it ) {

    string dirpath = (*it);
    if ( dirpath[ dirpath.size() - 1 ] == '/' ) {
      dirpath.erase ( dirpath.size() - 1, 1 );
    }

    find_dlls ( dirpath, files_ );
    
  }

  // go through every plugin file and call the init() method

  for ( vector < string > :: iterator it = files_.begin(); it != files_.end(); ++it ) {
    handles_[ (*it) ] = dlopen ( it->c_str(), 0 );

    char * err = dlerror();
    if ( err != NULL ) {
      ostringstream o;
      o << "error opening plugin \"" << (*it) << "\": " << err;
      HARP_THROW( o.str().c_str() );
    }

    string initname = "initialize";
    #ifdef NEED_USCORE
    initname = "_initialize";
    #endif

    void (*init)();
    *(void **) (&init) = dlsym ( handles_[ (*it) ], initname.c_str() );

    err = dlerror();
    if ( err != NULL ) {
      ostringstream o;
      o << "error loading symbols from plugin \"" << (*it) << "\": " << err;
      HARP_THROW( o.str().c_str() );
    }

    (*init)();

  }

#endif

}


harp::plugin_registry::~plugin_registry ( ) {

  // clear all registered plugins

  image_plugins_.clear();
  spec_plugins_.clear();
  psf_plugins_.clear();

  // dlclose all plugin files

  #ifdef USE_PLUGINS

  for ( std::map < std::string, void * > :: iterator it = handles_.begin(); it != handles_.end(); ++it ) {
    int ret = dlclose ( it->second );
    char * err = dlerror();

    if ( ret != 0 ) {
      ostringstream o;
      o << "error closing plugin \"" << it->first << "\": " << err;
      HARP_THROW( o.str().c_str() );
    }

  }

  #endif

}


image * harp::plugin_registry::create_image ( std::string const & name, boost::property_tree::ptree const & props ) {
  if ( image_plugins_.count ( name ) == 0 ) {
    ostringstream o;
    o << "image plugin \"" << name << "\" is not registered";
    HARP_THROW( o.str().c_str() );
  }
  return (*image_plugins_[name])( props );
}


psf * harp::plugin_registry::create_psf ( std::string const & name, boost::property_tree::ptree const & props ) {
  if ( psf_plugins_.count ( name ) == 0 ) {
    ostringstream o;
    o << "psf plugin \"" << name << "\" is not registered";
    HARP_THROW( o.str().c_str() );
  }
  return (*psf_plugins_[name])( props );
}


spec * harp::plugin_registry::create_spec ( std::string const & name, boost::property_tree::ptree const & props ) {
  if ( spec_plugins_.count ( name ) == 0 ) {
    ostringstream o;
    o << "spec plugin \"" << name << "\" is not registered";
    HARP_THROW( o.str().c_str() );
  }
  return (*spec_plugins_[name])( props );
}




