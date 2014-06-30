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

#include <harp/static_plugins.hpp>


using namespace std;
using namespace harp;


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
  shared_ext = LT_MODULE_EXT;
#endif

  boost::cmatch what;
  boost::regex expr;
  expr = "harp_plugin_.*.so";

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


void harp::plugin_registry::register_image ( std::string const & type, image_factory create, std::string const & version ) {
  if ( image_plugins_.count ( type ) > 0 ) {
    ostringstream o;
    o << "image plugin \"" << type << "\" is already registered";
    HARP_THROW( o.str().c_str() );
  }

  string const & internal = source_version();
  if ( version != internal ) {
    ostringstream o;
    o << "image plugin \"" << type << "\" was compiled with a different version HARP, refusing to load";
    HARP_THROW( o.str().c_str() );
  }

  image_plugins_[ type ] = create;

  return;
}


void harp::plugin_registry::register_spec ( std::string const & type, spec_factory create, std::string const & version ) {
  if ( spec_plugins_.count ( type ) > 0 ) {
    ostringstream o;
    o << "spec plugin \"" << type << "\" is already registered";
    HARP_THROW( o.str().c_str() );
  }

  string const & internal = source_version();
  if ( version != internal ) {
    ostringstream o;
    o << "spec plugin \"" << type << "\" was compiled with a different version HARP, refusing to load";
    HARP_THROW( o.str().c_str() );
  }

  spec_plugins_[ type ] = create;

  return;
}


void harp::plugin_registry::register_psf ( std::string const & type, psf_factory create, std::string const & version ) {
  if ( psf_plugins_.count ( type ) > 0 ) {
    ostringstream o;
    o << "psf plugin \"" << type << "\" is already registered";
    HARP_THROW( o.str().c_str() );
  }

  string const & internal = source_version();
  if ( version != internal ) {
    ostringstream o;
    o << "psf plugin \"" << type << "\" was compiled with a different version HARP, refusing to load";
    HARP_THROW( o.str().c_str() );
  }

  psf_plugins_[ type ] = create;

  return;
}


void harp::plugin_registry::register_targets ( std::string const & type, targets_factory create, std::string const & version ) {
  if ( targets_plugins_.count ( type ) > 0 ) {
    ostringstream o;
    o << "targets plugin \"" << type << "\" is already registered";
    HARP_THROW( o.str().c_str() );
  }

  string const & internal = source_version();
  if ( version != internal ) {
    ostringstream o;
    o << "targets plugin \"" << type << "\" was compiled with a different version HARP, refusing to load";
    HARP_THROW( o.str().c_str() );
  }

  targets_plugins_[ type ] = create;

  return;
}


harp::plugin_registry::plugin_registry ( bool debug ) {

  string prefix = "  harp plugin_registry:  ";

  // register internal, static plugins

  #include "../plugins/plugin_register.cpp"

  // parse the plugin path and get the list of files to use

  path_ = "";
  files_.clear();

#ifdef USE_PLUGINS

  char * envval = getenv ( plugin_path );
  if ( envval ) {
    path_ = envval;
  }

  if ( debug ) {
    cerr << prefix << "using " << plugin_path << " = " << path_ << endl;
  }

  size_t offset = 0;
  size_t len;
  size_t split = 0;

  vector < string > dirs;

  if ( path_ != "" ) {

    while ( split != string::npos ) {

      split = path_.find ( ':', offset );

      if ( split != string::npos ) {
        len = split - offset;
        dirs.push_back ( path_.substr ( offset, len ) );
        offset = split + 1;
      }

    }

    dirs.push_back ( path_.substr ( offset, string::npos ) );

  }

  for ( vector < string > :: iterator it = dirs.begin(); it != dirs.end(); ++it ) {

    string dirpath = (*it);
    if ( dirpath[ dirpath.size() - 1 ] == '/' ) {
      dirpath.erase ( dirpath.size() - 1, 1 );
    }

    if ( debug ) {
      cerr << prefix << "searching " << dirpath << " for plugins" << endl;
    }

    find_dlls ( dirpath, files_ );
    
  }

  // go through every plugin file and call the init() method

  // flush error buffer
  char * err = dlerror();

  for ( vector < string > :: iterator it = files_.begin(); it != files_.end(); ++it ) {
    handles_[ (*it) ] = dlopen ( it->c_str(), RTLD_NOW | RTLD_LOCAL );

    if ( debug ) {
      cerr << prefix << "loading DLL " << (*it) << endl;
    }

    err = dlerror();
    if ( err != NULL ) {
      ostringstream o;
      o << "error opening plugin \"" << (*it) << "\": " << err;
      HARP_THROW( o.str().c_str() );
    }

    // do global library symbols have an underscore prefix?
    string initname = "initialize";
    #ifdef NEED_USCORE
    initname = "_initialize";
    #endif

    // function pointer to the initialize function in the DLL
    void (*init)( void * );

    // this casting technique is recommended from the dlsym manpage...
    *(void **) (&init) = dlsym ( handles_[ (*it) ], initname.c_str() );

    err = dlerror();
    if ( err != NULL ) {
      ostringstream o;
      o << "error loading symbols from plugin \"" << (*it) << "\": " << err;
      HARP_THROW( o.str().c_str() );
    }

    void * thisvoid = static_cast < void * > (this);

    (*init)( thisvoid );

  }

#endif

}


harp::plugin_registry::~plugin_registry ( ) {

  // clear all registered plugins

  image_plugins_.clear();
  spec_plugins_.clear();
  psf_plugins_.clear();
  targets_plugins_.clear();

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


image * harp::plugin_registry::create_image ( std::string const & type, boost::property_tree::ptree const & props ) {
  if ( image_plugins_.count ( type ) == 0 ) {
    cerr << "harp:  image plugin \"" << type << "\" is not registered- did you forget to set HARP_PLUGIN_PATH?" << endl;
    ostringstream o;
    o << "image plugin \"" << type << "\" is not registered";
    HARP_THROW( o.str().c_str() );
  }
  return (*image_plugins_[type])( props );
}


psf * harp::plugin_registry::create_psf ( std::string const & type, boost::property_tree::ptree const & props ) {
  if ( psf_plugins_.count ( type ) == 0 ) {
    cerr << "harp:  psf plugin \"" << type << "\" is not registered- did you forget to set HARP_PLUGIN_PATH?" << endl;
    ostringstream o;
    o << "psf plugin \"" << type << "\" is not registered";
    HARP_THROW( o.str().c_str() );
  }
  return (*psf_plugins_[type])( props );
}


spec * harp::plugin_registry::create_spec ( std::string const & type, boost::property_tree::ptree const & props ) {
  if ( spec_plugins_.count ( type ) == 0 ) {
    cerr << "harp:  spec plugin \"" << type << "\" is not registered- did you forget to set HARP_PLUGIN_PATH?" << endl;
    ostringstream o;
    o << "spec plugin \"" << type << "\" is not registered";
    HARP_THROW( o.str().c_str() );
  }
  return (*spec_plugins_[type])( props );
}


targets * harp::plugin_registry::create_targets ( std::string const & type, boost::property_tree::ptree const & props ) {
  if ( targets_plugins_.count ( type ) == 0 ) {
    cerr << "harp:  targets plugin \"" << type << "\" is not registered- did you forget to set HARP_PLUGIN_PATH?" << endl;
    ostringstream o;
    o << "targets plugin \"" << type << "\" is not registered";
    HARP_THROW( o.str().c_str() );
  }
  return (*targets_plugins_[type])( props );
}

