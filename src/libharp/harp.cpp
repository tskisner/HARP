// @COPYRIGHT@

#include <harp_internal.hpp>

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



using namespace std;
using namespace harp;


// Include the generated file that contains the git revision
#include "git-version.cpp"

// List of directories to search for plugins.  This should be colon-separated, like
// the standard PATH variable.  In each directory, we will look for files with
// names like "harp_plugin_*.so".
static const char plugin_path[] = "HARP_PLUGIN_PATH";


harp::exception::exception ( char const * msg, char const * file, int line ) {
  int ret;
  ret = snprintf ( msg_, BIGSTRLEN, "Exception at line %d of file %s:  %s", line, file, msg );
}


harp::exception::~exception ( ) throw() { }


const char* harp::exception::what() const throw() { 
  return msg_;
}


// utilities

std::string harp::ptree_quote ( const std::string & s ) {
  return "\"" + s + "\"";
}

void harp::ptree_print ( const boost::property_tree::ptree & pt, int level ) {
  const std::string sep ( 2 * level, ' ' );
  BOOST_FOREACH ( const boost::property_tree::ptree::value_type & v, pt ) {
    std::cerr << sep << ptree_quote ( v.first ) << " : " << ptree_quote ( v.second.data() ) << "\n";
    ptree_print ( v.second, level + 1 );
  }
  return;
}

void harp::ptree_print ( const boost::property_tree::ptree & pt ) {
  ptree_print ( pt, 0 );
  return;
}


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


harp::plugin_registry::plugin_registry ( ) {

  // register internal, static plugins





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
      dirpath.pop_back();
    }

    find_dlls ( dirpath, files_ );
    
  }



#endif

}


harp::plugin_registry::~plugin_registry ( ) {

  // clear all registered plugins



  // dlclose all plugin files


}




