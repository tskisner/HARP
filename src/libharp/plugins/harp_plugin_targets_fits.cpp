/*
  High Performance Astrophysical Reconstruction and Processing (HARP)
  (c) 2014-2015, The Regents of the University of California, 
  through Lawrence Berkeley National Laboratory.  See top
  level LICENSE file for details.
*/


#include <harp/static_plugins.hpp>

using namespace std;
using namespace harp;


static const char * targets_fits_key_path = "path";
static const char * targets_fits_key_hdu = "hdu";


harp::targets_fits::targets_fits ( ) : targets () {
  nobjects_ = 0;
  objects_.clear();
  hdu_ = 1;
  path_ = "";
}


harp::targets_fits::targets_fits ( boost::property_tree::ptree const & props ) : targets ( "fits", props ) {

  path_ = props.get < string > ( targets_fits_key_path, "" );

  hdu_ = props.get < int > ( targets_fits_key_hdu, 2 );

  if ( path_ != "" ) {

    fitsfile * fp;

    fits::open_read ( fp, path_ );

    fits::bin_seek ( fp, hdu_ );

    vector < string > colnames;
    fits::bin_info ( fp, nobjects_, colnames );

    vector < string > objcolnames(1);
    objcolnames[0] = "OBJTYPE";

    vector < int > objcols;
    objcols = fits::bin_columns ( fp, objcolnames );

    vector < string > objtypes;
    fits::bin_read_column_strings ( fp, 0, nobjects_ - 1, objcols[0], objtypes );

    objects_.clear();

    for ( size_t i = 0; i < nobjects_; ++i ) {
      object_type type = object_str2type ( objtypes[i] );
      objects_.push_back ( object_p ( new object ( type, "" ) ) );
    }

  }
  
}


harp::targets_fits::~targets_fits ( ) {
  
}


size_t harp::targets_fits::n_objects ( ) const {
  return nobjects_;
}


std::vector < object_p > harp::targets_fits::objects ( ) const {
  return objects_;
}


void harp::targets_fits::write ( std::string const & path, int hdu, std::vector < object_p > const & objects ) {

  fitsfile * fp;
    
  fits::create ( fp, path );

  int cur = 1;

  while ( cur < hdu ) {
    fits::img_append < double > ( fp, 0, 0 );
    ++cur;
  }

  vector < string > colnames ( 1 );
  vector < string > coltypes ( 1 );
  vector < string > colunits ( 1 );

  colnames[0] = "OBJTYPE";
  coltypes[0] = "8A";
  colunits[0] = "None";

  fits::bin_create ( fp, string("TARGETINFO"), objects.size(), colnames, coltypes, colunits );

  vector < string > objnames ( objects.size() );

  for ( size_t i = 0; i < objects.size(); ++i ) {
    objnames[i] = object_type2str ( objects[i]->type() );
  }

  fits::bin_write_column_strings ( fp, 0, objects.size() - 1, 0, objnames );

  fits::close ( fp );

  return;
}


BOOST_CLASS_EXPORT(harp::targets_fits)

targets * harp::targets_fits_create ( boost::property_tree::ptree const & props ) {
  return new targets_fits ( props );
}
