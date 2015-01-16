/*
  High Performance Astrophysical Reconstruction and Processing (HARP)
  (c) 2014-2015, The Regents of the University of California, 
  through Lawrence Berkeley National Laboratory.  See top
  level LICENSE file for details.
*/


#include <harp/static_plugins.hpp>

using namespace std;
using namespace harp;


harp::object_desi::object_desi ( ) : object () {

}

harp::object_desi::object_desi ( object_type type,
        std::string in_targetcat,
        int64_t in_targetid,
        int64_t in_targetmask,
        float * in_mag,
        std::string in_filter,
        int64_t in_spectroid,
        int64_t in_positioner,
        int32_t in_fiber,
        float in_lambdaref,
        double in_ra_target,
        double in_dec_target,
        double in_ra_obs,
        double in_dec_obs,
        double in_x_target,
        double in_y_target,
        double in_x_fvcobs,
        double in_y_fvcobs,
        float in_x_fvcerr,
        float in_y_fvcerr ) : object () {

  set_type ( type );

  targetcat = in_targetcat;
  targetid = in_targetid;
  targetmask = in_targetmask;
  mag[0] = in_mag[0];
  mag[1] = in_mag[1];
  mag[2] = in_mag[2];
  mag[3] = in_mag[3];
  mag[4] = in_mag[4];

  filter = in_filter;
  spectroid = in_spectroid;
  positioner = in_positioner;
  fiber = in_fiber;
  lambdaref = in_lambdaref;
  ra_target = in_ra_target;
  dec_target = in_dec_target;
  ra_obs = in_ra_obs;
  dec_obs = in_dec_obs;
  x_target = in_x_target;
  y_target = in_y_target;
  x_fvcobs = in_x_fvcobs;
  y_fvcobs = in_y_fvcobs;
  x_fvcerr = in_x_fvcerr;
  y_fvcerr = in_y_fvcerr;

  ostringstream o;
  o << targetid;
  set_name ( o.str() );

}


static const char * targets_desi_key_path = "path";


harp::targets_desi::targets_desi ( ) : targets () {
  nobjects_ = 0;
  objects_.clear();
  hdu_ = 1;
  path_ = "";
  tileid = -1;
  tilera = -1.0;
  tiledec = -1.0;
  expid = -1;
  night = "";
  vdmodel = "";
  voptics = "";
  vfibvcam = "";
  hexpdrot = -1.0;
  dateobs = "";
}


harp::targets_desi::targets_desi ( boost::property_tree::ptree const & props ) : targets ( "desi", props ) {

  path_ = props.get < string > ( targets_desi_key_path, "" );

  hdu_ = 2;

  // seek to table

  fitsfile * fp;

  fits::open_read ( fp, path_ );

  fits::bin_seek ( fp, hdu_ );

  // read keywords



  // read table and populate object list

  /*
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

  */
  
}


harp::targets_desi::~targets_desi ( ) {
  
}


size_t harp::targets_desi::n_objects ( ) const {
  return nobjects_;
}


vector < object_p > harp::targets_desi::objects ( ) const {
  vector < object_p > temp;
  for ( vector < object_desi_p > :: const_iterator it = objects_.begin(); it != objects_.end(); ++it ) {
    //temp.push_back ( (*it)-> )
  }
  return temp;
}


vector < object_desi_p > harp::targets_desi::desi_objects ( ) const {
  return objects_;
}


void harp::targets_desi::write ( std::string const & path, boost::property_tree::ptree const & meta, vector < object_desi_p > const & objects ) {

  fitsfile * fp;
    
  fits::create ( fp, path );

  /*
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
  */

  fits::close ( fp );

  return;
}


BOOST_CLASS_EXPORT(harp::targets_desi)

targets * harp::targets_desi_create ( boost::property_tree::ptree const & props ) {
  return new targets_desi ( props );
}
