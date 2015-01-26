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
static const char * targets_desi_key_tileid = "TILEID";
static const char * targets_desi_key_telera = "TELERA";
static const char * targets_desi_key_teledec = "TELEDEC";
static const char * targets_desi_key_expid = "EXPID";
static const char * targets_desi_key_night = "NIGHT";
static const char * targets_desi_key_vdmodel = "VDMODEL";
static const char * targets_desi_key_voptics = "VOPTICS";
static const char * targets_desi_key_vfibvcam = "VFIBVCAM";
static const char * targets_desi_key_hexpdrot = "HEXPDROT";
static const char * targets_desi_key_dateobs = "DATE-OBS";




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

  fits::key_read ( fp, string(targets_desi_key_tileid), tileid );

  fits::key_read ( fp, string(targets_desi_key_telera), tilera );

  fits::key_read ( fp, string(targets_desi_key_teledec), tiledec );

  fits::key_read ( fp, string(targets_desi_key_expid), expid );

  fits::key_read ( fp, string(targets_desi_key_night), night );

  fits::key_read ( fp, string(targets_desi_key_vdmodel), vdmodel );

  fits::key_read ( fp, string(targets_desi_key_voptics), voptics );

  fits::key_read ( fp, string(targets_desi_key_vfibvcam), vfibvcam );

  fits::key_read ( fp, string(targets_desi_key_hexpdrot), hexpdrot );

  fits::key_read ( fp, string(targets_desi_key_dateobs), dateobs );

  // read table and populate object list

  vector < string > fnames;
  fits::bin_info ( fp, nobjects_, fnames );

  vector < string > cols = colnames();

  vector < int > colidx;
  colidx = fits::bin_columns ( fp, cols );

  // get optimal number of rows to process at once

  int status = 0;
  long optimal;
  int ret = fits_get_rowsize ( fp, &optimal, &status );
  fits::check ( status );

  // read data in buffered way

  long long rows = (long long)nobjects_;
  long long offset = 0;
  long long n = optimal;
  int anynul;

  objects_.clear();

  while ( n == optimal ) {

    if ( offset + optimal > rows ) {
      n = rows - offset;
    }

    if ( n > 0 ) {

      char ** objstr;
      char ** targetcat;
      long long targetid[n];
      long long targetmask[n];
      float mag[5*n];
      char ** filter;
      long long spectroid[n];
      long long positioner[n];
      int fiber[n];
      float lambdaref[n];
      double ra_target[n];
      double dec_target[n];
      double ra_obs[n];
      double dec_obs[n];
      double x_target[n];
      double y_target[n];
      double x_fvcobs[n];
      double y_fvcobs[n];
      float x_fvcerr[n];
      float y_fvcerr[n];

      char empty[5];
      strcpy ( empty, "" );

      objstr = (char**) malloc ( n * sizeof(char*) );
      targetcat = (char**) malloc ( n * sizeof(char*) );
      filter = (char**) malloc ( n * sizeof(char*) );

      if ( ! ( objstr && targetcat && filter ) ) {
        ostringstream o;
        o << "cannot allocate column info for desi targets";
        HARP_THROW( o.str().c_str() );
      }

      for ( size_t c = 0; c < n; ++c ) {
        objstr[c] = (char*) malloc ( FLEN_VALUE * sizeof(char) );
        targetcat[c] = (char*) malloc ( FLEN_VALUE * sizeof(char) );
        filter[c] = (char*) malloc ( FLEN_VALUE * sizeof(char) );
        if ( ! ( objstr[c] && targetcat[c] && filter[c] ) ) {
          ostringstream o;
          o << "cannot allocate column info for desi targets";
          HARP_THROW( o.str().c_str() );
        }
      }

      ret = fits_read_col_str ( fp, colidx[0], offset + 1, 1, n, empty, objstr, &anynul, &status );
      fits::check ( status );

      ret = fits_read_col_str ( fp, colidx[1], offset + 1, 1, n, empty, targetcat, &anynul, &status );
      fits::check ( status );

      ret = fits_read_col_lnglng ( fp, colidx[2], offset + 1, 1, n, 0, targetid, &anynul, &status );
      fits::check ( status );

      ret = fits_read_col_lnglng ( fp, colidx[3], offset + 1, 1, n, 0, targetmask, &anynul, &status );
      fits::check ( status );

      ret = fits_read_col_flt ( fp, colidx[4], offset + 1, 1, n, 0, mag, &anynul, &status );
      fits::check ( status );

      ret = fits_read_col_str ( fp, colidx[5], offset + 1, 1, n, empty, filter, &anynul, &status );
      fits::check ( status );

      ret = fits_read_col_lnglng ( fp, colidx[6], offset + 1, 1, n, 0, spectroid, &anynul, &status );
      fits::check ( status );

      ret = fits_read_col_lnglng ( fp, colidx[7], offset + 1, 1, n, 0, positioner, &anynul, &status );
      fits::check ( status );

      ret = fits_read_col_int ( fp, colidx[8], offset + 1, 1, n, 0, fiber, &anynul, &status );
      fits::check ( status );

      ret = fits_read_col_flt ( fp, colidx[9], offset + 1, 1, n, 0, lambdaref, &anynul, &status );
      fits::check ( status );

      ret = fits_read_col_dbl ( fp, colidx[10], offset + 1, 1, n, 0, ra_target, &anynul, &status );
      fits::check ( status );

      ret = fits_read_col_dbl ( fp, colidx[11], offset + 1, 1, n, 0, dec_target, &anynul, &status );
      fits::check ( status );

      ret = fits_read_col_dbl ( fp, colidx[12], offset + 1, 1, n, 0, ra_obs, &anynul, &status );
      fits::check ( status );

      ret = fits_read_col_dbl ( fp, colidx[13], offset + 1, 1, n, 0, dec_obs, &anynul, &status );
      fits::check ( status );

      ret = fits_read_col_dbl ( fp, colidx[14], offset + 1, 1, n, 0, x_target, &anynul, &status );
      fits::check ( status );

      ret = fits_read_col_dbl ( fp, colidx[15], offset + 1, 1, n, 0, y_target, &anynul, &status );
      fits::check ( status );

      ret = fits_read_col_dbl ( fp, colidx[16], offset + 1, 1, n, 0, x_fvcobs, &anynul, &status );
      fits::check ( status );

      ret = fits_read_col_dbl ( fp, colidx[17], offset + 1, 1, n, 0, y_fvcobs, &anynul, &status );
      fits::check ( status );

      ret = fits_read_col_flt ( fp, colidx[18], offset + 1, 1, n, 0, x_fvcerr, &anynul, &status );
      fits::check ( status );

      ret = fits_read_col_flt ( fp, colidx[19], offset + 1, 1, n, 0, y_fvcerr, &anynul, &status );
      fits::check ( status );

      for ( size_t i = 0; i < n; ++i ) {

        object_desi_p obj ( new object_desi ( desi2harp ( objstr[i] ),
          targetcat[i],
          targetid[i],
          targetmask[i],
          &(mag[5*i]),
          filter[i],
          spectroid[i],
          positioner[i],
          fiber[i],
          lambdaref[i],
          ra_target[i],
          dec_target[i],
          ra_obs[i],
          dec_obs[i],
          x_target[i],
          y_target[i],
          x_fvcobs[i],
          y_fvcobs[i],
          x_fvcerr[i],
          y_fvcerr[i] ) );

        objects_.push_back ( boost::dynamic_pointer_cast < object_desi, object > ( obj->shared_from_this() ) );

      }

      for ( size_t c = 0; c < n; ++c ) {
        free ( objstr[c] );
        free ( targetcat[c] );
        free ( filter[c] );
      }
      free ( objstr );
      free ( targetcat );
      free ( filter );

    }

    offset += n;
  }

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

  // append empty image

  fits::img_append < double > ( fp, 0, 0 );

  // append table

  vector < string > names = colnames();
  vector < string > types ( names.size() );
  vector < string > units ( names.size() );

  types[0] = "10A";   //OBJTYPE
  units[0] = "None";

  types[1] = "20A";   //TARGETCAT
  units[1] = "None";
  
  types[2] = "1K";    //TARGETID
  units[2] = "None";
  
  types[3] = "1K";    //TARGET_MASK0
  units[3] = "None";

  types[4] = "5E";    //MAG
  units[4] = "None";

  types[5] = "50A";   //FILTER
  units[5] = "None";

  types[6] = "1K";    //SPECTROID
  units[6] = "None";

  types[7] = "1K";    //POSITIONER
  units[7] = "None";

  types[8] = "1J";   //FIBER
  units[8] = "None";

  types[9] = "1E";   //LAMBDAREF
  units[9] = "None";

  types[10] = "1D";   //RA_TARGET
  units[10] = "None";

  types[11] = "1D";   //DEC_TARGET
  units[11] = "None";

  types[0] = "1D";   //RA_OBS
  units[0] = "None";

  types[0] = "1D";   //DEC_OBS
  units[0] = "None";

  types[0] = "1D";   //X_TARGET
  units[0] = "None";

  types[0] = "1D";   //Y_TARGET
  units[0] = "None";

  types[0] = "1D";   //X_FVCOBS
  units[0] = "None";

  types[0] = "1D";   //Y_FVCOBS
  units[0] = "None";

  types[0] = "1E";   //Y_FVCERR
  units[0] = "None";

  types[0] = "1E";   //X_FVCERR
  units[0] = "None";

  fits::bin_create ( fp, string("FIBERMAP"), objects.size(), names, types, units );

  // write keywords

  if ( meta.count ( targets_desi_key_tileid ) ) {
    fits::key_write ( fp, targets_desi_key_tileid, meta.get < int > ( targets_desi_key_tileid ), "Tile ID" );
  } else {
    HARP_THROW( "you must specify the Tile ID" );
  }

  if ( meta.count ( targets_desi_key_telera ) ) {
    fits::key_write ( fp, targets_desi_key_telera, meta.get < float > ( targets_desi_key_telera ), "Telescope central RA [deg]" );
  } else {
    HARP_THROW( "you must specify the telescope RA" );
  }

  if ( meta.count ( targets_desi_key_teledec ) ) {
    fits::key_write ( fp, targets_desi_key_teledec, meta.get < float > ( targets_desi_key_teledec ), "Telescope central dec [deg]" );
  } else {
    HARP_THROW( "you must specify the telescope DEC" );
  }

  if ( meta.count ( targets_desi_key_expid ) ) {
    fits::key_write ( fp, targets_desi_key_expid, meta.get < int > ( targets_desi_key_expid ), "Exposure number" );
  } else {
    HARP_THROW( "you must specify the exposure number" );
  }

  if ( meta.count ( targets_desi_key_night ) ) {
    fits::key_write ( fp, targets_desi_key_night, meta.get < string > ( targets_desi_key_night ), "Night YEARMMDD" );
  } else {
    HARP_THROW( "you must specify the night" );
  }

  if ( meta.count ( targets_desi_key_vdmodel ) ) {
    fits::key_write ( fp, targets_desi_key_vdmodel, meta.get < string > ( targets_desi_key_vdmodel ), "desimodel version" );
  } else {
    HARP_THROW( "you must specify the desimodel version" );
  }

  if ( meta.count ( targets_desi_key_voptics ) ) {
    fits::key_write ( fp, targets_desi_key_voptics, meta.get < string > ( targets_desi_key_voptics ), "optics model version" );
  } else {
    HARP_THROW( "you must specify the optics model version" );
  }

  if ( meta.count ( targets_desi_key_vfibvcam ) ) {
    fits::key_write ( fp, targets_desi_key_vfibvcam, meta.get < string > ( targets_desi_key_vfibvcam ), "fiber view code version" );
  } else {
    HARP_THROW( "you must specify the fiber view code version" );
  }

  if ( meta.count ( targets_desi_key_hexpdrot ) ) {
    fits::key_write ( fp, targets_desi_key_hexpdrot, meta.get < float > ( targets_desi_key_hexpdrot ), "hexapod rotation [deg]" );
  } else {
    HARP_THROW( "you must specify the hexapod rotation angle" );
  }

  if ( meta.count ( targets_desi_key_dateobs ) ) {
    fits::key_write ( fp, targets_desi_key_dateobs, meta.get < string > ( targets_desi_key_dateobs ), "Date of observation in UTC" );
  } else {
    HARP_THROW( "you must specify the date of observation" );
  }

  // write data

  // get optimal number of rows to process at once

  int status = 0;
  long optimal;
  int ret = fits_get_rowsize ( fp, &optimal, &status );
  fits::check ( status );

  // write data in buffered way

  long long rows = (long long)objects.size();
  long long offset = 0;
  long long n = optimal;

  while ( n == optimal ) {

    if ( offset + optimal > rows ) {
      n = rows - offset;
    }

    if ( n > 0 ) {

      char ** objstr;
      char ** targetcat;
      long long targetid[n];
      long long targetmask[n];
      float mag[5*n];
      char ** filter;
      long long spectroid[n];
      long long positioner[n];
      int fiber[n];
      float lambdaref[n];
      double ra_target[n];
      double dec_target[n];
      double ra_obs[n];
      double dec_obs[n];
      double x_target[n];
      double y_target[n];
      double x_fvcobs[n];
      double y_fvcobs[n];
      float x_fvcerr[n];
      float y_fvcerr[n];

      objstr = (char**) malloc ( n * sizeof(char*) );
      targetcat = (char**) malloc ( n * sizeof(char*) );
      filter = (char**) malloc ( n * sizeof(char*) );

      if ( ! ( objstr && targetcat && filter ) ) {
        ostringstream o;
        o << "cannot allocate column info for desi targets";
        HARP_THROW( o.str().c_str() );
      }

      for ( size_t c = 0; c < n; ++c ) {
        objstr[c] = (char*) malloc ( FLEN_VALUE * sizeof(char) );
        targetcat[c] = (char*) malloc ( FLEN_VALUE * sizeof(char) );
        filter[c] = (char*) malloc ( FLEN_VALUE * sizeof(char) );
        if ( ! ( objstr[c] && targetcat[c] && filter[c] ) ) {
          ostringstream o;
          o << "cannot allocate column info for desi targets";
          HARP_THROW( o.str().c_str() );
        }
      }

      for ( size_t i = 0; i < n; ++i ) {

        strncpy ( objstr[i], harp2desi( objects[i]->type() ).c_str(), 10 );
        strncpy ( targetcat[i], objects[i]->targetcat.c_str(), 20 );

        targetid[i] = objects[i]->targetid;
        targetmask[i] = objects[i]->targetmask;

        mag[ 5*i ] = objects[i]->mag[0];
        mag[ 5*i + 1 ] = objects[i]->mag[1];
        mag[ 5*i + 2 ] = objects[i]->mag[2];
        mag[ 5*i + 3 ] = objects[i]->mag[3];
        mag[ 5*i + 4 ] = objects[i]->mag[4];

        strncpy ( filter[i], objects[i]->filter.c_str(), 50 );

        spectroid[i] = objects[i]->spectroid;

        positioner[i] = objects[i]->positioner;

        fiber[i] = objects[i]->fiber;

        lambdaref[i] = objects[i]->lambdaref;

        ra_target[i] = objects[i]->ra_target;
        dec_target[i] = objects[i]->dec_target;

        ra_obs[i] = objects[i]->ra_obs;
        dec_obs[i] = objects[i]->dec_obs;

        x_target[i] = objects[i]->x_target;
        y_target[i] = objects[i]->y_target;

        x_fvcobs[i] = objects[i]->x_fvcobs;
        y_fvcobs[i] = objects[i]->y_fvcobs;

        x_fvcerr[i] = objects[i]->x_fvcerr;
        y_fvcerr[i] = objects[i]->y_fvcerr;

      }

      ret = fits_write_col_str ( fp, 1, offset + 1, 1, n, objstr, &status );
      fits::check ( status );

      ret = fits_write_col_str ( fp, 2, offset + 1, 1, n, targetcat, &status );
      fits::check ( status );

      ret = fits_write_col_lnglng ( fp, 3, offset + 1, 1, n, targetid, &status );
      fits::check ( status );

      ret = fits_write_col_lnglng ( fp, 4, offset + 1, 1, n, targetmask, &status );
      fits::check ( status );

      ret = fits_write_col_flt ( fp, 5, offset + 1, 1, n, mag, &status );
      fits::check ( status );

      ret = fits_write_col_str ( fp, 6, offset + 1, 1, n, filter, &status );
      fits::check ( status );

      ret = fits_write_col_lnglng ( fp, 7, offset + 1, 1, n, spectroid, &status );
      fits::check ( status );

      ret = fits_write_col_lnglng ( fp, 8, offset + 1, 1, n, positioner, &status );
      fits::check ( status );

      ret = fits_write_col_int ( fp, 9, offset + 1, 1, n, fiber, &status );
      fits::check ( status );

      ret = fits_write_col_flt ( fp, 10, offset + 1, 1, n, lambdaref, &status );
      fits::check ( status );

      ret = fits_write_col_dbl ( fp, 11, offset + 1, 1, n, ra_target, &status );
      fits::check ( status );

      ret = fits_write_col_dbl ( fp, 12, offset + 1, 1, n, dec_target, &status );
      fits::check ( status );

      ret = fits_write_col_dbl ( fp, 13, offset + 1, 1, n, ra_obs, &status );
      fits::check ( status );

      ret = fits_write_col_dbl ( fp, 14, offset + 1, 1, n, dec_obs, &status );
      fits::check ( status );

      ret = fits_write_col_dbl ( fp, 15, offset + 1, 1, n, x_target, &status );
      fits::check ( status );

      ret = fits_write_col_dbl ( fp, 16, offset + 1, 1, n, y_target, &status );
      fits::check ( status );

      ret = fits_write_col_dbl ( fp, 17, offset + 1, 1, n, x_fvcobs, &status );
      fits::check ( status );

      ret = fits_write_col_dbl ( fp, 18, offset + 1, 1, n, y_fvcobs, &status );
      fits::check ( status );

      ret = fits_write_col_flt ( fp, 19, offset + 1, 1, n, x_fvcerr, &status );
      fits::check ( status );

      ret = fits_write_col_flt ( fp, 20, offset + 1, 1, n, y_fvcerr, &status );
      fits::check ( status );

      for ( size_t c = 0; c < n; ++c ) {
        free ( objstr[c] );
        free ( targetcat[c] );
        free ( filter[c] );
      }
      free ( objstr );
      free ( targetcat );
      free ( filter );

    }

    offset += n;
  }

  fits::close ( fp );

  return;
}


vector < string > harp::targets_desi::colnames () {
  vector < string > ret;

  ret.push_back ( "OBJTYPE" );
  ret.push_back (  "TARGETCAT" );
  ret.push_back ( "TARGETID" );
  ret.push_back ( "TARGET_MASK0" );
  ret.push_back ( "MAG" );
  ret.push_back ( "FILTER" );
  ret.push_back ( "SPECTROID" );
  ret.push_back ( "POSITIONER" );
  ret.push_back ( "FIBER" );
  ret.push_back ( "LAMBDAREF" );
  ret.push_back ( "RA_TARGET" );
  ret.push_back ( "DEC_TARGET" );
  ret.push_back ( "RA_OBS" );
  ret.push_back ( "DEC_OBS" );
  ret.push_back ( "X_TARGET" );
  ret.push_back ( "Y_TARGET" );
  ret.push_back ( "X_FVCOBS" );
  ret.push_back ( "Y_FVCOBS" );
  ret.push_back ( "X_FVCERR" );
  ret.push_back ( "Y_FVCERR" );

  return ret;
}


object_type harp::targets_desi::desi2harp ( std::string const & name ) {
  object_type ret;
  if ( name.compare ( "STD" ) == 0 ) {
    ret = OBJECT_CALIB;
  } else if ( name.compare ( "SKY" ) == 0 ) {
    ret = OBJECT_SKY;
  } else if ( name.compare ( "LRG" ) == 0 ) {
    ret = OBJECT_LRG;
  } else if ( name.compare ( "ELG" ) == 0 ) {
    ret = OBJECT_ELG;
  } else if ( name.compare ( "QSO" ) == 0 ) {
    ret = OBJECT_QSO;
  } else if ( name.compare ( "STAR" ) == 0 ) {
    ret = OBJECT_STAR;
  } else {
    ret = OBJECT_UNKNOWN;
  }
  return ret;
}


std::string harp::targets_desi::harp2desi ( object_type type ) {
  string ret;
  switch ( type ) {
    case OBJECT_CALIB :
      ret = "STD";
      break;
    case OBJECT_SKY :
      ret = "SKY";
      break;
    case OBJECT_LRG :
      ret = "LRG";
      break;
    case OBJECT_ELG :
      ret = "ELG";
      break;
    case OBJECT_QSO :
      ret = "QSO";
      break;
    case OBJECT_STAR :
      ret = "STAR";
      break;
    default :
      ret = "UNKNOWN";
  }
  return ret;
}



BOOST_CLASS_EXPORT(harp::targets_desi)

targets * harp::targets_desi_create ( boost::property_tree::ptree const & props ) {
  return new targets_desi ( props );
}
