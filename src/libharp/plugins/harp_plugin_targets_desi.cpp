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
  telera = -1.0;
  teledec = -1.0;
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

  meta_ = fits::key_read_all ( fp );

  fits::key_parse ( meta_, targets_desi_key_tileid, tileid );

  fits::key_parse ( meta_, targets_desi_key_telera, telera );

  fits::key_parse ( meta_, targets_desi_key_teledec, teledec );

  fits::key_parse ( meta_, targets_desi_key_expid, expid );

  fits::key_parse ( meta_, targets_desi_key_night, night );

  fits::key_parse ( meta_, targets_desi_key_vdmodel, vdmodel );

  fits::key_parse ( meta_, targets_desi_key_voptics, voptics );

  fits::key_parse ( meta_, targets_desi_key_vfibvcam, vfibvcam );

  fits::key_parse ( meta_, targets_desi_key_hexpdrot, hexpdrot );

  fits::key_parse ( meta_, targets_desi_key_dateobs, dateobs );


  // read table and populate object list

  vector < string > fnames;
  fits::bin_info ( fp, nobjects_, fnames );

  vector < string > cnames = colnames();

  vector < int > colidx = fits::bin_columns ( fp, cnames );

  // get optimal number of rows to process at once

  int status = 0;
  long optimal;
  int ret = fits_get_rowsize ( fp, &optimal, &status );
  fits::check ( status );

  //cerr << "target optimal rows = " << optimal << endl;

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

      //cerr << "ctor reading rows " << offset << " - " << offset + n - 1 << endl;

      ret = fits_read_col_str ( fp, 1 + colidx[0], offset + 1, 1, n, empty, objstr, &anynul, &status );
      fits::check ( status );

      ret = fits_read_col_str ( fp, 1 + colidx[1], offset + 1, 1, n, empty, targetcat, &anynul, &status );
      fits::check ( status );

      ret = fits_read_col_lnglng ( fp, 1 + colidx[2], offset + 1, 1, n, 0, targetid, &anynul, &status );
      fits::check ( status );

      ret = fits_read_col_lnglng ( fp, 1 + colidx[3], offset + 1, 1, n, 0, targetmask, &anynul, &status );
      fits::check ( status );

      ret = fits_read_col_flt ( fp, 1 + colidx[4], offset + 1, 1, 5*n, 0, mag, &anynul, &status );
      fits::check ( status );

      ret = fits_read_col_str ( fp, 1 + colidx[5], offset + 1, 1, n, empty, filter, &anynul, &status );
      fits::check ( status );

      ret = fits_read_col_lnglng ( fp, 1 + colidx[6], offset + 1, 1, n, 0, spectroid, &anynul, &status );
      fits::check ( status );

      ret = fits_read_col_lnglng ( fp, 1 + colidx[7], offset + 1, 1, n, 0, positioner, &anynul, &status );
      fits::check ( status );

      ret = fits_read_col_int ( fp, 1 + colidx[8], offset + 1, 1, n, 0, fiber, &anynul, &status );
      fits::check ( status );

      ret = fits_read_col_flt ( fp, 1 + colidx[9], offset + 1, 1, n, 0, lambdaref, &anynul, &status );
      fits::check ( status );

      ret = fits_read_col_dbl ( fp, 1 + colidx[10], offset + 1, 1, n, 0, ra_target, &anynul, &status );
      fits::check ( status );

      ret = fits_read_col_dbl ( fp, 1 + colidx[11], offset + 1, 1, n, 0, dec_target, &anynul, &status );
      fits::check ( status );

      ret = fits_read_col_dbl ( fp, 1 + colidx[12], offset + 1, 1, n, 0, ra_obs, &anynul, &status );
      fits::check ( status );

      ret = fits_read_col_dbl ( fp, 1 + colidx[13], offset + 1, 1, n, 0, dec_obs, &anynul, &status );
      fits::check ( status );

      ret = fits_read_col_dbl ( fp, 1 + colidx[14], offset + 1, 1, n, 0, x_target, &anynul, &status );
      fits::check ( status );

      ret = fits_read_col_dbl ( fp, 1 + colidx[15], offset + 1, 1, n, 0, y_target, &anynul, &status );
      fits::check ( status );

      ret = fits_read_col_dbl ( fp, 1 + colidx[16], offset + 1, 1, n, 0, x_fvcobs, &anynul, &status );
      fits::check ( status );

      ret = fits_read_col_dbl ( fp, 1 + colidx[17], offset + 1, 1, n, 0, y_fvcobs, &anynul, &status );
      fits::check ( status );

      ret = fits_read_col_flt ( fp, 1 + colidx[18], offset + 1, 1, n, 0, x_fvcerr, &anynul, &status );
      fits::check ( status );

      ret = fits_read_col_flt ( fp, 1 + colidx[19], offset + 1, 1, n, 0, y_fvcerr, &anynul, &status );
      fits::check ( status );

      for ( size_t i = 0; i < n; ++i ) {

        objects_.push_back ( object_desi_p ( new object_desi ( desi2harp ( objstr[i] ),
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
          y_fvcerr[i] ) ) );

        /*
        cerr << offset + i << " : " << endl;
        cerr << "  " << targetcat[i] << endl;
        cerr << "  " << targetid[i] << endl;
        cerr << "  " << targetmask[i] << endl;
        cerr << "  mag = " << endl;
        cerr << "   " << mag[5*i+0] << endl;
        cerr << "   " << mag[5*i+1] << endl;
        cerr << "   " << mag[5*i+2] << endl;
        cerr << "   " << mag[5*i+3] << endl;
        cerr << "   " << mag[5*i+4] << endl;
        cerr << "  " << filter[i] << endl;
        cerr << "  " << spectroid[i] << endl;
        cerr << "  " << positioner[i] << endl;
        cerr << "  " << fiber[i] << endl;
        */

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
    temp.push_back ( boost::dynamic_pointer_cast < object_desi, object > ( *it ) );
  }
  return temp;
}


boost::property_tree::ptree harp::targets_desi::meta () const {
  return meta_;
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

  types[12] = "1D";   //RA_OBS
  units[12] = "None";

  types[13] = "1D";   //DEC_OBS
  units[13] = "None";

  types[14] = "1D";   //X_TARGET
  units[14] = "None";

  types[15] = "1D";   //Y_TARGET
  units[15] = "None";

  types[16] = "1D";   //X_FVCOBS
  units[16] = "None";

  types[17] = "1D";   //Y_FVCOBS
  units[17] = "None";

  types[18] = "1E";   //Y_FVCERR
  units[18] = "None";

  types[19] = "1E";   //X_FVCERR
  units[19] = "None";

  fits::bin_create ( fp, string("FIBERMAP"), objects.size(), names, types, units );

  // check required keywords and write

  if ( meta.count ( targets_desi_key_tileid ) ) {
    fits::key_require ( meta, targets_desi_key_tileid, "I" );
  } else {
    HARP_THROW( "you must specify the Tile ID" );
  }

  if ( meta.count ( targets_desi_key_telera ) ) {
    fits::key_require ( meta, targets_desi_key_telera, "F" );
  } else {
    HARP_THROW( "you must specify the telescope RA" );
  }

  if ( meta.count ( targets_desi_key_teledec ) ) {
    fits::key_require ( meta, targets_desi_key_teledec, "F" );
  } else {
    HARP_THROW( "you must specify the telescope DEC" );
  }

  if ( meta.count ( targets_desi_key_expid ) ) {
    fits::key_require ( meta, targets_desi_key_expid, "I" );
  } else {
    HARP_THROW( "you must specify the exposure number" );
  }

  if ( meta.count ( targets_desi_key_night ) ) {
    fits::key_require ( meta, targets_desi_key_night, "C" );
  } else {
    HARP_THROW( "you must specify the night" );
  }

  if ( meta.count ( targets_desi_key_vdmodel ) ) {
    fits::key_require ( meta, targets_desi_key_vdmodel, "C" );
  } else {
    HARP_THROW( "you must specify the desimodel version" );
  }

  if ( meta.count ( targets_desi_key_voptics ) ) {
    fits::key_require ( meta, targets_desi_key_voptics, "C" );
  } else {
    HARP_THROW( "you must specify the optics model version" );
  }

  if ( meta.count ( targets_desi_key_vfibvcam ) ) {
    fits::key_require ( meta, targets_desi_key_vfibvcam, "C" );
  } else {
    HARP_THROW( "you must specify the fiber view code version" );
  }

  if ( meta.count ( targets_desi_key_hexpdrot ) ) {
    fits::key_require ( meta, targets_desi_key_hexpdrot, "F" );
  } else {
    HARP_THROW( "you must specify the hexapod rotation angle" );
  }

  if ( meta.count ( targets_desi_key_dateobs ) ) {
    fits::key_require ( meta, targets_desi_key_dateobs, "C" );
  } else {
    HARP_THROW( "you must specify the date of observation" );
  }

  fits::key_write_all ( fp, meta );

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

        strncpy ( objstr[i], harp2desi( objects[offset+i]->type() ).c_str(), 10 );
        strncpy ( targetcat[i], objects[offset+i]->targetcat.c_str(), 20 );

        targetid[i] = objects[offset+i]->targetid;
        targetmask[i] = objects[offset+i]->targetmask;

        mag[ 5*i ] = objects[offset+i]->mag[0];
        mag[ 5*i + 1 ] = objects[offset+i]->mag[1];
        mag[ 5*i + 2 ] = objects[offset+i]->mag[2];
        mag[ 5*i + 3 ] = objects[offset+i]->mag[3];
        mag[ 5*i + 4 ] = objects[offset+i]->mag[4];

        strncpy ( filter[i], objects[offset+i]->filter.c_str(), 50 );

        spectroid[i] = objects[offset+i]->spectroid;

        positioner[i] = objects[offset+i]->positioner;

        fiber[i] = objects[offset+i]->fiber;

        lambdaref[i] = objects[offset+i]->lambdaref;

        ra_target[i] = objects[offset+i]->ra_target;
        dec_target[i] = objects[offset+i]->dec_target;

        ra_obs[i] = objects[offset+i]->ra_obs;
        dec_obs[i] = objects[offset+i]->dec_obs;

        x_target[i] = objects[offset+i]->x_target;
        y_target[i] = objects[offset+i]->y_target;

        x_fvcobs[i] = objects[offset+i]->x_fvcobs;
        y_fvcobs[i] = objects[offset+i]->y_fvcobs;

        x_fvcerr[i] = objects[offset+i]->x_fvcerr;
        y_fvcerr[i] = objects[offset+i]->y_fvcerr;

      }

      //cerr << "writing rows " << offset << " - " << offset + n - 1 << endl;

      ret = fits_write_col_str ( fp, 1, offset + 1, 1, n, objstr, &status );
      fits::check ( status );

      ret = fits_write_col_str ( fp, 2, offset + 1, 1, n, targetcat, &status );
      fits::check ( status );

      ret = fits_write_col_lnglng ( fp, 3, offset + 1, 1, n, targetid, &status );
      fits::check ( status );

      ret = fits_write_col_lnglng ( fp, 4, offset + 1, 1, n, targetmask, &status );
      fits::check ( status );

      ret = fits_write_col_flt ( fp, 5, offset + 1, 1, 5*n, mag, &status );
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
