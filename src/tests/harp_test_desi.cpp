#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include <harp_test.hpp>


extern "C" {
  #include <sys/stat.h>
  #include <unistd.h>
}

using namespace std;
using namespace harp;


void harp::test_desi ( string const & datadir ) {

  plugin_registry & reg = plugin_registry::get();

  string desidir = datadir + "/desi";

  struct stat statbuf;
  int ret;
  
  ret = stat ( desidir.c_str(), &statbuf );
  
  if ( ret != 0 ) {
    cerr << "Skipping DESI tests (test data not found)" << endl;
    return;
  }

  cout << "Testing DESI file formats..." << endl;

  string desipix = desidir + "/pix-r0-00000000.fits";
  string desitarg = desidir + "/fibermap-00000002.fits";

  boost::property_tree::ptree props;

  // read, write, read image
  
  props.clear();
  props.put ( "path", desipix );

  boost::shared_ptr < image_desi > img ( new image_desi( props ) );

  string chk_desipix = datadir + "/desi_check_pix.fits.out";

  vector_double data;
  vector_double invvar;
  vector_mask msk;

  img->values ( data );
  img->inv_variance ( invvar );
  img->mask ( msk );

  image_desi::write ( chk_desipix, img->meta(), img->n_rows(), data, invvar, msk );

  vector_double check_data;
  vector_double check_invvar;
  vector_mask check_msk;

  props.clear();
  props.put ( "path", chk_desipix );

  boost::shared_ptr < image_desi > check_img ( new image_desi( props ) );

  check_img->values ( check_data );
  check_img->inv_variance ( check_invvar );
  check_img->mask ( check_msk );

  for ( size_t i = 0; i < check_msk.size(); ++i ) {
    if ( check_msk[i] != msk[i] ) {
      cerr << "FAIL on mask[" << i << "], " << check_msk[i] << " != " << msk[i];
      exit(1);
    }
  }

  for ( size_t i = 0; i < check_data.size(); ++i ) {
    if ( fabs ( ( check_data[i] - data[i] ) / data[i] ) > 1.0e-8 ) {
      cerr << "FAIL on data[" << i << "], " << check_data[i] << " != " << data[i];
      exit(1);
    }
  }

  for ( size_t i = 0; i < check_invvar.size(); ++i ) {
    if ( fabs ( ( check_invvar[i] - invvar[i] ) / invvar[i] ) > 1.0e-8 ) {
      cerr << "FAIL on invvar[" << i << "], " << check_invvar[i] << " != " << invvar[i];
      exit(1);
    }
  }


  props.clear();
  props.put ( "path", desitarg );

  boost::shared_ptr < targets_desi > targ ( new targets_desi( props ) );

  std::vector < object_desi_p > objs = targ->desi_objects();

  string chk_desitarg = datadir + "/desi_check_target.fits.out";

  targets_desi::write ( chk_desitarg, targ->meta(), objs );

  props.clear();
  props.put ( "path", chk_desitarg );

  boost::shared_ptr < targets_desi > check_targ ( new targets_desi( props ) );

  std::vector < object_desi_p > check_objs = check_targ->desi_objects();  

  for ( size_t i = 0; i < objs.size(); ++i ) {

    /*
    cerr << i << " : " << endl;
    cerr << "  " << objs[i]->targetcat << " : " << check_objs[i]->targetcat << endl;
    cerr << "  " << objs[i]->targetid << " : " << check_objs[i]->targetid << endl;
    cerr << "  " << objs[i]->targetmask << " : " << check_objs[i]->targetmask << endl;
    cerr << "  mag = " << endl;
    cerr << "   " << objs[i]->mag[0] << " : " << check_objs[i]->mag[0] << endl;
    cerr << "   " << objs[i]->mag[1] << " : " << check_objs[i]->mag[1] << endl;
    cerr << "   " << objs[i]->mag[2] << " : " << check_objs[i]->mag[2] << endl;
    cerr << "   " << objs[i]->mag[3] << " : " << check_objs[i]->mag[3] << endl;
    cerr << "   " << objs[i]->mag[4] << " : " << check_objs[i]->mag[4] << endl;
    cerr << "  " << objs[i]->filter << " : " << check_objs[i]->filter << endl;
    cerr << "  " << objs[i]->spectroid << " : " << check_objs[i]->spectroid << endl;
    cerr << "  " << objs[i]->positioner << " : " << check_objs[i]->positioner << endl;
    cerr << "  " << objs[i]->fiber << " : " << check_objs[i]->fiber << endl;
    */

    /*
    if ( check_objs[i]->targetcat != objs[i]->targetcat ) {
      cerr << "obj[" << i << "] targetcat \"" << check_objs[i]->targetcat << "\" != \"" << objs[i]->targetcat << "\"" << endl;
      exit(1);
    }
    */

    if ( check_objs[i]->targetid != objs[i]->targetid ) {
      cerr << "obj[" << i << "] targetid \"" << check_objs[i]->targetid << "\" != \"" << objs[i]->targetid << "\"" << endl;
      exit(1);
    }

    if ( check_objs[i]->targetmask != objs[i]->targetmask ) {
      cerr << "obj[" << i << "] targetmask \"" << check_objs[i]->targetmask << "\" != \"" << objs[i]->targetmask << "\"" << endl;
      exit(1);
    }

    for ( size_t j = 0; j < 5; ++j ) {
      if ( fabs ( ( check_objs[i]->mag[j] - objs[i]->mag[j] ) / objs[i]->mag[j] ) > 1.0e-5 )  {
        cerr << "obj[" << i << "] mag[" << j << "] \"" << check_objs[i]->mag[j] << "\" != \"" << objs[i]->mag[j] << "\"" << endl;
        exit(1);
      }
    }

    /*
    if ( check_objs[i]->filter != objs[i]->filter ) {
      cerr << "obj[" << i << "] filter \"" << check_objs[i]->filter << "\" != \"" << objs[i]->filter << "\"" << endl;
      exit(1);
    }
    */

    if ( check_objs[i]->spectroid != objs[i]->spectroid ) {
      cerr << "obj[" << i << "] spectroid \"" << check_objs[i]->spectroid << "\" != \"" << objs[i]->spectroid << "\"" << endl;
      exit(1);
    }

    if ( check_objs[i]->positioner != objs[i]->positioner ) {
      cerr << "obj[" << i << "] positioner \"" << check_objs[i]->positioner << "\" != \"" << objs[i]->positioner << "\"" << endl;
      exit(1);
    }

    if ( check_objs[i]->fiber != objs[i]->fiber ) {
      cerr << "obj[" << i << "] fiber \"" << check_objs[i]->fiber << "\" != \"" << objs[i]->fiber << "\"" << endl;
      exit(1);
    }

    if ( fabs ( ( check_objs[i]->lambdaref - objs[i]->lambdaref ) / objs[i]->lambdaref ) > 1.0e-5 )  {
      cerr << "obj[" << i << "] lambdaref \"" << check_objs[i]->lambdaref << "\" != \"" << objs[i]->lambdaref << "\"" << endl;
      exit(1);
    }

    if ( fabs ( ( check_objs[i]->ra_target - objs[i]->ra_target ) / objs[i]->ra_target ) > 1.0e-5 )  {
      cerr << "obj[" << i << "] ra_target \"" << check_objs[i]->ra_target << "\" != \"" << objs[i]->ra_target << "\"" << endl;
      exit(1);
    }

    if ( fabs ( ( check_objs[i]->dec_target - objs[i]->dec_target ) / objs[i]->dec_target ) > 1.0e-5 )  {
      cerr << "obj[" << i << "] dec_target \"" << check_objs[i]->dec_target << "\" != \"" << objs[i]->dec_target << "\"" << endl;
      exit(1);
    }

    if ( fabs ( ( check_objs[i]->ra_obs - objs[i]->ra_obs ) / objs[i]->ra_obs ) > 1.0e-5 )  {
      cerr << "obj[" << i << "] ra_obs \"" << check_objs[i]->ra_obs << "\" != \"" << objs[i]->ra_obs << "\"" << endl;
      exit(1);
    }

    if ( fabs ( ( check_objs[i]->dec_obs - objs[i]->dec_obs ) / objs[i]->dec_obs ) > 1.0e-5 )  {
      cerr << "obj[" << i << "] dec_obs \"" << check_objs[i]->dec_obs << "\" != \"" << objs[i]->dec_obs << "\"" << endl;
      exit(1);
    }

    if ( fabs ( ( check_objs[i]->x_target - objs[i]->x_target ) / objs[i]->x_target ) > 1.0e-5 )  {
      cerr << "obj[" << i << "] x_target \"" << check_objs[i]->x_target << "\" != \"" << objs[i]->x_target << "\"" << endl;
      exit(1);
    }

    if ( fabs ( ( check_objs[i]->y_target - objs[i]->y_target ) / objs[i]->y_target ) > 1.0e-5 )  {
      cerr << "obj[" << i << "] y_target \"" << check_objs[i]->y_target << "\" != \"" << objs[i]->y_target << "\"" << endl;
      exit(1);
    }

    if ( fabs ( ( check_objs[i]->x_fvcobs - objs[i]->x_fvcobs ) / objs[i]->x_fvcobs ) > 1.0e-5 )  {
      cerr << "obj[" << i << "] x_fvcobs \"" << check_objs[i]->x_fvcobs << "\" != \"" << objs[i]->x_fvcobs << "\"" << endl;
      exit(1);
    }

    if ( fabs ( ( check_objs[i]->y_fvcobs - objs[i]->y_fvcobs ) / objs[i]->y_fvcobs ) > 1.0e-5 )  {
      cerr << "obj[" << i << "] y_fvcobs \"" << check_objs[i]->y_fvcobs << "\" != \"" << objs[i]->y_fvcobs << "\"" << endl;
      exit(1);
    }

    if ( fabs ( ( check_objs[i]->x_fvcerr - objs[i]->x_fvcerr ) / objs[i]->x_fvcerr ) > 1.0e-5 )  {
      cerr << "obj[" << i << "] x_fvcerr \"" << check_objs[i]->x_fvcerr << "\" != \"" << objs[i]->x_fvcerr << "\"" << endl;
      exit(1);
    }

    if ( fabs ( ( check_objs[i]->y_fvcerr - objs[i]->y_fvcerr ) / objs[i]->y_fvcerr ) > 1.0e-5 )  {
      cerr << "obj[" << i << "] y_fvcerr \"" << check_objs[i]->y_fvcerr << "\" != \"" << objs[i]->y_fvcerr << "\"" << endl;
      exit(1);
    }


  }


  // write, read spec

  string desispec = datadir + "/desi_check_spec.fits.out";

  boost::property_tree::ptree child;

  props.clear();

  child.clear();
  child.put ( "VAL", 8000.0 );
  child.put ( "TYPE", "F" );
  child.put ( "COM", "Starting wavelength [Angstroms]" );
  props.put_child ( "CRVAL1", child );

  child.clear();
  child.put ( "VAL", 0.6 );
  child.put ( "TYPE", "F" );
  child.put ( "COM", "Wavelength step [Angstroms]" );
  props.put_child ( "CDELT1", child );

  child.clear();
  child.put ( "VAL", "vac" );
  child.put ( "TYPE", "C" );
  child.put ( "COM", "Vacuum wavelengths" );
  props.put_child ( "AIRORVAC", child );

  child.clear();
  child.put ( "VAL", 0 );
  child.put ( "TYPE", "I" );
  child.put ( "COM", "linear wavelength steps, not log10" );
  props.put_child ( "LOGLAM", child );

  child.clear();
  child.put ( "VAL", desipix );
  child.put ( "TYPE", "C" );
  child.put ( "COM", "Input simulation file" );
  props.put_child ( "SIMFILE", child );

  child.clear();
  child.put ( "VAL", "r0" );
  child.put ( "TYPE", "C" );
  child.put ( "COM", "Spectograph Camera" );
  props.put_child ( "CAMERA", child );

  child.clear();
  child.put ( "VAL", 1000.0 );
  child.put ( "TYPE", "F" );
  child.put ( "COM", "Exposure time [sec]" );
  props.put_child ( "EXPTIME", child );

  child.clear();
  child.put ( "VAL", 2.9 );
  child.put ( "TYPE", "F" );
  child.put ( "COM", "Read noise [electrons]" );
  props.put_child ( "RDNOISE", child );

  child.clear();
  child.put ( "VAL", "arc" );
  child.put ( "TYPE", "C" );
  child.put ( "COM", "Exposure type (arc, flat, science)" );
  props.put_child ( "FLAVOR", child );

  child.clear();
  child.put ( "VAL", "NULL" );
  child.put ( "TYPE", "C" );
  child.put ( "COM", "Input spectral PSF" );
  props.put_child ( "IN_PSF", child );

  child.clear();
  child.put ( "VAL", desipix );
  child.put ( "TYPE", "C" );
  child.put ( "COM", "Input image" );
  props.put_child ( "IN_IMG", child );

  size_t n_lambda = 100;
  size_t n_spec = 25;
  size_t n_bin = n_spec * n_lambda;

  vector_double spec_data ( n_bin );
  vector_double spec_invvar ( n_bin );
  vector_double spec_lambda ( n_lambda );

  for ( size_t i = 0; i < n_spec; ++i ) {
    for ( size_t j = 0; j < n_lambda; ++j ) {
      spec_data[i * n_lambda + j] = 100.0 * (double)i + (double)j;
      spec_invvar[i * n_lambda + j] = 1.0 / spec_data[i * n_lambda + j];
    }
  }

  for ( size_t j = 0; j < n_lambda; ++j ) {
    spec_lambda[j] = 8000.0 + 0.6 * (double)j;
  }

  spec_desi::write ( desispec, props, spec_data, spec_invvar, spec_lambda );

  
  props.clear();
  props.put ( "path", desispec );

  boost::shared_ptr < spec_desi > spc ( new spec_desi( props ) );

  vector_double chk_spec_data;
  vector_double chk_spec_invvar;
  vector_double chk_spec_lambda;

  spc->values( chk_spec_data );
  spc->inv_variance( chk_spec_invvar );
  spc->lambda( chk_spec_lambda );

  for ( size_t i = 0; i < chk_spec_data.size(); ++i ) {
    if ( fabs ( ( chk_spec_data[i] - spec_data[i] ) / spec_data[i] ) > 1.0e-8 ) {
      cerr << "FAIL on spec_data[" << i << "], " << chk_spec_data[i] << " != " << spec_data[i];
      exit(1);
    }
  }

  for ( size_t i = 0; i < chk_spec_invvar.size(); ++i ) {
    if ( fabs ( ( chk_spec_invvar[i] - spec_invvar[i] ) / spec_invvar[i] ) > 1.0e-8 ) {
      cerr << "FAIL on spec_invvar[" << i << "], " << chk_spec_invvar[i] << " != " << spec_invvar[i];
      exit(1);
    }
  }

  for ( size_t i = 0; i < chk_spec_lambda.size(); ++i ) {
    if ( fabs ( ( chk_spec_lambda[i] - spec_lambda[i] ) / spec_lambda[i] ) > 1.0e-8 ) {
      cerr << "FAIL on spec_lambda[" << i << "], " << chk_spec_lambda[i] << " != " << spec_lambda[i];
      exit(1);
    }
  }


  cout << "  (PASSED)" << endl;
     
  return;
}



