#include <iostream>

#include <fstream>

#include <harp_test.hpp>

#include <boost/random.hpp>

extern "C" {
#include <unistd.h>
}

using namespace std;
using namespace harp;



extern "C" void dgesdd_ ( char *, int *, int *, double *, int *, double *, double *, int *, double *, int *, double *, int *, int *, int * );

extern "C" void dsyevd_ ( char *, char *, int *, double *, int *, double *, double *, int *, int *, int *, int * );


void lapack_svd ( int dim, double * mat, double * eigen ) {
  char jobz = 'O';
  int lda = dim;
  int ldu = dim;
  double * U = NULL;
  double * VT = moat::double_alloc ( dim * dim );
  int lwork;
  double * work;
  int * iwork = moat::int_alloc ( 8 * dim );
  int info;

  lwork = dim * ( 5 * dim + 7 );
  work = moat::double_alloc ( lwork );

  dgesdd_ ( &jobz, &dim, &dim, mat, &lda, eigen, U, &ldu, VT, &dim, work, &lwork, iwork, &info );

  free ( iwork );
  free ( work );
  free ( VT );

  return;
}


int lapack_ev ( int dim, double * mat, double * eigen ) {
  char jobz = 'V';
  char uplo = 'L';
  int lwork;
  int liwork;
  double * work;
  int * iwork;
  int info;

  lwork = 1 + 6 * dim + 2 * dim * dim;
  work = moat::double_alloc ( lwork );

  liwork = 3 + 5 * dim;
  iwork = moat::int_alloc ( liwork );

  dsyevd_ ( &jobz, &uplo, &dim, mat, &dim, eigen, work, &lwork, iwork, &liwork, &info );

  free ( iwork );
  free ( work );

  return info;
}




void toy_profile ( string const & name, string const & desc, double & totaltime, double & opencltime, map < string, long long int > & papi ) {
  
  cerr << "Profiling:   " << desc << ":  " << totaltime << " seconds" << endl;
  
  return;
}


void harp::test_toy ( string const & datadir ) {
  
  cerr << "Testing toy image construction..." << endl;
  
  std::map < std::string, std::string > params;
  params[ "path" ] = datadir + "/test_medium_input_image.fits";
  params[ "signal" ] = "1";
  params[ "noise" ] = "2";
  image_p testimg ( image::create ( string("toy"), params ) );
  
  cerr << "  (PASSED)" << endl;
  
  
  cerr << "Testing toy image read..." << endl;
  
  mat_denserow pix ( testimg->rows(), testimg->cols() );
  mat_denserow noise ( testimg->rows(), testimg->cols() );
  
  testimg->read ( 0, 0, pix );
  testimg->read_noise ( 0, 0, noise );
  
  cerr << "  (PASSED)" << endl;
  
  
  cerr << "Testing toy spectrum construction..." << endl;
  
  params.clear();
  params[ "path" ] = datadir + "/test_medium_input_spectra.fits";
  params[ "hdu" ] = "1";
  spec_p testspec ( spec::create ( string("toy"), params ) );
  
  cerr << "  (PASSED)" << endl;
  
  
  cerr << "Testing toy spectrum read..." << endl;
  
  vec_dense spec ( testspec->spectrum_size( 3 ) );
  
  testspec->read_spectrum ( 3, spec );

  //for ( size_t i = 0; i < spec.size(); ++i ) {
  //  cout << spec[i] << endl;
  //}
  
  cerr << "  (PASSED)" << endl;
  
  
  cerr << "Testing toy PSF construction..." << endl;
  
  params.clear();
  params[ "path" ] = datadir + "/test_medium_psf.fits";
  params[ "corr" ] = "10";
  psf_p testpsf ( psf::create ( string("toy"), params ) );
  
  cerr << "  found " << testpsf->nspec() << " spectra" << endl;
  
  cerr << "  each with " << testpsf->specsize(0) << " flux bins" << endl;
  
  cerr << "  (PASSED)" << endl;
  
  
  cerr << "Testing toy PSF lambda read..." << endl;
  
  vec_dense lambda;
  
  for ( size_t i = 0; i < testpsf->nspec(); ++i ) {
    testpsf->lambda ( i, lambda );
  }
  
  cerr << "  (PASSED)" << endl;
  
  
  cerr << "Testing toy PSF calculation of sparse projection matrix..." << endl;
  
  moat::profile * prof = moat::profile::get ( );

  prof->reg ( "PCG_PSF", "compute projection matrix" );
  prof->reg ( "PCG_REMAP", "remap projection matrix" );
  
  size_t nbins = testpsf->nspec() * testpsf->specsize(0);
  size_t npix = testimg->rows() * testimg->cols();
  
  mat_compcol projmat ( npix, nbins );

  testpsf->projection ( string("PCG_PSF"), string("PCG_REMAP"), 0, testpsf->nspec() - 1, 0, testpsf->specsize(0) - 1, (size_t)0, testimg->cols() - 1, (size_t)0, testimg->rows() - 1, projmat );
  
  cerr << "  (PASSED)" << endl;

  
  cerr << "Testing toy solve..." << endl << endl;
  
  vec_dense inspec ( nbins );
  vec_dense outspec ( nbins );
  vec_flag flags ( nbins );
  
  vec_dense measured ( npix );
  vec_dense imgnoise ( npix );

  
  // read input spectrum

  testspec->read ( inspec );

  // read image data and noise

  testimg->read ( measured );
  testimg->read_noise ( imgnoise );
  
  for ( size_t b = 0; b < nbins; ++b ) {
    flags[b] = 0;
    /*
    //if ( ( b % 60 < 5 ) || ( b % 60 > 55 ) ) {
    //  flags[b] = 1;
    //}
    
    if ( ( b % 5 == 0 ) && ( b % testpsf->specsize(0) != 0 ) ) {
      inspec[b] = 2000.0;
    } else {
      inspec[b] = 0.0;
    }
    */
    outspec[b] = 0.0;
  }

  // write input spectra
  /*
  cerr << "Write input spectra..." << endl;

  std::ostringstream o;

  string delpath = datadir + "/test_input_spectra.fits.out";
  int ret = unlink ( delpath.c_str() );

  params.clear();
  params[ "hdu" ] = "1";
  o.str("");
  o << testpsf->nspec();
  params[ "nspec" ] = o.str();
  o.str("");
  o << testpsf->specsize(0);
  params[ "specsize" ] = o.str();
  spec_p specinput ( spec::create ( string("toy"), params ) );

  specinput->write ( datadir + "/test_input_spectra.fits.out", inspec );

  cerr << "  (PASSED)" << endl;

  // construct signal image

  cerr << "Construct input noisy image..." << endl;

  boost::numeric::ublas::axpy_prod ( projmat, inspec, measured, true );

  // add noise to image

  typedef boost::ecuyer1988 base_generator_type;
  base_generator_type generator(42u);
  
  vec_dense rms ( npix );
  
  for ( size_t i = 0; i < npix; ++i ) {
    //rms[i] = sqrt( 16.0 + measured[i] );
    rms[i] = sqrt( measured[i] / 16.0 );
    
    boost::normal_distribution < double > dist ( 0.0, rms[i] );
    
    boost::variate_generator < base_generator_type&, boost::normal_distribution < double > > gauss ( generator, dist );
    
    measured[i] += gauss();
    imgnoise[i] = 1.0 / (rms[i] * rms[i]);
  }

  cerr << "  (PASSED)" << endl;

  // write out image

  cerr << "Write input image and noise variance..." << endl;

  delpath = datadir + "/test_input_image.fits.out";
  ret = unlink ( delpath.c_str() );

  params.clear();

  params[ "signal" ] = "1";
  params[ "noise" ] = "2";  
  o.str("");
  o << testimg->rows();
  params[ "rows" ] = o.str();
  o.str("");
  o << testimg->cols();
  params[ "cols" ] = o.str();
  
  image_p outimage ( image::create ( string("toy"), params ) );
  
  outimage->write ( datadir + "/test_input_image.fits.out", measured );
  outimage->write_noise ( datadir + "/test_input_image.fits.out", imgnoise );

  cerr << "  (PASSED)" << endl;
  */

  
  // construct inverse pixel noise covariance

  cerr << "Construct pixel noise covariance..." << endl;
  
  mat_comprow invnoise ( npix, npix, npix );

  for ( size_t i = 0; i < npix; ++i ) {
    invnoise ( i, i ) = imgnoise[i];
  }
  
  cerr << "  (PASSED)" << endl;


  // construct the inverse spectral covariance matrix

  mat_comprow invcov ( nbins, nbins );

  mat_compcol temp ( npix, nbins );

  cerr << "Multiply N^-1 A..." << endl;

  boost::numeric::ublas::axpy_prod ( invnoise, projmat, temp, true );

  cerr << "  (PASSED)" << endl;


  cerr << "Multiply A^T (N^-1 A)..." << endl;

  boost::numeric::ublas::axpy_prod ( boost::numeric::ublas::trans ( projmat ), temp, invcov, true );

  cerr << "  (PASSED)" << endl;


  // remap to row-major?

  size_t i, j, k;



  
  // SVD
  
  cerr << "SVD" << endl;

  vec_dense eigenvals ( nbins );
  vector < vec_dense > eigenvecs ( nbins );

  vec_dense finvcov ( nbins * nbins );
  for ( j = 0; j < nbins; ++j ) {
    for ( i = 0; i < nbins; ++i ) {
      finvcov[ j * nbins + i ] = invcov ( i, j );
    }
  }

  int fret = lapack_ev ( (int)nbins, &(finvcov(0)), &(eigenvals[0]) );

  mat_densecol W ( nbins, nbins );
  for ( j = 0; j < nbins; ++j ) {
    for ( i = 0; i < nbins; ++i ) {
      W ( j, i ) = finvcov[ j * nbins + i ];
    }
  }


  string outdata = datadir + "/test_output_W.out";

  fstream out;
  out.open ( outdata.c_str(), ios::out );   

  for ( i = 0; i < nbins; ++i ) {
    for ( j = 0; j < 6; ++j ) {
      out << W(i,j) << " ";
    }
    out << " .... " << endl;
  }

  out.close();


  size_t neigen = (size_t)( nbins );

  cerr << "Condition number = " << eigenvals[neigen-1] / eigenvals[0] << endl;

  for ( i = 0; i < nbins; ++i ) {
    cerr << i << ": " << eigenvals[i] << " => " << W(0,i) << " " << W(1,i) << " " << W(2,i) << " ..." << endl;
  }

  /*
  svd < mat_comprow > ( invcov, eigenvals, eigenvecs );

  size_t firsteigen = eigenvecs.size() - eigenvals.size();
  size_t neigen = eigenvals.size() - firsteigen;
  */
  

  cerr << "  (PASSED)" << endl;


  // noise weighted spectra

  cerr << "compute Z" << endl;

  vec_dense z ( nbins );
  vec_dense ztemp ( npix );

  boost::numeric::ublas::axpy_prod ( invnoise, measured, ztemp, true );

  boost::numeric::ublas::axpy_prod ( ztemp, projmat, z, true ); 
  
  cerr << "  (PASSED)" << endl; 


  // compute norm vector

  cerr << "compute norm" << endl;

  vec_dense norm ( nbins );

  vec_dense sqrtvals ( nbins );

  moat::sf::sqrt ( neigen, &(eigenvals[0]), &(sqrtvals[0]) );

  mat_comprow diag ( nbins, nbins, nbins );
  for ( i = 0; i < neigen; ++i ) {
    diag ( i, i ) = sqrtvals[i];
  }
  for ( i = neigen; i < nbins; ++i ) {
    diag ( i, i ) = 0.0;
  }

  vec_dense vtemp1 ( nbins );
  vec_dense vtemp2 ( nbins );
  mat_denserow res ( nbins, nbins );
  mat_denserow cov ( nbins, nbins );
  mat_denserow mtemp1 ( nbins, nbins );
  //mat_denserow mtemp2 ( nbins, nbins );
  vec_dense rc ( nbins );
  vec_dense ones ( nbins );

  for ( i = 0; i < nbins; ++i ) {
    ones[i] = 1.0;
  }

  boost::numeric::ublas::axpy_prod ( W, ones, vtemp1, true );

  boost::numeric::ublas::axpy_prod ( diag, vtemp1, vtemp2, true );

  boost::numeric::ublas::axpy_prod ( vtemp2, W, norm, true );
  
  outdata = datadir + "/test_output_norm.out";

  out.open ( outdata.c_str(), ios::out );

  for ( i = 0; i < neigen; ++i ) {
    out << eigenvals[i] << " " << norm[i] << endl;
  }

  out.close(); 

  cerr << "  (PASSED)" << endl;

  // compute R explicitly

  boost::numeric::ublas::axpy_prod ( diag, W, mtemp1, true );

  boost::numeric::ublas::axpy_prod ( boost::numeric::ublas::trans ( W ), mtemp1, res, true );

  outdata = datadir + "/test_output_R.out";

  out.open ( outdata.c_str(), ios::out );   

  for ( i = 0; i < nbins; ++i ) {
    for ( j = 0; j < nbins; ++j ) {
      res ( i, j ) /= norm[i];
      if ( j < 6 ) {
        out << res(i,j) << " ";
      }
    }
    out << " .... " << endl;
  }

  out.close();

  vec_dense Rinspec ( nbins );

  boost::numeric::ublas::axpy_prod ( res, inspec, Rinspec, true );

  // compute C explicity

  for ( i = 0; i < neigen; ++i ) {
    diag ( i, i ) = 1.0 / eigenvals[i];
  }
  for ( i = neigen; i < nbins; ++i ) {
    diag ( i, i ) = 0.0;
  }

  boost::numeric::ublas::axpy_prod ( diag, W, mtemp1, true );

  boost::numeric::ublas::axpy_prod ( boost::numeric::ublas::trans ( W ), mtemp1, cov, true );

  outdata = datadir + "/test_output_cov.out";

  out.open ( outdata.c_str(), ios::out );

  for ( i = 0; i < nbins; ++i ) {
    for ( j = 0; j < 6; ++j ) {
      out << cov(i,j) << " ";
    }
    out << " .... " << endl;
  }

  out.close();


  /*
  // compute the RC product

  cerr << "compute RC" << endl;

  for ( i = 0; i < neigen; ++i ) {
    //diag ( i, i ) = 1.0 / eigenvals[i];
    diag ( i, i ) = 1.0 / sqrt ( eigenvals[i] );
  }
  for ( i = neigen; i < nbins; ++i ) {
    diag ( i, i ) = 0.0;
  }

  boost::numeric::ublas::axpy_prod ( finvcov, z, vtemp1, true );

  boost::numeric::ublas::axpy_prod ( diag, vtemp1, vtemp2, true );

  boost::numeric::ublas::axpy_prod ( vtemp2, finvcov, rc, true );
  */

  cerr << "  (PASSED)" << endl;


  cerr << "final spectra" << endl;

  boost::numeric::ublas::axpy_prod ( res, cov, mtemp1, true );

  boost::numeric::ublas::axpy_prod ( mtemp1, z, outspec, true );

  outdata = datadir + "/test_output_spectra.out";

  out.open ( outdata.c_str(), ios::out );
  
  for ( size_t i = 0; i < nbins; ++i ) {
    //outspec[i] = rc[i] / norm[i];
    //outspec[i] = rc[i];
    if ( flags[i] == 0 ) {
      out << i << " " << Rinspec[i] << " " << outspec[i] << " " << 1.0/norm[i] << " " << norm[i] * (Rinspec[i] - outspec[i]) << endl;
    }
  }
  
  out.close();


  cerr << "  (PASSED)" << endl;
     
  return;
}



