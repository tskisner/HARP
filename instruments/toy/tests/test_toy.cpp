#include <iostream>

#include <fstream>

#include <harp_test.hpp>

#include <boost/random.hpp>

using namespace std;
using namespace harp;


void toy_pcgmle_prec ( data_vec & in, data_vec & out, int_vec & flags, void * data ) {
  data_vec * prec = (data_vec *) data;
  
  data_vec :: const_iterator vit;
  for ( vit = in.begin(); vit != in.end(); ++vit ) {
    size_t pos = vit.index();
    if ( flags[ pos ] == 0 ) {
      out[ pos ] = (*prec)[pos] * (*vit);
    } else {
      out[ pos ] = 0.0;
    }
  }
  
  return;
}


void toy_pcgmle_report ( double const & norm, double const & deltazero, int const & iter, double const & alpha, double const & beta, double const & delta, double const & epsilon ) {
  
  double relerr = sqrt ( delta / deltazero );
  
  cerr << "  PCG iter " << iter << ": alpha = " << alpha << " beta = " << beta << " delta = " << delta << " epsilon = " << epsilon << " rel err = " << relerr << endl;
  
  return;
}


void toy_pcgmle_profile ( string const & name, string const & desc, double & totaltime, double & opencltime, map < string, long long int > & papi ) {
  
  cerr << "Profiling:   " << desc << ":  " << totaltime << " seconds" << endl;
  
  return;
}


void harp::test_toy ( string const & datadir ) {
  
  cerr << "Testing toy image construction..." << endl;
  
  std::map < std::string, std::string > params;
  
  params[ "path" ] = datadir + "/test_image.fits";
  params[ "signal" ] = "1";
  params[ "noise" ] = "2";
  
  image_p testimg ( image::create ( string("toy"), params ) );
  
  cerr << "  (PASSED)" << endl;
  
  
  cerr << "Testing toy image read..." << endl;
  
  dense_rowmat pix ( testimg->rows(), testimg->cols() );
  dense_rowmat noise ( testimg->rows(), testimg->cols() );
  
  dense_rowmat_view pixview ( pix, mv_range ( 0, pix.size1() ), mv_range ( 0, pix.size2() ) );
  dense_rowmat_view noiseview ( noise, mv_range ( 0, noise.size1() ), mv_range ( 0, noise.size2() ) );
  
  testimg->read ( 0, 0, pixview );
  testimg->read_noise ( 0, 0, noiseview );
  
  cerr << "  (PASSED)" << endl;
  
  
  cerr << "Testing toy spectrum construction..." << endl;
  
  params.clear();
  
  params[ "path" ] = datadir + "/test_spectra.fits";
  params[ "hdu" ] = "1";
  
  spec_p testspec ( spec::create ( string("toy"), params ) );
  
  cerr << "  (PASSED)" << endl;
  
  
  cerr << "Testing toy spectrum read..." << endl;
  
  data_vec spec ( testspec->spectrum_size( 3 ) );
  data_vec_view specview ( spec, mv_range ( 0, spec.size() ) );
  
  testspec->read_spectrum ( 3, specview );

  //for ( size_t i = 0; i < spec.size(); ++i ) {
  //  cout << spec[i] << endl;
  //}
  
  cerr << "  (PASSED)" << endl;
  
  
  cerr << "Testing toy PSF construction..." << endl;
  
  params.clear();
  
  params[ "path" ] = datadir + "/test_psf.fits";
  
  params[ "corr" ] = "20";
  
  params[ "binning" ] = "4";
  
  psf_p testpsf ( psf::create ( string("toy"), params ) );
  
  cerr << "  found " << testpsf->nspec() << " spectra" << endl;
  
  cerr << "  each with " << testpsf->specsize(0) << " flux bins" << endl;
  
  cerr << "  (PASSED)" << endl;
  
  
  cerr << "Testing toy PSF lambda read..." << endl;
  
  data_vec lambda;
  
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
  
  comp_rowmat projmat ( npix, nbins );
  
  testpsf->projection ( 0, testpsf->nspec() - 1, 0, testpsf->specsize(0) - 1, (size_t)0, testimg->cols() - 1, (size_t)0, testimg->rows() - 1, projmat );
  
  cerr << "  (PASSED)" << endl;

  
  cerr << "Testing toy PCG-MLE solve..." << endl;
  
  data_vec inspec ( nbins );
  data_vec outspec ( nbins );
  int_vec flags ( nbins );
  
  data_vec measured ( npix );
  
  for ( size_t b = 0; b < nbins; ++b ) {
    flags[b] = 0;
    //if ( ( b % 60 < 5 ) || ( b % 60 > 55 ) ) {
    //  flags[b] = 1;
    //}
    
    if ( ( b % 5 == 0 ) && ( b % testpsf->specsize(0) != 0 ) ) {
      inspec[b] = 1000.0;
    } else {
      inspec[b] = 0.0;
    }
    outspec[b] = 0.0;
  }
  
  // write input spectra
  
  
  boost::numeric::ublas::axpy_prod ( projmat, inspec, measured, true );
  
  dense_rowmat outmat ( testimg->rows(), testimg->cols() );
  
  dense_rowmat_view outview ( outmat, mv_range ( 0, testimg->rows() ), mv_range ( 0, testimg->cols() ) );
  
  size_t pixoff = 0;
  for ( size_t i = 0; i < testimg->rows(); ++i ) {
    for ( size_t j = 0; j < testimg->cols(); ++j ) {
      outmat( i, j ) = measured[pixoff];
      ++pixoff;
    }
  }
  
  
  params.clear();
  
  params[ "signal" ] = "1";
  params[ "noise" ] = "2";
  
  std::ostringstream o;
  o << testimg->rows();
  params[ "rows" ] = o.str();
  
  o.str("");
  o << testimg->cols();
  params[ "cols" ] = o.str();
  
  image_p outsigimage ( image::create ( string("toy"), params ) );
  
  outsigimage->write ( "!" + datadir + "/toy_MLE_inputs.fits.out", 0, 0, outview );
  
  
  typedef boost::ecuyer1988 base_generator_type;
  base_generator_type generator(42u);
  
  data_vec rms ( npix );
  
  for ( size_t i = 0; i < npix; ++i ) {
    rms[i] = sqrt( 16.0 + measured[i] );
    
    boost::normal_distribution < double > dist ( 0.0, rms[i] );
    
    boost::variate_generator < base_generator_type&, boost::normal_distribution < double > > gauss ( generator, dist );
    
    measured[i] += gauss();
  }
  
  
  pixoff = 0;
  for ( size_t i = 0; i < testimg->rows(); ++i ) {
    for ( size_t j = 0; j < testimg->cols(); ++j ) {
      outmat( i, j ) = measured[pixoff];
      ++pixoff;
    }
  }
  
  
  params.clear();
  
  params[ "signal" ] = "4";
  
  o.str("");
  o << testimg->rows();
  params[ "rows" ] = o.str();
  
  o.str("");
  o << testimg->cols();
  params[ "cols" ] = o.str();
  
  image_p outsnimage ( image::create ( string("toy"), params ) );
  
  outsnimage->write ( datadir + "/toy_MLE_inputs.fits.out", 0, 0, outview );
  
  // construct inverse noise covariance and preconditioner
  

  comp_rowmat invnoise ( npix, npix );
  data_vec precdata ( nbins );
  
  for ( size_t i = 0; i < npix; ++i ) {
    invnoise( i, i ) = 1.0 / ( rms[i] * rms[i] );
  }
  
  for ( size_t i = 0; i < nbins; ++i ) {
    precdata[i] = 0.0;
    for ( size_t j = 0; j < npix; ++j ) {
      precdata[i] += projmat( j, i ) * projmat( j, i ) * invnoise( j, j );
    }
  }
  
  for ( size_t i = 0; i < nbins; ++i ) {
    if ( precdata[i] > 0.0 ) {
      precdata[i] = 1.0 / precdata[i];
    } else {
      precdata[i] = 1.0;
    }
  }
  
   
  prof->reg ( "PCG_PREC", "applying preconditioner" );
  prof->reg ( "PCG_PMV", "projection matrix-vector multiply" );
  prof->reg ( "PCG_NMV", "N^-1 matrix-vector multiply" );
  prof->reg ( "PCG_VEC", "vector ops time" );
  prof->reg ( "PCG_TOT", "Total PCG time" );
  
  
  data_vec rhs ( nbins );
  data_vec q ( nbins );
  data_vec r ( nbins );
  data_vec s ( nbins );
  data_vec d ( nbins );
  
  double err = moat::la::pcg_mle < comp_rowmat, comp_rowmat, data_vec, int_vec > ( true, true, projmat, invnoise, measured, outspec, q, r, s, d, flags, rhs, 100, 1.0e-12, toy_pcgmle_prec, (void*)&precdata, toy_pcgmle_report, "PCG_TOT", "PCG_VEC", "PCG_PMV", "PCG_NMV", "PCG_PREC" );
  
  prof->stop_all();
  
  prof->query ( toy_pcgmle_profile );
  
  prof->unreg ( "PCG_PREC" );
  prof->unreg ( "PCG_PMV" );
  prof->unreg ( "PCG_NMV" );
  prof->unreg ( "PCG_VEC" );
  prof->unreg ( "PCG_TOT" );
  prof->unreg ( "PCG_PSF" );
  prof->unreg ( "PCG_REMAP" );
  
  string outdata = datadir + "/toy_MLE_spectra.out";
  
  fstream out;
  out.open ( outdata.c_str(), ios::out );
  
  for ( size_t i = 0; i < nbins; ++i ) {
    if ( flags[i] == 0 ) {
      out << i << " " << inspec[i] << " " << outspec[i] << " " << sqrt(precdata[i]) << endl;
    }
  }
  
  out.close();
  
  cerr << "  (PASSED)" << endl;
  
  return;
}



