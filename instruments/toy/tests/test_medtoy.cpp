#include <iostream>

#include <fstream>


#include <harp_test.hpp>

#include <boost/random.hpp>

using namespace std;
using namespace harp;


void medtoy_pcgmle_prec ( data_vec_view & in, data_vec_view & out, int_vec_view & flags, void * data ) {
  data_vec * prec = (data_vec *) data;
  
  data_vec_view :: const_iterator vit;
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


void medtoy_pcgmle_report ( double const & norm, double const & deltazero, int const & iter, double const & alpha, double const & beta, double const & delta, double const & epsilon ) {
  
  double relerr = sqrt ( delta / deltazero );
  
  cerr << "  PCG iter " << iter << ": alpha = " << alpha << " beta = " << beta << " delta = " << delta << " epsilon = " << epsilon << " relative err = " << relerr << endl;
  
  return;
}


void medtoy_pcgmle_profile ( string const & name, string const & desc, double & totaltime, double & opencltime, map < string, long long int > & papi ) {
  
  cerr << "Profiling:   " << desc << ":  " << totaltime << " seconds" << endl;
  
  return;
}


void harp::test_medtoy ( string const & datadir ) {
  
  
  string psffile = datadir + "/psf_gauss2d_500.fits";
  
  cerr << "Testing medium toy format spectral extraction" << endl;


  cerr << "  Reading input PSF..." << endl;
  
  std::map < std::string, std::string > params;
  
  params.clear();
  
  params[ "path" ] = psffile;
  
  params[ "corr" ] = "20";
  
  params[ "binning" ] = "4";
  
  psf_p resp ( psf::create ( string("toy"), params ) );
  
  
  cerr << "  Generating input spectra..." << endl;

  size_t rows = 500;
  size_t cols = 500;
  size_t npix = rows * cols;

  size_t nspec = resp->nspec();
  size_t specbins = resp->specsize(0);
  size_t nbins = nspec * specbins;
  
  data_vec truth ( nbins );
  
  size_t b = 0;
  for ( size_t i = 0; i < nspec; ++i ) {
    for ( size_t j = 0; j < specbins; ++j ) {
      if ( ( b % 10 == 0 ) && ( b % specbins != 0 ) ) {
        truth[b] = 2000.0;
      } else {
        truth[b] = 0.0;
      }
      ++b;
    }
  }

  
  cerr << "  Sampling PSF to build sparse projection matrix..." << endl;
  
  moat::profile * prof = moat::profile::get ( );

  prof->reg ( "PCG_PSF", "build projection matrix" );
  
  comp_rowmat projmat ( npix, nbins );
  
  prof->start ( "PCG_PSF" );
  
  resp->projection ( 0, nspec - 1, 0, specbins - 1, 0, cols - 1, 0, rows - 1, projmat );
  
  prof->stop ( "PCG_PSF" );
  
  
  cerr << "  Computing projected signal image..." << endl;
  
  data_vec measured ( npix );
  
  boost::numeric::ublas::axpy_prod ( projmat, truth, measured, true );
  
  dense_rowmat tempmat ( rows, cols );
  
  dense_rowmat_view outview ( tempmat, mv_range ( 0, rows ), mv_range ( 0, cols ) );
  
  size_t pixoff = 0;
  for ( size_t i = 0; i < rows; ++i ) {
    for ( size_t j = 0; j < cols; ++j ) {
      outview( i, j ) = measured[pixoff];
      ++pixoff;
    }
  }
  
  params.clear();
  
  params[ "hdu" ] = "1";
  
  std::ostringstream o;
  o << rows;
  params[ "rows" ] = o.str();
  
  o.str("");
  o << cols;
  params[ "cols" ] = o.str();
  
  image_p outsigimage ( image::create ( string("toy"), params ) );
  
  outsigimage->write ( "!" + datadir + "/medtoy_MLE_inputs.fits.out", 0, 0, outview );
  
  
  cerr << "  Generating noise realization..." << endl;
  
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
  for ( size_t i = 0; i < rows; ++i ) {
    for ( size_t j = 0; j < cols; ++j ) {
      tempmat( i, j ) = measured[pixoff];
      ++pixoff;
    }
  }
  
  params.clear();
  
  params[ "hdu" ] = "2";
  
  o.str("");
  o << rows;
  params[ "rows" ] = o.str();
  
  o.str("");
  o << cols;
  params[ "cols" ] = o.str();
  
  image_p outsnimage ( image::create ( string("toy"), params ) );
  
  outsnimage->write ( datadir + "/medtoy_MLE_inputs.fits.out", 0, 0, outview );
  
  
  cerr << "  Computing preconditioner..." << endl;
  
  data_vec outspec ( nbins );
  int_vec flags ( nbins );
  
  for ( size_t b = 0; b < nbins; ++b ) {
    flags[b] = 0;
    outspec[b] = 0.0;
  }
  
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
    precdata[i] = 1.0 / precdata[i];
  }
  
  
  cerr << "  Solving PCG..." << endl;
   
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
  
  double err = moat::la::pcg_mle < comp_rowmat, comp_rowmat, data_vec, int_vec > ( true, true, projmat, invnoise, measured, outspec, q, r, s, d, flags, rhs, 1000, 1.0e-12, medtoy_pcgmle_prec, (void*)&precdata, medtoy_pcgmle_report, "PCG_TOT", "PCG_VEC", "PCG_PMV", "PCG_NMV", "PCG_PREC" );
  
  prof->stop_all();
  
  prof->query ( medtoy_pcgmle_profile );
  
  prof->unreg ( "PCG_PREC" );
  prof->unreg ( "PCG_PMV" );
  prof->unreg ( "PCG_NMV" );
  prof->unreg ( "PCG_VEC" );
  prof->unreg ( "PCG_TOT" );
  prof->unreg ( "PCG_PSF" );
  
  string outdata = datadir + "/medtoy_MLE_spectra.out";
  
  fstream out;
  out.open ( outdata.c_str(), ios::out );
  
  for ( size_t i = 0; i < nbins; ++i ) {
    if ( flags[i] == 0 ) {
      out << i << " " << truth[i] << " " << outspec[i] << " " << sqrt(precdata[i]) << endl;
    }
  }
  
  out.close();
  
  cerr << "  (PASSED)" << endl;
  
  return;
}



