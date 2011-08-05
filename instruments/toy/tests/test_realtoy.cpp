#include <iostream>

#include <fstream>


#include <harp_test.hpp>

#include <boost/random.hpp>

using namespace std;
using namespace harp;


void realtoy_pcgmle_prec ( data_vec & in, data_vec & out, int_vec & flags, void * data ) {
  data_vec * prec = (data_vec *) data;
  
  size_t vit;
  size_t n = in.size();
  
  #ifdef _OPENMP
  #pragma omp parallel for default(none) private(vit) shared(n, flags, in, out, prec) schedule(static)
  #endif
  for ( vit = 0; vit < n; ++vit ) {
    if ( flags[ vit ] == 0 ) {
      out[ vit ] = (*prec)[ vit ] * in[ vit ];
    } else {
      out[ vit ] = 0.0;
    }
  }
  
  return;
}


void realtoy_pcgmle_report ( double const & norm, double const & deltazero, int const & iter, double const & alpha, double const & beta, double const & delta, double const & epsilon ) {
  
  double relerr = sqrt ( delta / deltazero );
  
  cerr << "  PCG iter " << iter << ": alpha = " << alpha << " beta = " << beta << " delta = " << delta << " epsilon = " << epsilon << " relative err = " << relerr << endl;
  
  return;
}


void realtoy_pcgmle_profile ( string const & name, string const & desc, double & totaltime, double & opencltime, map < string, long long int > & papi ) {
  
  cerr << "Profiling:   " << desc << ":  " << totaltime << " seconds" << endl;
  
  return;
}


void harp::test_realtoy ( string const & datadir ) {
  
  moat::profile * prof = moat::profile::get ( );
  
  prof->reg ( "PCG_PSF", "compute projection matrix" );
  prof->reg ( "PCG_REMAP", "remap projection matrix" );
  prof->reg ( "PCG_PRECALC", "compute preconditioner" );

  prof->reg ( "PCG_PREC", "applying preconditioner" );
  prof->reg ( "PCG_PMV", "projection matrix-vector multiply" );
  prof->reg ( "PCG_NMV", "N^-1 matrix-vector multiply" );
  prof->reg ( "PCG_VEC", "vector ops time" );
  prof->reg ( "PCG_TOT", "Total PCG time" );
  
  string psffile = datadir + "/realtoy_psf.fits";
  
  cerr << "Testing realistic toy format spectral extraction" << endl;


  cerr << "  Reading input PSF..." << endl;
  
  std::map < std::string, std::string > params;
  
  params.clear();
  
  params[ "path" ] = psffile;
  
  params[ "corr" ] = "10";
  
  params[ "binning" ] = "1";
  
  psf_p resp ( psf::create ( string("toy"), params ) );
  
  size_t nspec = resp->nspec();
  size_t specbins = resp->specsize(0);
  size_t nbins = nspec * specbins;
  
  
  cerr << "  Reading input true spectra..." << endl;
  
  data_vec truth ( nbins );
  
  for ( size_t s = 0; s < nspec; ++s ) {
  
    params.clear();
  
    params[ "path" ] = datadir + "/realtoy_input_spectra.fits";
    params[ "hdu" ] = "1";
    
    ostringstream o;
    o << s;
    params[ "pos" ] = o.str();
  
    spectrum_p truthspec ( spectrum::create ( string("toy"), params ) );
  
    data_vec_view truthview ( truth, mv_range ( s * specbins, (s+1) * specbins ) );
  
    truthspec->read ( truthview );
  
  }
  
  
  cerr << "  Reading data image..." << endl;
  
  params.clear();
  
  params[ "path" ] = datadir + "/realtoy_input_image.fits";
  params[ "hdu" ] = "1";
  
  image_p dataimg ( image::create ( string("toy"), params ) );
  
  dense_rowmat img ( dataimg->rows(), dataimg->cols() );
  
  dense_rowmat_view imgview ( img, mv_range ( 0, img.size1() ), mv_range ( 0, img.size2() ) );
  
  dataimg->read ( 0, 0, imgview );
  
  size_t rows = dataimg->rows();
  size_t cols = dataimg->cols();
  size_t npix = rows * cols;
  
  
  cerr << "  Reading pure signal image..." << endl;
  
  params.clear();
  
  params[ "path" ] = datadir + "/realtoy_input_image.fits";
  params[ "hdu" ] = "3";
  
  image_p sigimg ( image::create ( string("toy"), params ) );
  
  dense_rowmat sig ( sigimg->rows(), sigimg->cols() );
  
  dense_rowmat_view sigview ( sig, mv_range ( 0, sig.size1() ), mv_range ( 0, sig.size2() ) );
  
  sigimg->read ( 0, 0, sigview );
  
  data_vec rms ( npix );
  data_vec measured ( npix );
  
  size_t pixoff = 0;
  for ( size_t i = 0; i < rows; ++i ) {
    for ( size_t j = 0; j < cols; ++j ) {
      rms[i*cols + j] = sqrt( 16.0 + sig(i,j) );
      measured[pixoff] = img(i,j);
      ++pixoff;
    }
  }
  
  
  cerr << "  Reading noise covariance image..." << endl;
  
  params.clear();
  
  params[ "path" ] = datadir + "/realtoy_input_image.fits";
  params[ "hdu" ] = "2";
  
  image_p nseimg ( image::create ( string("toy"), params ) );
  
  dense_rowmat nse ( sigimg->rows(), sigimg->cols() );
  
  dense_rowmat_view nseview ( nse, mv_range ( 0, nse.size1() ), mv_range ( 0, nse.size2() ) );
  
  nseimg->read ( 0, 0, nseview );
  
  data_vec trueinv ( npix );
  
  pixoff = 0;
  for ( size_t i = 0; i < rows; ++i ) {
    for ( size_t j = 0; j < cols; ++j ) {
      trueinv[pixoff] = nse(i,j);
      ++pixoff;
    }
  }
  
  
  cerr << "  Computing PSF projection..." << endl;
  
  comp_rowmat projmat ( npix, nbins );
  
  resp->projection ( 0, nspec - 1, 0, specbins - 1, 0, cols - 1, 0, rows - 1, projmat );
  
  string outpsffile = datadir + "/realtoy_compare_psf.out";
  
  fstream out;
  
  /*
  out.open ( outpsffile.c_str(), ios::out );
  
  comp_rowmat :: iterator1 rowit;
  comp_rowmat :: iterator2 colit;
  
  for ( rowit = projmat.begin1(); rowit != projmat.end1(); ++rowit ) {
    for ( colit = rowit.begin(); colit != rowit.end(); ++colit ) {
      out << rowit.index1() << " " << colit.index2() << " " << (*colit) << endl;
    }
  }

  out.close();
  */
  
  
  cerr << "  Testing PSF convolution..." << endl;
  
  data_vec compare ( npix );
  
  string outtruthfile = datadir + "/realtoy_compare_truth.out";
  
  out.open ( outtruthfile.c_str(), ios::out );
  
  for ( size_t i = 0; i < truth.size(); ++i ) {
    out << i << " " << truth(i) << endl;
  }

  out.close();
  
  boost::numeric::ublas::axpy_prod ( projmat, truth, compare, true );
  
  dense_rowmat tempmat ( rows, cols );
  
  dense_rowmat_view outview ( tempmat, mv_range ( 0, rows ), mv_range ( 0, cols ) );
  
  pixoff = 0;
  for ( size_t i = 0; i < rows; ++i ) {
    for ( size_t j = 0; j < cols; ++j ) {
      outview( i, j ) = compare[pixoff];
      ++pixoff;
    }
  }
  
  params.clear();
  
  params[ "hdu" ] = "1";
  
  ostringstream o;
  o << rows;
  params[ "rows" ] = o.str();
  
  o.str("");
  o << cols;
  params[ "cols" ] = o.str();
  
  image_p outsigimage ( image::create ( string("toy"), params ) );
  
  outsigimage->write ( "!" + datadir + "/realtoy_psf_convolved.fits.out", 0, 0, outview );
  
  string outimgfile = datadir + "/realtoy_compare_sigmap.out";
  
  out.open ( outimgfile.c_str(), ios::out );
  
  double ratio;
  pixoff = 0;
  for ( size_t i = 0; i < rows; ++i ) {
    for ( size_t j = 0; j < cols; ++j ) {
      ratio = compare[pixoff] - sig(i,j);
      out << pixoff << " " << i << " " << j << " " << sig(i,j) << " " << compare[pixoff] << " " << ratio << endl;
      ++pixoff;
    }
  }

  out.close();
  
  
  
  cerr << "  Computing preconditioner..." << endl;
  
  prof->start ( "PCG_PRECALC" );
  
  data_vec outspec ( nbins );
  int_vec flags ( nbins );
  
  for ( size_t b = 0; b < nbins; ++b ) {
    flags[b] = 0;
    outspec[b] = 0.0;
  }
  
  comp_rowmat invnoise ( npix, npix );
  data_vec precdata ( nbins );
  
  for ( size_t i = 0; i < npix; ++i ) {
    invnoise( i, i ) = trueinv[i];
    //invnoise( i, i ) = 1.0 / ( rms[i] * rms[i] );
    //cout << measured[i] << endl;
  }
  
  for ( size_t i = 0; i < nbins; ++i ) {
    precdata[i] = 0.0;
    for ( size_t j = 0; j < npix; ++j ) {
      precdata[i] += projmat( j, i ) * projmat( j, i ) * invnoise( j, j );
    }
  }
  
  for ( size_t i = 0; i < nbins; ++i ) {
    precdata[i] = 1.0 / precdata[i];
    //cout << precdata[i] << endl;
  }
  
  prof->stop ( "PCG_PRECALC" );
  
  
  cerr << "  Solving PCG..." << endl;
  
  data_vec rhs ( nbins );
  data_vec q ( nbins );
  data_vec r ( nbins );
  data_vec s ( nbins );
  data_vec d ( nbins );
  
  double err = moat::la::pcg_mle < comp_rowmat, comp_rowmat, data_vec, int_vec > ( true, true, projmat, invnoise, measured, outspec, q, r, s, d, flags, rhs, 100, 1.0e-12, realtoy_pcgmle_prec, (void*)&precdata, realtoy_pcgmle_report, "PCG_TOT", "PCG_VEC", "PCG_PMV", "PCG_NMV", "PCG_PREC" );
  
  prof->stop_all();
  
  prof->query ( realtoy_pcgmle_profile );
  
  prof->unreg ( "PCG_PREC" );
  prof->unreg ( "PCG_PMV" );
  prof->unreg ( "PCG_NMV" );
  prof->unreg ( "PCG_VEC" );
  prof->unreg ( "PCG_TOT" );
  prof->unreg ( "PCG_PSF" );
  prof->unreg ( "PCG_REMAP" );
  prof->unreg ( "PCG_PRECALC" );
  
  
  string outspecfile = datadir + "/realtoy_solved_spectra.fits.out";
  
  string rmcom = "rm -f " + outspecfile;
  
  system( rmcom.c_str() );
  
  params.clear();
  
  o.str("");
  o << specbins;

  fitsfile * fp;
  fits::open_readwrite ( fp, outspecfile );
  
  fits::img_append ( fp, nspec, specbins );
  
  fits::close ( fp );


  //params[ "path" ] = outspecfile;
  params[ "hdu" ] = "1";
  params[ "size" ] = o.str();
  
  for ( size_t s = 0; s < nspec; ++s ) {
  
    o.str("");
    o << s;
    params[ "pos" ] = o.str();
  
    spectrum_p solvespec ( spectrum::create ( string("toy"), params ) );
  
    data_vec_view solvespecview ( outspec, mv_range ( s * specbins, (s+1) * specbins ) );
  
    solvespec->write ( outspecfile, solvespecview );
  
  }
  
  
  outspecfile = datadir + "/realtoy_solved_spectra.out";
  
  out.open ( outspecfile.c_str(), ios::out );
  
  for ( size_t i = 0; i < nbins; ++i ) {
    if ( flags[i] == 0 ) {
      out << i << " " << truth[i] << " " << outspec[i] << " " << sqrt(precdata[i]) << endl;
    }
  }
  
  out.close();
  
  
  cerr << "  (PASSED)" << endl;
  
  return;
}



