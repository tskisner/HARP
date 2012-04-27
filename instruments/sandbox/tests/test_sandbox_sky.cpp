#include <iostream>
#include <fstream>

#include <harp_test.hpp>

#include <boost/random.hpp>

#include <moat.hpp>

extern "C" {
#include <unistd.h>
}

using namespace std;
using namespace harp;





void sandbox_sky_profile ( string const & name, string const & desc, double & totaltime, double & opencltime, map < string, long long int > & papi ) {
  
  cerr << "Profiling:   " << desc << ":  " << totaltime << " seconds" << endl;
  
  return;
}


void harp::test_sandbox_sky ( string const & datadir ) {
  
  cerr << "Testing sandbox sky subtraction..." << endl;
  
  boost::property_tree::ptree props;
  props.clear();
  props.put ( "path", datadir + "/test_medium_psf.fits" );
  props.put ( "corr", 10 );
  psf_p testpsf ( psf::create ( string("sandbox"), props ) );

  size_t nspec = testpsf->nspec();
  size_t specsize = testpsf->specsize(0);

  size_t nbins = nspec * specsize;

  size_t nx = 67;
  size_t ny = 105;

  size_t npix = nx * ny;

  mat_compcol projmat ( npix, nbins );

  testpsf->projection ( string("SANDBOX_PSF"), string("SANDBOX_REMAP"), 0, nspec - 1, 0, specsize - 1, (size_t)0, nx - 1, (size_t)0, ny - 1, projmat );

  // construct fake spectra

  // for now, use single component...
  vec_dense truesky ( nbins );
  for ( size_t i = 0; i < nspec; ++i ) {
    for ( size_t j = 0; j < specsize; ++j ) {
      if ( (j + 4) % 20 == 0 ) {
        truesky[ specsize * i + j ] = 20000.0;
      } else {
        truesky[ specsize * i + j ] = 1000.0;
      }
    }
  }

  // we use 3 spectra as "sky fibers"
  vec_dense trueobj ( nbins );
  for ( size_t i = 0; i < nspec; ++i ) {
    for ( size_t j = 0; j < specsize; ++j ) {
      if ( (i + 1) % 3 == 0 ) {
        trueobj[ specsize * i + j ] = 0.0;
      } else {
        if ( (j + 1) % (10 + i) == 0 ) {
          trueobj[ specsize * i + j ] = 2000.0;
        } else {
          trueobj[ specsize * i + j ] = 0.0;
        }
      }
    }
  }

  vec_dense truespec ( nbins );
  for ( size_t i = 0; i < nspec; ++i ) {
    for ( size_t j = 0; j < specsize; ++j ) {
      truespec[ specsize * i + j ] = truesky[ specsize * i + j ] + trueobj[ specsize * i + j ];
    }
  }

  // compute noiseless image

  cerr << "  compute signal image..." << endl;

  vec_dense noiseless ( npix );

  boost::numeric::ublas::axpy_prod ( projmat, truespec, noiseless, true );

  // compute noise realization

  cerr << "  generate pixel noise..." << endl;

  vec_dense imgnoise ( npix );
  vec_dense measured ( npix );

  typedef boost::ecuyer1988 base_generator_type;
  base_generator_type generator(42u);
  
  vec_dense rms ( npix );
  
  for ( size_t i = 0; i < npix; ++i ) {
    rms[i] = sqrt( 16.0 + noiseless[i] );
    
    boost::normal_distribution < double > dist ( 0.0, rms[i] );
    
    boost::variate_generator < base_generator_type&, boost::normal_distribution < double > > gauss ( generator, dist );
    
    imgnoise[i] = gauss();

    measured[i] = noiseless[i] + imgnoise[i];
  }

  // construct inverse pixel noise covariance

  cerr << "  construct pixel inverse covariance..." << endl;

  mat_comprow invnoise ( npix, npix, npix );

  for ( size_t i = 0; i < npix; ++i ) {
    invnoise ( i, i ) = 1.0 / ( rms[i] * rms[i] );
  }

  cerr << "  compute noise weighted spectra..." << endl;

  vec_dense z ( nbins );

  noise_weighted_spec ( projmat, invnoise, measured, z );


  // construct the inverse spectral covariance matrix

  mat_comprow invcov ( nbins, nbins );

  inverse_covariance ( projmat, invnoise, invcov );


  // extraction

  vec_dense outspec ( nbins );

  cerr << "  extract..." << endl;

  extract_dense ( invcov, z, outspec );

  // write output  

  string outdata = datadir + "/test_skysub_output.out";

  fstream out;
  out.open ( outdata.c_str(), ios::out );
  
  for ( size_t i = 0; i < nbins; ++i ) {
    out << i << " " << truespec[i] << " " << outspec[i] << " " << 1.0 / sqrt( invcov(i,i) )<< endl;
  }
  
  out.close();

  // compute PCA of sky fibers as starting point for sky subtraction

  // we use 3 spectra as "sky fibers"

  mat_denserow skydata ( 3, specsize );
  
  size_t k = 0;
  for ( size_t i = 0; i < nspec; ++i ) {
    for ( size_t j = 0; j < specsize; ++j ) {
      if ( (i + 1) % 3 == 0 ) {
        skydata ( j, k ) = outspec[ specsize * i + j ];
        ++k;
      }
    }
  }

  boost::numeric::ublas::SingularValueDecomposition < double > svd ( skydata, false, true, true );

  mat_denserow pcomp = svd.getV();

  // restrict to just 2 components and use as the starting point for iterative
  // sky subtraction.  build component-appended design matrix.

  typedef boost::numeric::ublas::matrix_range < mat_denserow > mdr_view;

  mdr_view skycomp ( pcomp, mv_range ( 0, 2 ), mv_range ( 0, specsize ) );

  mat_compcol projmat_skycomp;

  append_sky_comp ( projmat, skycomp, nspec, projmat_skycomp );


  // do extraction to 




  /*

  outdata = datadir + "/test_skysub_est-sky0.out";

  out.open ( outdata.c_str(), ios::out );

  for ( size_t j = 0; j < specsize; ++j ) {
    meansky[ j ] /= meanskyhits[ j ];
    out << j << " " << truesky[j] << " " << meansky[j] << endl;
  }

  out.close();

  outdata = datadir + "/test_skysub_est-obj0.out";

  out.open ( outdata.c_str(), ios::out );

  for ( size_t i = 0; i < nspec; ++i ) {
    for ( size_t j = 0; j < specsize; ++j ) {
      size_t bin = i * specsize + j;
      out << bin << " " << trueobj[bin] << " " << outspec[bin] - meansky[j] << endl;
    }
  }

  out.close();

  // Iterative sky subtraction

  mat_dense 

  */






  cerr << "  (PASSED)" << endl;
     
  return;
}


