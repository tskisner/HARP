#include <iostream>

#include <harp_test.hpp>

#include <boost/random.hpp>

using namespace std;
using namespace harp;


void harp::test_tinyKLT ( string const & datadir ) {

  cerr << "Testing tiny KLT transform..." << endl;

  cerr.precision(16);
  
  // construct synthetic design matrix

  size_t nbins = 25;
  size_t npix = 225;
  
  mat_compcol projmat ( npix, nbins );

  // loop over traces
  for ( size_t i = 0; i < 5; ++i ) {

    // loop over bins
    for ( size_t j = 0; j < 5; ++j ) {
      size_t pixrow = 3 * j + 1;
      size_t pixcol = 3 * i + 1;

      // center pixel = 100
      projmat( 15 * pixrow + pixcol, 5 * i + j ) = 100.0;

      // 4 neigbors = 10
      projmat( 15 * (pixrow - 1) + pixcol, 5 * i + j ) = 10.0;
      projmat( 15 * (pixrow + 1) + pixcol, 5 * i + j ) = 10.0;
      projmat( 15 * pixrow + (pixcol - 1), 5 * i + j ) = 10.0;
      projmat( 15 * pixrow + (pixcol + 1), 5 * i + j ) = 10.0;
    }
  }

  // construct fake spectra

  vec_dense truespec ( nbins );
  for ( size_t i = 0; i < 5; ++i ) {
    for ( size_t j = 0; j < 5; ++j ) {
      truespec( 5 * i + j ) = (double)( 4 * (i + 1) * (j + 1) );
    }
  }

  // compute noiseless image

  vec_dense noiseless ( npix );

  boost::numeric::ublas::axpy_prod ( projmat, truespec, noiseless, true );

  // compute noise realization

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

  mat_comprow invnoise ( npix, npix, npix );

  for ( size_t i = 0; i < npix; ++i ) {
    invnoise ( i, i ) = 1.0 / ( rms[i] * rms[i] );
  }

  // construct rhs

  vec_dense z ( nbins );

  noise_weighted_spec < mat_compcol, mat_comprow, vec_dense > ( projmat, invnoise, measured, z );

  // construct the inverse spectral covariance matrix

  mat_comprow invcov ( nbins, nbins );

  mat_dynrow builder ( nbins, nbins );

  mat_compcol temp ( npix, nbins );

  boost::numeric::ublas::axpy_prod ( invnoise, projmat, temp, true );

  boost::numeric::ublas::axpy_prod ( boost::numeric::ublas::trans ( projmat ), temp, builder, true );

  mat_dynrow::iterator2 itcol;
  mat_dynrow::iterator1 itrow;

  for ( itcol = builder.begin2(); itcol != builder.end2(); ++itcol ) {
    for ( itrow = itcol.begin(); itrow != itcol.end(); ++itrow ) {
      invcov ( itrow.index1(), itrow.index2() ) = (*itrow);
    }
  }
  
  // extraction

  vec_dense outspec ( nbins );

  extract_dense < mat_comprow > ( invcov, z, outspec );

  for ( size_t i = 0; i < 5; ++i ) {
    for ( size_t j = 0; j < 5; ++j ) {
      cerr << "(" << i << "," << j << ") " << truespec( 5 * i + j ) << " " << outspec( 5 * i + j) << " +/- " << 1.0 / sqrt( invcov(5*i+j, 5*i+j) ) << endl;
    }
  }
  
  cerr << "  (PASSED)" << endl;

  return;
}

