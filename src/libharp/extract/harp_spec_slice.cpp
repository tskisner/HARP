/*
  High Performance Astrophysical Reconstruction and Processing (HARP)
  (c) 2014-2015, The Regents of the University of California, 
  through Lawrence Berkeley National Laboratory.  See top
  level LICENSE file for details.
*/

#include <harp_internal.hpp>

using namespace std;
using namespace harp;


spec_slice_region const & harp::spec_slice::full_region ( ) const {
  return full_region_;
}

// helper function for 1-D decomposition

void harp::spec_slice::calc ( size_t n, size_t chunk, size_t overlap, size_t offset, vector < size_t > & start, vector < size_t > & stop, vector < size_t > & good_start, vector < size_t > & good_stop ) {

  if ( offset < overlap ) {
    HARP_THROW( "for slice calculation, the offset must be at least as much as the overlap" );
  }

  start.clear();
  stop.clear();
  good_start.clear();
  good_stop.clear();

  size_t nchunk = (size_t) ( n / chunk );

  size_t chunkstart;
  size_t chunkstop;

  //cout << "calc:  " << n << " / " << chunk << " off = " << offset << " overlap = " << overlap << endl;

  for ( size_t i = 0; i < nchunk; ++i ) {
    chunkstart = i * chunk + offset;
    chunkstop = (i + 1) * chunk - 1 + offset;

    good_start.push_back ( chunkstart );
    good_stop.push_back ( chunkstop );

    start.push_back ( chunkstart - overlap );
    stop.push_back ( chunkstop + overlap );
    //cout << "calc:    " << i << ":  [" << (chunkstart - overlap) << " (" << chunkstart << " " << chunkstop << ") " << (chunkstop + overlap) << "]" << endl;
  }

  // do we need an extra chunk at the end?

  if ( n > nchunk * chunk ) {
    chunkstart = nchunk * chunk + offset;
    chunkstop = n - 1 + offset;

    good_start.push_back ( chunkstart );
    good_stop.push_back ( chunkstop );
    start.push_back ( chunkstart - overlap );
    stop.push_back ( chunkstop + overlap );

    //cout << "calc:    " << nchunk << ":  [" << (chunkstart - overlap) << " (" << chunkstart << " " << chunkstop << ") " << (chunkstop + overlap) << "]" << endl;
  } 

  return;
}


// the chunk size specified is the desired *output* size.  The input will have the additional overlaps factored in.

harp::spec_slice::spec_slice ( size_t nworker, size_t first_spec, size_t first_lambda, size_t nspec, size_t nlambda, size_t chunk_nspec, size_t chunk_nlambda, size_t overlap_spec, size_t overlap_lambda ) {

  nworker_ = nworker;
  overlap_spec_ = overlap_spec;
  overlap_lambda_ = overlap_lambda;
  first_spec_ = first_spec;
  first_lambda_ = first_lambda;
  nspec_ = nspec;
  nlambda_ = nlambda;
  chunk_nspec_ = chunk_nspec;
  chunk_nlambda_ = chunk_nlambda;

  // Determine spectral chunk boundaries

  vector < size_t > spec_start;
  vector < size_t > spec_stop;
  vector < size_t > good_spec_start;
  vector < size_t > good_spec_stop;

  calc ( nspec_, chunk_nspec_, overlap_spec_, first_spec_, spec_start, spec_stop, good_spec_start, good_spec_stop );

  size_t nchunk_spec = spec_start.size();

  // Determine lambda chunk boundaries

  vector < size_t > lambda_start;
  vector < size_t > lambda_stop;
  vector < size_t > good_lambda_start;
  vector < size_t > good_lambda_stop;

  calc ( nlambda_, chunk_nlambda_, overlap_lambda_, first_lambda_, lambda_start, lambda_stop, good_lambda_start, good_lambda_stop );

  size_t nchunk_lambda = lambda_start.size();

  // Define a region which contains all bins

  full_region_.overlap_spec = 0;
  full_region_.overlap_lambda = 0;
  full_region_.first_spec = spec_start[0];
  full_region_.first_lambda = lambda_start[0];
  full_region_.first_good_spec = full_region_.first_spec;
  full_region_.first_good_lambda = full_region_.first_lambda;
  full_region_.n_spec = spec_stop[nchunk_spec - 1] - full_region_.first_spec + 1;
  full_region_.n_lambda = lambda_stop[nchunk_lambda - 1] - full_region_.first_lambda + 1;
  full_region_.n_good_spec = full_region_.n_spec;
  full_region_.n_good_lambda = full_region_.n_lambda;

  // Assign ranges of chunks to workers

  size_t total_chunks = nchunk_spec * nchunk_lambda;

  size_t min = (size_t)( total_chunks / nworker_ );
  size_t leftover = total_chunks % nworker_;

  vector < size_t > worker_first ( nworker_ );
  vector < size_t > worker_n ( nworker_ );

  for ( size_t w = 0; w < nworker_; ++w ) {
    worker_n [ w ] = min;

    if ( w < leftover ) {
      ++worker_n [ w ];
      worker_first [ w ] = w * worker_n [ w ];
    } else {
      worker_first [ w ] = ( (worker_n[w] + 1) * leftover ) + ( worker_n[w] * (w - leftover) );
    }
  }

  // for each worker, construct the assigned regions

  for ( size_t w = 0; w < nworker_; ++w ) {

    //cerr << "DBG:  worker " << w << " has " << worker_n[w] << " chunks" << endl;

    for ( size_t chunk = 0; chunk < worker_n[w]; ++chunk ) {

      size_t abs_chunk = worker_first[w] + chunk;

      // chunks are distributed in wavelength-major order

      size_t abs_lambda = (size_t)( abs_chunk / nchunk_spec );
      size_t abs_spec = abs_chunk - ( abs_lambda * nchunk_spec );

      spec_slice_region reg;
      reg.overlap_spec = overlap_spec_;
      reg.overlap_lambda = overlap_lambda_;
      reg.first_spec = spec_start[ abs_spec ];
      reg.first_lambda = lambda_start[ abs_lambda ];
      reg.first_good_spec = good_spec_start[ abs_spec ];
      reg.first_good_lambda = good_lambda_start[ abs_lambda ];
      reg.n_spec = spec_stop[ abs_spec ] - reg.first_spec + 1;
      reg.n_lambda = lambda_stop[ abs_lambda ] - reg.first_lambda + 1;
      reg.n_good_spec = good_spec_stop[ abs_spec ] - reg.first_good_spec + 1;
      reg.n_good_lambda = good_lambda_stop[ abs_lambda ] - reg.first_good_lambda + 1;

      (regions_[ w ]).push_back ( reg );

      /*
      cerr << "DBG:  worker " << w << " global chunk " << abs_chunk << " :" << endl;
      cerr << "DBG:     n_spec = " << regions_[w][ regions_[w].size() - 1 ].n_spec << endl;
      cerr << "DBG:     n_good_spec = " << regions_[w][ regions_[w].size() - 1 ].n_good_spec << endl;
      cerr << "DBG:     first_spec = " << regions_[w][ regions_[w].size() - 1 ].first_spec << endl;
      cerr << "DBG:     first_good_spec = " << regions_[w][ regions_[w].size() - 1 ].first_good_spec << endl;
      cerr << "DBG:     overlap_spec = " << regions_[w][ regions_[w].size() - 1 ].overlap_spec << endl;
      cerr << "DBG:     n_lambda = " << regions_[w][ regions_[w].size() - 1 ].n_lambda << endl;
      cerr << "DBG:     n_good_lambda = " << regions_[w][ regions_[w].size() - 1 ].n_good_lambda << endl;
      cerr << "DBG:     first_lambda = " << regions_[w][ regions_[w].size() - 1 ].first_lambda << endl;
      cerr << "DBG:     first_good_lambda = " << regions_[w][ regions_[w].size() - 1 ].first_good_lambda << endl;
      cerr << "DBG:     overlap_lambda = " << regions_[w][ regions_[w].size() - 1 ].overlap_lambda << endl;
      */

    }

  }

}


std::vector < spec_slice_region > harp::spec_slice::regions ( size_t const & worker ) const {
  
  if ( worker > nworker_ - 1 ) {
    HARP_THROW( "worker rank is out of range" );
  }

  std::map < size_t, std::vector < spec_slice_region > > :: const_iterator rit = regions_.find ( worker );

  if ( rit == regions_.end() ) {
    return std::vector < spec_slice_region > ();
  } else {
    return rit->second;
  }
}


