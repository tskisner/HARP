/*
  High Performance Astrophysical Reconstruction and Processing (HARP)
  (c) 2014-2015, The Regents of the University of California, 
  through Lawrence Berkeley National Laboratory.  See top
  level LICENSE file for details.
*/

#include <harp_mpi_internal.hpp>

using namespace std;
using namespace harp;



// the chunk size specified is the desired *output* size.  The input will have the additional overlaps factored in.

harp::mpi_spec_slice::mpi_spec_slice ( boost::mpi::communicator const & rcomm, boost::mpi::communicator const & gcomm, size_t nspec, size_t nlambda, size_t chunk_nspec, size_t chunk_nlambda, size_t overlap_spec, size_t overlap_lambda ) {

  rcomm_ = rcomm;
  gcomm_ = gcomm;

  local_.reset ( new spec_slice ( rcomm_.size(), nspec, nlambda, chunk_nspec, chunk_nlambda, overlap_spec, overlap_lambda ) );

}


std::vector < spec_slice_region > harp::mpi_spec_slice::regions ( ) const {
  return local_->regions ( rcomm_.rank() );
}


