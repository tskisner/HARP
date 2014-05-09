// @COPYRIGHT@

#include <harp_mpi_internal.hpp>


using namespace std;
using namespace harp;


harp::mpi_psf::mpi_psf ( boost::mpi::communicator const & comm, std::string const & type, boost::property_tree::ptree const & props ) {

  comm_ = comm;

  int rank = comm.rank();

  bool reg_mpi = true;
  plugin_registry & reg = plugin_registry::get( reg_mpi );

  if ( rank == 0 ) {
    // instantiate
    local_.reset ( reg.create_psf ( type, props ) );
  }
    
  // broadcast to all processes

  mpi_comm_buffer_type buf;

  if ( rank == 0 ) {
    mpi_comm_pack ( local_, buf );
  }

  boost::mpi::broadcast ( comm_, buf, 0 );

  if ( rank != 0 ) {
    mpi_comm_unpack ( buf, local_ );
  }

}


mpi_psf * harp::mpi_psf::redistribute ( boost::mpi::communicator const & comm ) {

  mpi_psf * ret = new mpi_psf();

  ret->comm_ = comm;

  if ( comm.rank() == 0 ) {
    ret->local_ = local_;
  }

  mpi_comm_buffer_type buf;

  if ( comm.rank() == 0 ) {
    mpi_comm_pack ( ret->local_, buf );
  }

  boost::mpi::broadcast ( comm, buf, 0 );

  if ( comm.rank() != 0 ) {
    mpi_comm_unpack ( buf, ret->local_ );
  }

  return ret;
}


size_t harp::mpi_psf::n_spec ( ) const {
  return local_->n_spec();
}


size_t harp::mpi_psf::n_lambda ( ) const {
  return local_->n_lambda();
}


size_t harp::mpi_psf::img_rows ( ) const {
  return local_->img_rows();
}


size_t harp::mpi_psf::img_cols ( ) const {
  return local_->img_cols();
}


vector_double harp::mpi_psf::lambda ( ) const {
  return local_->lambda();
}


void harp::mpi_psf::extent ( size_t spec, size_t lambda, size_t & x_offset, size_t & y_offset, size_t & n_x, size_t & n_y ) const {
  local_->extent ( spec, lambda, x_offset, y_offset, n_x, n_y );
  return;
}


void harp::mpi_psf::project_transpose ( std::map < size_t, std::set < size_t > > const & speclambda, mpi_matrix_sparse & AT ) const {

  // check for consistent dimensions

  size_t total = psf::total_bins ( speclambda );

  if ( total != AT.rows() ) {
    HARP_MPI_ABORT( AT.comm().rank(), "Design matrix number of rows does not match length of bin list" );
  }

  // first, select the bins that fall within our local range of rows

  mpi_matrix_sparse_block & block = AT.block();

  std::map < size_t, std::set < size_t > > local_speclambda;

  size_t binoffset = 0;

  for ( std::map < size_t, std::set < size_t > > :: const_iterator itspec = speclambda.begin(); itspec != speclambda.end(); ++itspec ) {

    for ( std::set < size_t > :: const_iterator itlambda = itspec->second.begin(); itlambda != itspec->second.end(); ++itlambda ) {

      if ( ( binoffset >= block.firstrow ) && ( binoffset < ( block.firstrow + block.rows ) ) ) {

        local_speclambda[ itspec->first ].insert ( *itlambda );

      }

      ++binoffset;

    }

  }

  // if we have any rows, get the projection for our local block, and populate the elements of A^T

  if ( ! local_speclambda.empty() ) {

    matrix_double_sparse ( local_AT );

    local_->project_transpose ( local_speclambda, local_AT );

    // get the number of non-zeros in the sparse boost matrix, and use that to reserve space in the
    // local block of the distributed sparse matrix.

    size_t local_nnz = local_AT.nnz();

    block.vals = local_nnz;
    block.row.resize ( local_nnz );
    block.col.resize ( local_nnz );
    block.row_offset.resize ( block.rows );
    block.row_nnz.resize ( block.rows );
    block.data.resize ( local_nnz );

    // iterate over non-zeros of boost sparse matrix and populate local piece of 
    // distributed matrix.  Our local block stores the offsets and nnz of rows that are empty,
    // but the boost row iterator skips empty rows.  So we need fill in the values for these
    // empty rows.

    size_t dense_row = 0;

    size_t elem = 0;

    for ( matrix_double_sparse :: const_iterator1 rowit = local_AT.begin1(); rowit != local_AT.end1(); ++rowit ) {

      size_t row = rowit.index1();

      while ( dense_row < row ) {
        if ( dense_row == 0 ) {
          block.row_offset [ dense_row ] = 0;
          block.row_nnz [ dense_row ] = 0;
        } else {
          block.row_offset [ dense_row ] = block.row_offset [ dense_row - 1 ];
          block.row_nnz [ dense_row ] = 0;
        }
        ++dense_row;
      }

      block.row_offset [ row ] = elem;
      block.row_nnz [ row ] = 0;

      for ( matrix_double_sparse :: const_iterator2 colit = rowit.begin(); colit != rowit.end(); ++colit ) {

        size_t col = colit.index2();

        block.row [ elem ] = row;
        block.col [ elem ] = col;
        block.data [ elem ] = (*colit);

        ++block.row_nnz [ row ];
        ++elem;

      }

      ++dense_row;

    }

  }

  return;
}


void harp::mpi_psf::project_transpose ( mpi_matrix_sparse & AT ) const {

  std::map < size_t, std::set < size_t > > speclambda;

  for ( size_t spec = 0; spec < n_spec(); ++spec ) {
    for ( size_t lambda = 0; lambda < n_lambda(); ++lambda ) {
      speclambda[ spec ].insert ( lambda );
    }
  }

  project_transpose ( speclambda, AT );

  return;
}


std::string harp::mpi_psf::type ( ) const {
  return local_->type();
}


