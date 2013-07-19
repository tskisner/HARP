// @COPYRIGHT@

#include <harp_internal.hpp>


using namespace std;
using namespace harp;


// Produce a design matrix suitable for simultaneous extraction
// and sky subtraction

void harp::sky_design ( matrix_sparse const & AT, std::vector < bool > const & sky, matrix_sparse & skyAT, bool skysub ) {

  int np;
  int myp;

  MPI_Comm_size ( AT.Comm(), &np );
  MPI_Comm_rank ( AT.Comm(), &myp );

  size_t nbins = AT.Height();
  size_t npix = AT.Width();
  size_t nspec = sky.size();

  size_t nsky = 0;
  for ( size_t i = 0; i < nspec; ++i ) {
    if ( sky[i] ) {
      ++nsky;
    }
  }

  size_t nlambda = (size_t)( nbins / nspec );

  if ( nspec * nlambda != nbins ) {
    std::ostringstream o;
    o << "sky vector size (" << nspec << ") does not divide evenly into PSF spectral dimension (" << nbins << ")";
    HARP_THROW( o.str().c_str() );
  }

  size_t sky_nspec = nspec - nsky + 1;
  size_t sky_nbins = sky_nspec * nlambda;

  // build the full mapping of old bins to new bins

  vector < int64_t > old_to_new ( nbins );
  vector < int64_t > old_to_sky ( nbins );

  size_t skystart = (sky_nspec - 1) * nlambda;
  size_t newoff = 0;
  size_t oldoff = 0;

  for ( size_t i = 0; i < nspec; ++i ) {
    for ( size_t j = 0; j < nlambda; ++j ) {
      if ( sky[i] ) {
        old_to_new[ oldoff ] = -1;
      } else {
        old_to_new[ oldoff ] = newoff;
        ++newoff;
      }
      old_to_sky[ oldoff ] = skystart + j;
      ++oldoff;
    }
  }

  /*
  for ( size_t i = 0; i < nbins; ++i ) {
    cerr << "bin map " << i << " --> " << old_to_new[i] << endl;
  }
  */

  // set up new design matrix and get our local range

  skyAT.ResizeTo ( sky_nbins, npix );

  size_t sky_firstrow = skyAT.FirstLocalRow();
  size_t sky_rows = skyAT.LocalHeight();

  size_t firstrow = AT.FirstLocalRow();
  size_t rows = AT.LocalHeight();

  // without some tedious calculations, we do not know how many
  // non-zeros each process will have in the new matrix.  So as a 
  // first guess, we just reserve the number of non-zeros in the
  // old matrix.

  skyAT.Reserve ( AT.NumLocalEntries() );

  double * sky_buf = skyAT.ValueBuffer();

  for ( size_t i = 0; i < skyAT.NumLocalEntries(); ++i ) {
    sky_buf[i] = 0.0;
  }

  // accumulate our local block of the original design matrix.

  size_t offset;
  size_t nnz;
  size_t global_row;
  size_t global_col;
  size_t sky_global_row;

  skyAT.StartAssembly();

  for ( size_t row = 0; row < rows; ++row ) {

    global_row = firstrow + row;
    offset = AT.LocalEntryOffset ( row );
    nnz = AT.NumConnections ( row );

    // fill in object elements
    
    if ( old_to_new[ global_row ] >= 0 ) {

      sky_global_row = old_to_new [ global_row ];

      if ( ( sky_global_row >= sky_firstrow ) && ( sky_global_row < sky_firstrow + sky_rows ) ) {

        for ( size_t col = 0; col < nnz; ++col ) {

          global_col = AT.Col ( offset + col );

          skyAT.Update ( sky_global_row, global_col, AT.Value ( offset + col ) );

        }
      }

    }

    // fill in sky elements

    // FIXME: this currently disables simultaneous splitting of object/sky power
    // in object fibers

    if ( ( old_to_new[ global_row ] < 0 ) || skysub ) {

      sky_global_row = old_to_sky [ global_row ];

      if ( ( sky_global_row >= sky_firstrow ) && ( sky_global_row < sky_firstrow + sky_rows ) ) {

        for ( size_t col = 0; col < nnz; ++col ) {

          global_col = AT.Col ( offset + col );

          skyAT.Update ( sky_global_row, global_col, AT.Value ( offset + col ) );

        }
      }

    }

  }

  skyAT.StopAssembly();

  // In order to build up the new matrix, send our data to the previous rank process
  // and receive from the next rank process.  Repeat this (number of process times)
  // to guarantee that all rows have been filled in.  

  int nshift = np - 1;

  int to_proc;
  if ( myp > 0 ) {
    to_proc = myp - 1;
  } else {
    to_proc = np - 1;
  }

  int from_proc;
  if ( myp == np - 1 ) {
    from_proc = 0;
  } else {
    from_proc = myp + 1;
  }

  char * sendbuf = NULL;
  size_t sendbytes;
  char * recvbuf = NULL;
  size_t recvbytes;
  MPI_Request send_size_request;
  MPI_Request send_request;
  MPI_Status status;

  for ( int shift = 0; shift < nshift; ++shift ) {

    if ( shift == 0 ) {
      // first shift, send our own data

      sparse_block * myblock = new sparse_block ( AT );

      sendbuf = myblock->pack ( sendbytes );

      delete ( myblock );

    } else {
      // pass along the buffer

      sendbytes = recvbytes;

      sendbuf = (char*)malloc ( sendbytes );
      if ( ! sendbuf ) {
        HARP_THROW( "cannot allocate send buffer" );
      }

      memcpy ( (void*)sendbuf, (void*)recvbuf, sendbytes );

      free ( recvbuf );
      
    }

    int temp = MPI_Barrier ( AT.Comm() );

    int send_size_key = (shift * 2 * np) + 2 * myp;
    int send_data_key = (shift * 2 * np) + 2 * myp + 1;
    int recv_size_key = (shift * 2 * np) + 2 * from_proc;
    int recv_data_key = (shift * 2 * np) + 2 * from_proc + 1;

    int ret = MPI_Isend ( (void*)(&sendbytes), 1, MPI_UNSIGNED_LONG, to_proc, send_size_key, AT.Comm(), &send_size_request );
    mpi_check ( AT.Comm(), ret );

    ret = MPI_Isend ( (void*)sendbuf, sendbytes, MPI_CHAR, to_proc, send_data_key, AT.Comm(), &send_request );
    mpi_check ( AT.Comm(), ret );

    // receive block from sender

    ret = MPI_Recv ( (void*)(&recvbytes), 1, MPI_UNSIGNED_LONG, from_proc, recv_size_key, AT.Comm(), &status );
    mpi_check ( AT.Comm(), ret );

    recvbuf = (char*)malloc ( recvbytes );
    if ( ! recvbuf ) {
      HARP_THROW( "cannot allocate receive buffer" );
    }

    ret = MPI_Recv ( (void*)recvbuf, recvbytes, MPI_CHAR, from_proc, recv_data_key, AT.Comm(), &status );
    mpi_check ( AT.Comm(), ret );

    // reconstruct sparse_block

    sparse_block * other_block = new sparse_block ( recvbuf, recvbytes );

    // accumulate block

    skyAT.StartAssembly();

    for ( size_t row = 0; row < other_block->local_rows; ++row ) {

      global_row = other_block->local_firstrow + row;
      offset = other_block->local_row_offset[ row ];
      nnz = other_block->local_row_nnz[ row ];

      // fill in object elements
      
      if ( old_to_new[ global_row ] >= 0 ) {

        sky_global_row = old_to_new [ global_row ];

        if ( ( sky_global_row >= sky_firstrow ) && ( sky_global_row < sky_firstrow + sky_rows ) ) {

          for ( size_t col = 0; col < nnz; ++col ) {

            global_col = other_block->local_col[ offset + col ];

            skyAT.Update ( sky_global_row, global_col, other_block->data[ offset + col ] );

          }

        }

      }
      
      // fill in sky elements

      if ( ( old_to_new[ global_row ] < 0 ) || skysub ) {

        sky_global_row = old_to_sky [ global_row ];

        if ( ( sky_global_row >= sky_firstrow ) && ( sky_global_row < sky_firstrow + sky_rows ) ) {

          for ( size_t col = 0; col < nnz; ++col ) {

            global_col = other_block->local_col[ offset + col ];

            skyAT.Update ( sky_global_row, global_col, other_block->data[ offset + col ] );

          }

        }
      }

    }

    delete other_block;

    skyAT.StopAssembly();

    // free send buffer

    ret = MPI_Wait ( &send_size_request, &status );
    mpi_check ( AT.Comm(), ret );

    ret = MPI_Wait ( &send_request, &status );
    mpi_check ( AT.Comm(), ret );

    free ( sendbuf );

  }

  return;
}


// in case of solving for mean sky, this function subtracts this
// and adds error in quadrature

void harp::sky_subtract ( size_t nobj, matrix_dist const & Rf_orig, matrix_dist const & err_orig, matrix_dist & Rf, matrix_dist & err ) {

  size_t nspec = nobj + 1;
  size_t nlambda = (size_t)( Rf_orig.Height() / nspec );

  if ( nlambda * nspec != Rf_orig.Height() ) {
    std::ostringstream o;
    o << "number of spectra does not divide evenly into vector dimension" ;
    HARP_THROW( o.str().c_str() );
  }

  // get local copy of sky and error bars

  matrix_local local_sky ( nlambda, 1 );
  local_matrix_zero ( local_sky );

  matrix_local local_sky_err ( nlambda, 1 );
  local_matrix_zero ( local_sky_err );

  elem::AxpyInterface < double > globloc;
  globloc.Attach( elem::GLOBAL_TO_LOCAL, Rf_orig );
  globloc.Axpy ( 1.0, local_sky, (nobj * nlambda), 0 );
  globloc.Detach();

  globloc.Attach( elem::GLOBAL_TO_LOCAL, err_orig );
  globloc.Axpy ( 1.0, local_sky_err, (nobj * nlambda), 0 );
  globloc.Detach();

  /*
  for ( size_t d = 0; d < nlambda; ++d ) {
    cerr << d << " " << local_sky.Get(d,0) << " " << local_sky_err.Get(d,0) << endl;
  }
  */

  // go through local spectra and subtract sky, adding errors
  // in quadrature.

  size_t hlocal = Rf_orig.LocalHeight();
  size_t wlocal = Rf_orig.LocalWidth();

  size_t rowoff = Rf_orig.ColShift();
  size_t rowstride = Rf_orig.ColStride();
  size_t row;

  double oval;
  double sval;
  double oerr;
  double serr;
  size_t spec;
  size_t lambda;

  if ( wlocal > 0 ) {
    for ( size_t j = 0; j < hlocal; ++j ) {
      row = rowoff + j * rowstride;
      spec = (size_t)( row / nlambda );
      lambda = row - spec * nlambda;
      //cerr << j << ": row " << row << ", spec " << spec << ", lambda " << lambda << endl;

      oval = Rf_orig.GetLocal ( j, 0 );
      oerr = err_orig.GetLocal ( j, 0 );

      cerr << "    " << oval << " " << oerr << endl;

      if ( row < (nobj * nlambda) ) {
        
        sval = local_sky.Get ( lambda, 0 );
        serr = local_sky_err.Get ( lambda, 0 );

        //cerr << "    new = " << oval << " - " << sval << endl;
        //cerr << "    err = sqrt (" << oerr << "^2 + " << serr << "^2 )" << endl;

        oval -= sval;
        oerr = sqrt ( oerr * oerr + serr * serr );

        Rf.SetLocal ( j, 0, oval );
        err.SetLocal ( j, 0, oerr );

      } else {
        // just copy sky information

        Rf.SetLocal ( j, 0, oval );
        err.SetLocal ( j, 0, oerr );

      }

    }
  }

  return;
}


void harp::sky_design2 ( matrix_sparse const & AT, std::vector < bool > const & sky, matrix_sparse & skyAT ) {

  int np;
  int myp;

  MPI_Comm_size ( AT.Comm(), &np );
  MPI_Comm_rank ( AT.Comm(), &myp );

  size_t nbins = AT.Height();
  size_t npix = AT.Width();
  size_t nspec = sky.size();

  size_t nlambda = (size_t)( nbins / nspec );

  if ( nspec * nlambda != nbins ) {
    std::ostringstream o;
    o << "sky vector size (" << nspec << ") does not divide evenly into PSF spectral dimension (" << nbins << ")";
    HARP_THROW( o.str().c_str() );
  }

  size_t sky_nspec = nspec + 1;
  size_t sky_nbins = sky_nspec * nlambda;
  size_t sky_offset = nspec * nlambda;

  // set up new design matrix and get our local range

  skyAT.ResizeTo ( sky_nbins, npix );

  size_t sky_firstrow = skyAT.FirstLocalRow();
  size_t sky_rows = skyAT.LocalHeight();

  size_t firstrow = AT.FirstLocalRow();
  size_t rows = AT.LocalHeight();

  // without some tedious calculations, we do not know how many
  // non-zeros each process will have in the new matrix.  So as a 
  // first guess, we just reserve the number of non-zeros in the
  // old matrix.

  skyAT.Reserve ( AT.NumLocalEntries() );

  double * sky_buf = skyAT.ValueBuffer();

  for ( size_t i = 0; i < skyAT.NumLocalEntries(); ++i ) {
    sky_buf[i] = 0.0;
  }

  // accumulate our local block of the original design matrix.

  size_t offset;
  size_t nnz;
  size_t global_row;
  size_t global_col;
  size_t sky_global_row;

  skyAT.StartAssembly();

  for ( size_t row = 0; row < rows; ++row ) {

    global_row = firstrow + row;
    offset = AT.LocalEntryOffset ( row );
    nnz = AT.NumConnections ( row );

    // fill in object elements

    if ( ( global_row >= sky_firstrow ) && ( global_row < sky_firstrow + sky_rows ) ) {

      for ( size_t col = 0; col < nnz; ++col ) {

        global_col = AT.Col ( offset + col );

        skyAT.Update ( global_row, global_col, AT.Value ( offset + col ) );

      }
    }

    // fill in sky elements

    sky_global_row = ( global_row % nlambda ) + sky_offset;

    if ( ( sky_global_row >= sky_firstrow ) && ( sky_global_row < sky_firstrow + sky_rows ) ) {

      for ( size_t col = 0; col < nnz; ++col ) {

        global_col = AT.Col ( offset + col );

        skyAT.Update ( sky_global_row, global_col, AT.Value ( offset + col ) );

      }
    }

  }

  skyAT.StopAssembly();

  // In order to build up the new matrix, send our data to the previous rank process
  // and receive from the next rank process.  Repeat this (number of process times)
  // to guarantee that all rows have been filled in.  

  int nshift = np - 1;

  int to_proc;
  if ( myp > 0 ) {
    to_proc = myp - 1;
  } else {
    to_proc = np - 1;
  }

  int from_proc;
  if ( myp == np - 1 ) {
    from_proc = 0;
  } else {
    from_proc = myp + 1;
  }

  char * sendbuf = NULL;
  size_t sendbytes;
  char * recvbuf = NULL;
  size_t recvbytes;
  MPI_Request send_size_request;
  MPI_Request send_request;
  MPI_Status status;

  for ( int shift = 0; shift < nshift; ++shift ) {

    if ( shift == 0 ) {
      // first shift, send our own data

      sparse_block * myblock = new sparse_block ( AT );

      sendbuf = myblock->pack ( sendbytes );

      delete ( myblock );

    } else {
      // pass along the buffer

      sendbytes = recvbytes;

      sendbuf = (char*)malloc ( sendbytes );
      if ( ! sendbuf ) {
        HARP_THROW( "cannot allocate send buffer" );
      }

      memcpy ( (void*)sendbuf, (void*)recvbuf, sendbytes );

      free ( recvbuf );
      
    }

    int temp = MPI_Barrier ( AT.Comm() );

    int send_size_key = (shift * 2 * np) + 2 * myp;
    int send_data_key = (shift * 2 * np) + 2 * myp + 1;
    int recv_size_key = (shift * 2 * np) + 2 * from_proc;
    int recv_data_key = (shift * 2 * np) + 2 * from_proc + 1;

    int ret = MPI_Isend ( (void*)(&sendbytes), 1, MPI_UNSIGNED_LONG, to_proc, send_size_key, AT.Comm(), &send_size_request );
    mpi_check ( AT.Comm(), ret );

    ret = MPI_Isend ( (void*)sendbuf, sendbytes, MPI_CHAR, to_proc, send_data_key, AT.Comm(), &send_request );
    mpi_check ( AT.Comm(), ret );

    // receive block from sender

    ret = MPI_Recv ( (void*)(&recvbytes), 1, MPI_UNSIGNED_LONG, from_proc, recv_size_key, AT.Comm(), &status );
    mpi_check ( AT.Comm(), ret );

    recvbuf = (char*)malloc ( recvbytes );
    if ( ! recvbuf ) {
      HARP_THROW( "cannot allocate receive buffer" );
    }

    ret = MPI_Recv ( (void*)recvbuf, recvbytes, MPI_CHAR, from_proc, recv_data_key, AT.Comm(), &status );
    mpi_check ( AT.Comm(), ret );

    // reconstruct sparse_block

    sparse_block * other_block = new sparse_block ( recvbuf, recvbytes );

    // accumulate block

    skyAT.StartAssembly();

    for ( size_t row = 0; row < other_block->local_rows; ++row ) {

      global_row = other_block->local_firstrow + row;
      offset = other_block->local_row_offset[ row ];
      nnz = other_block->local_row_nnz[ row ];

      // fill in object elements

      if ( ( global_row >= sky_firstrow ) && ( global_row < sky_firstrow + sky_rows ) ) {

        for ( size_t col = 0; col < nnz; ++col ) {

          global_col = other_block->local_col[ offset + col ];

          skyAT.Update ( global_row, global_col, other_block->data[ offset + col ] );

        }

      }
      
      // fill in sky elements

      sky_global_row = ( global_row % nlambda ) + sky_offset;

      if ( ( sky_global_row >= sky_firstrow ) && ( sky_global_row < sky_firstrow + sky_rows ) ) {

        for ( size_t col = 0; col < nnz; ++col ) {

          global_col = other_block->local_col[ offset + col ];

          skyAT.Update ( sky_global_row, global_col, other_block->data[ offset + col ] );

        }

      }

    }

    delete other_block;

    skyAT.StopAssembly();

    // free send buffer

    ret = MPI_Wait ( &send_size_request, &status );
    mpi_check ( AT.Comm(), ret );

    ret = MPI_Wait ( &send_request, &status );
    mpi_check ( AT.Comm(), ret );

    free ( sendbuf );

  }

  return;
}

