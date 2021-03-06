#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include <harp_mpi_test.hpp>


using namespace std;
using namespace harp;



void mpi_test_extract_subaccum ( mpi_spec_slice_p slice, mpi_matrix & full_data, mpi_matrix & check_data ) {

  // inter-gang communicator
  boost::mpi::communicator rcomm = slice->rank_comm();

  // intra-gang communicator
  boost::mpi::communicator gcomm = slice->gang_comm();

  // the process rank within this gang
  int np = gcomm.size();
  int myp = gcomm.rank();

  spec_slice_region full_region = slice->full_region();

  std::vector < spec_slice_region > myslice = slice->regions();

  size_t indx = 0;

  mpi_matrix_zero ( check_data );

  size_t hlocal;
  size_t wlocal;

  size_t rowoff;
  size_t rowstride;
  size_t row;

  for ( std::vector < spec_slice_region > :: const_iterator sit = myslice.begin(); sit != myslice.end(); ++sit ) {

    if ( gcomm.rank() == 0 ) {
      ostringstream dbg;
      dbg.str("");
    
      dbg << "Gang " << rcomm.rank() << " processing " << sit->first_spec << " (" << sit->first_good_spec << ") - " << (sit->first_spec + sit->n_spec - 1) << " (" << (sit->first_good_spec + sit->n_good_spec - 1) << ") X " << sit->first_lambda << " (" << sit->first_good_lambda << ") - " << (sit->first_lambda + sit->n_lambda - 1) << " (" << (sit->first_good_lambda + sit->n_good_lambda - 1) << ")" << endl;
      
      cerr << dbg.str();
    }

    size_t sub_nbin = sit->n_spec * sit->n_lambda;

    mpi_matrix sub_data ( sub_nbin, full_data.Width(), full_data.Grid() );
    mpi_matrix_zero ( sub_data );

    mpi_sub_spec ( full_region, (*sit), full_data, false, sub_data );

    hlocal = sub_data.LocalHeight();
    wlocal = sub_data.LocalWidth();

    rowoff = sub_data.ColShift();
    rowstride = sub_data.ColStride();
    size_t row;

    size_t coloff = sub_data.RowShift();
    size_t colstride = sub_data.RowStride();
    size_t col;

    for ( size_t i = 0; i < wlocal; ++i ) {

      for ( size_t j = 0; j < hlocal; ++j ) {

        row = rowoff + j * rowstride;
        col = coloff + i * colstride;
        double specval = sub_data.GetLocal ( j, i );

        size_t spec_offset = (size_t)( row / sit->n_lambda );
        size_t cur_spec = sit->first_spec + spec_offset;

        size_t lambda_offset = row - ( spec_offset * sit->n_lambda );
        size_t cur_lambda = sit->first_lambda + lambda_offset;

        size_t cur_bin = cur_spec * full_region.n_lambda + cur_lambda;

        if ( (size_t)specval != cur_bin ) {
          cerr << "FAIL on sub_spec gang " << rcomm.rank() << ", proc " << gcomm.rank() << ", region " << indx << ", bin " << cur_bin << ", column " << col << ", " << (size_t)specval << " != " << cur_bin << endl;
          exit(1);
        }

      }
    }

    mpi_accum_spec ( (*sit), full_region, sub_data, true, check_data );

    ++indx;

  }

  return;
}




void harp::mpi_test_extract ( string const & datadir ) {

  // global communicator and grid

  boost::mpi::communicator comm;

  int np = comm.size();
  int myp = comm.rank();

  El::Grid grid ( El::mpi::COMM_WORLD );

  // fake rank communicator for use with the global one

  boost::mpi::communicator selfcomm = comm.split ( myp, 0 );

  // intra-gang communicator

  int gangsize = (int)( np / 2 );
  if ( gangsize < 1 ) {
    gangsize = 1;
  }

  int ngang = (int)( np / gangsize );
  int gangtot = ngang * gangsize;
  if ( myp == 0 ) {
    cout << "  Using " << ngang << " gangs of " << gangsize << " processes each" << endl;
  }
  if ( gangtot < np ) {
    if ( myp == 0 ) {
      cout << "  WARNING: " << (np-gangtot) << " processes are idle" << endl;
    }
  }
  int gang = (int)( myp / gangsize );
  int grank = myp % gangsize;
  if ( gang >= ngang ) {
    gang = MPI_UNDEFINED;
    grank = MPI_UNDEFINED;
  }

  boost::mpi::communicator gcomm = comm.split ( gang, grank );

  El::Grid gang_grid ( (MPI_Comm)gcomm );

  // inter-gang communicator

  boost::mpi::communicator rcomm = comm.split ( grank, gang );

  // initialize plugin registry

  bool reg_mpi = true;
  plugin_registry & reg = plugin_registry::get( reg_mpi );



  if ( myp == 0 ) {
    cout << "Testing extraction spectral sub/accum functions..." << endl;
  }

  size_t nspec = 3;
  size_t nlambda = 10;
  size_t chunk_spec = 3;
  size_t overlap_spec = 0;
  size_t chunk_lambda = 2;
  size_t overlap_lambda = 3;

  size_t Rwidth = 2;
  size_t Rband = 2 * Rwidth + 1;

  size_t nbin = nspec * nlambda;

  size_t ncoltest = 10;

  mpi_matrix full_data ( nbin, ncoltest, grid );
  mpi_matrix_zero ( full_data );

  mpi_matrix check_data ( nbin, ncoltest, grid );

  for ( size_t i = 0; i < nbin; ++i ) {
    for ( size_t j = 0; j < ncoltest; ++j ) {
      full_data.Set ( i, j, (double)i );
    }
  }
  
  mpi_spec_slice_p slice ( new mpi_spec_slice ( selfcomm, comm, overlap_spec, overlap_lambda, nspec, nlambda-2*overlap_lambda, chunk_spec, chunk_lambda, overlap_spec, overlap_lambda ) );

  spec_slice_region full_region = slice->full_region();

  mpi_test_extract_subaccum ( slice, full_data, check_data );

  size_t hlocal = check_data.LocalHeight();
  size_t wlocal = check_data.LocalWidth();

  size_t rowoff = check_data.ColShift();
  size_t rowstride = check_data.ColStride();

  size_t coloff = check_data.RowShift();
  size_t colstride = check_data.RowStride();

  for ( size_t i = 0; i < wlocal; ++i ) {
    for ( size_t j = 0; j < hlocal; ++j ) {

      size_t row = rowoff + j * rowstride;
      size_t col = coloff + i * colstride;
      double specval = check_data.GetLocal ( j, i );

      size_t spec_offset = (size_t)( row / full_region.n_lambda );
      size_t cur_spec = full_region.first_spec + spec_offset;

      size_t lambda_offset = row - ( spec_offset * full_region.n_lambda );
      size_t cur_lambda = full_region.first_lambda + lambda_offset;

      size_t cur_bin = cur_spec * full_region.n_lambda + cur_lambda;

      if ( ( cur_spec >= overlap_spec ) && ( cur_spec < (nspec - overlap_spec) ) && ( cur_lambda >= overlap_lambda ) && ( cur_lambda < (nlambda - overlap_lambda) ) ) {

        if ( (size_t)specval != cur_bin ) {
          cerr << "FAIL on accum_spec world " << selfcomm.rank() << ", proc " << comm.rank() << ", bin " << cur_bin << ", column " << col << ", " << (size_t)specval << " != " << cur_bin << endl;
          exit(1);
        }

      }

    }

  }

  if ( myp == 0 ) {
    cout << "  (PASSED)" << endl;
  }

  
  if ( myp == 0 ) {
    cout << "Testing gang-parallel slice and accum..." << endl;
  }

  // setup data

  mpi_matrix_zero ( check_data );

  mpi_matrix gang_full_data ( nbin, ncoltest, gang_grid );
  mpi_matrix gang_check_data ( nbin, ncoltest, gang_grid );

  // world --> gang

  mpi_gang_distribute ( full_data, gang_full_data );

  // define slices

  mpi_spec_slice_p gang_slice ( new mpi_spec_slice ( rcomm, gcomm, overlap_spec, overlap_lambda, nspec, nlambda-2*overlap_lambda, chunk_spec, chunk_lambda, overlap_spec, overlap_lambda ) );

  full_region = slice->full_region();

  // do sub / accum within each gang

  mpi_test_extract_subaccum ( gang_slice, gang_full_data, gang_check_data );

  // accum from all gangs

  mpi_gang_accum ( gang_check_data, check_data );

  hlocal = check_data.LocalHeight();
  wlocal = check_data.LocalWidth();

  rowoff = check_data.ColShift();
  rowstride = check_data.ColStride();

  coloff = check_data.RowShift();
  colstride = check_data.RowStride();

  for ( size_t i = 0; i < wlocal; ++i ) {
    for ( size_t j = 0; j < hlocal; ++j ) {

      size_t row = rowoff + j * rowstride;
      size_t col = coloff + i * colstride;
      double specval = check_data.GetLocal ( j, i );

      size_t spec_offset = (size_t)( row / full_region.n_lambda );
      size_t cur_spec = full_region.first_spec + spec_offset;

      size_t lambda_offset = row - ( spec_offset * full_region.n_lambda );
      size_t cur_lambda = full_region.first_lambda + lambda_offset;

      size_t cur_bin = cur_spec * full_region.n_lambda + cur_lambda;

      if ( ( cur_spec >= overlap_spec ) && ( cur_spec < (nspec - overlap_spec) ) && ( cur_lambda >= overlap_lambda ) && ( cur_lambda < (nlambda - overlap_lambda) ) ) {

        if ( (size_t)specval != cur_bin ) {
          cerr << "FAIL on accum_spec gang " << rcomm.rank() << ", proc " << gcomm.rank() << ", bin " << cur_bin << ", column " << col << ", " << (size_t)specval << " != " << cur_bin << endl;
          exit(1);
        }

      }

    }

  }

  if ( myp == 0 ) {
    cout << "  (PASSED)" << endl;
  }


  if ( myp == 0 ) {
    cout << "Testing MPI inverse covariance construction..." << endl;
  }

  // create spec and read

  double first_lambda = 8000.0;
  double last_lambda = 8200.0;

  boost::property_tree::ptree spec_props;
  spec_props.clear();
  spec_props.put ( "nspec", nspec );
  spec_props.put ( "lambda_n", nlambda );
  spec_props.put ( "lambda_start", first_lambda );
  spec_props.put ( "lambda_stop", last_lambda );
  spec_props.put ( "back", 500.0 );
  spec_props.put ( "atm", 0.0 );
  spec_props.put ( "obj", 0.0 );
  spec_props.put ( "atmspace", 3 );
  spec_props.put ( "skymod", nspec );
  mpi_spec_p testspec ( new mpi_spec ( comm, "sim", spec_props ) );

  mpi_matrix truth ( nbin, 1, grid );
  vector_double lambda;

  testspec->values ( truth );
  testspec->lambda ( lambda );

  // instantiate the PSF

  boost::property_tree::ptree gauss_props;
  gauss_props.put ( "lambda_spec_type", "sim" );
  gauss_props.put_child ( "lambda_spec", spec_props );
  gauss_props.put ( "bundle_size", nspec );
  gauss_props.put ( "nbundle", 1 );
  gauss_props.put ( "fwhm", 6.0 );

  mpi_psf_p gauss_psf ( new mpi_psf ( comm, "gauss_sim", gauss_props ) );

  mpi_psf_p gang_gauss_psf ( gauss_psf->redistribute( gcomm ) );

  // instantiate image and read

  boost::property_tree::ptree img_props;
  img_props.put ( "spec_type", "sim" );
  img_props.put_child ( "spec", spec_props );
  img_props.put ( "psf_type", "gauss_sim" );
  img_props.put_child ( "psf", gauss_props );

  mpi_image_p img ( new mpi_image ( comm, "sim", img_props ) );

  elem_matrix_local img_data;
  elem_matrix_local img_inv;

  img->values ( img_data );
  img->inv_variance ( img_inv );

  // For this test, each gang will use its
  // first region from the slice.

  spec_slice_region test_region = gang_slice->regions()[0];
  size_t test_nbin = test_region.n_spec * test_region.n_lambda;

  map < size_t, set < size_t > > speclambda;

  for ( size_t s = 0; s < test_region.n_spec; ++s ) {
    for ( size_t l = 0; l < test_region.n_lambda; ++l ) {
      speclambda[ s + test_region.first_spec ].insert ( l + test_region.first_lambda );
    }
  }

  // no masking for this test...

  vector_mask mask ( img_inv.Height() );
  for ( size_t i = 0; i < mask.size(); ++i ) {
    mask[i] = 1;
  }

  // get the distributed AT

  mpi_matrix_sparse AT ( gcomm, test_nbin, img_data.Height() );

  gang_gauss_psf->project_transpose ( speclambda, AT );

  // in order to check the sparse distributed matrix, we need to reduce it to the root process

  matrix_double_sparse check_AT;
  check_AT.clear();

  if ( grank != 0 ) {

    gcomm.send ( 0, grank, AT.block() );

  } else {

    ostringstream outAT;
    outAT << datadir << "/test_AT_" << gangsize << "-" << gang << ".out";
    fstream out ( outAT.str().c_str(), ios::out );

    check_AT.resize( AT.rows(), AT.cols(), false );
    check_AT.clear();

    mpi_matrix_sparse_block other_block;

    for ( size_t p = 0; p < gangsize; ++p ) {

      if ( p != 0 ) {
        gcomm.recv ( p, p, other_block );
      }

      mpi_matrix_sparse_block & cur = ( p == 0 ) ? AT.block() : other_block;

      size_t block_firstrow = cur.firstrow;
      size_t block_nrows = cur.rows;
      size_t block_nvals = cur.vals;

      size_t lastrow = 0;
      size_t lastnnz = 0;

      for ( size_t v = 0; v < block_nvals; ++v ) {
        size_t brow = cur.row[ v ];
        if ( brow != lastrow ) {
          if ( v != 0 ) {
            if ( cur.row_nnz[ lastrow ] != lastnnz ) {
              cerr << "sparse matrix block for proc " << p << " has incorrect row_nnz[" << lastrow << "] != " << lastnnz << endl;
              exit(1);
            }
          }

          // first element of this new row, check that the row offset is set to our current element
          if ( cur.row_offset[ brow ] != v ) {
            cerr << "sparse matrix block for proc " << p << " has incorrect row_offset[" << brow << "] != " << v << endl;
            exit(1);
          }
          lastnnz = 0;
          lastrow = brow;
        }

        ++lastnnz;

        size_t bcol = cur.col[ v ];
        double dat = cur.data[ v ];

        check_AT( block_firstrow + brow, bcol ) = dat;
        out << p << " " << v << " " << (block_firstrow + brow) << " " << bcol << " " << dat << endl;
      }

    }

    out.close();

  }

  gcomm.barrier();

  // rank zero in each gang loads serial version of the classes and 
  // creates serial C^-1.  

  matrix_double serial_invC;

  vector_double serial_img_inv;

  vector_double serial_img_data;

  matrix_double_sparse serial_AT;

  if ( grank == 0 ) {

    // instantiate the serial PSF

    psf_p serial_gauss_psf ( reg.create_psf ( "gauss_sim", gauss_props ) );

    serial_AT.resize ( test_nbin, img_data.Height(), false );

    // instantiate image and read

    image_p serial_img ( reg.create_image ( "sim", img_props ) );

    serial_img->inv_variance ( serial_img_inv );
    serial_img->values ( serial_img_data );

    // get design matrix

    serial_gauss_psf->project_transpose ( speclambda, serial_AT );

    // compare to local copy of distributed AT

    for ( size_t i = 0; i < serial_AT.size1(); ++i ) {
      for ( size_t j = 0; j < serial_AT.size2(); ++j ) {

        double serval = serial_AT(i,j);
        double locval = check_AT(i,j);

        if ( fabs( serval ) > std::numeric_limits < double > :: epsilon() ) {
          double rel = fabs ( ( locval - serval ) / serval );
          if ( rel > std::numeric_limits < double > :: epsilon() ) {
            cerr << "FAIL on AT [ " << i << ", " << j << " ], " << locval << " != " << serval << endl;
            exit(1);
          }
        }

      }
    }

    // build matrix

    serial_invC.resize ( test_nbin, test_nbin );

    inverse_covariance ( serial_AT, serial_img_inv, mask, serial_invC );

  }

  // now each gang computes the same in parallel

  mpi_matrix invC ( test_nbin, test_nbin, gang_grid );

  mpi_inverse_covariance ( AT, img_inv, mask, invC );

  // now we get a local copy to the root process of every gang, for comparison purposes

  elem_matrix_local local_invC;

  El::AxpyInterface < double > globloc;
  globloc.Attach( El::GLOBAL_TO_LOCAL, invC );

  if ( grank == 0 ) {
    local_invC.Resize ( test_nbin, test_nbin );
    local_matrix_zero ( local_invC );
    globloc.Axpy ( 1.0, local_invC, 0, 0 );
  }

  globloc.Detach();

  // compare outputs- only the lower triangle!

  if ( grank == 0 ) {

    for ( size_t i = 0; i < test_nbin; ++i ) {
      for ( size_t j = 0; j <= i; ++j ) {

        //cerr << "(" << i << ", " << j << ") " << serial_invC(i,j) << " " << local_invC.Get(i,j) << endl;

        if ( fabs( serial_invC(i,j) ) > std::numeric_limits < double > :: epsilon() ) {
          double rel = fabs ( ( local_invC.Get(i,j) - serial_invC(i,j) ) / serial_invC(i,j) );
          if ( rel > std::numeric_limits < float > :: epsilon() ) {
            cerr << "FAIL on C^-1 [ " << i << ", " << j << " ], " << local_invC.Get(i,j) << " != " << serial_invC(i,j) << endl;
            exit(1);
          }
        }

      }
    }

  }

  if ( myp == 0 ) {
    cout << "  (PASSED)" << endl;
  }


  if ( myp == 0 ) {
    cout << "Testing MPI noise weighted spec construction..." << endl;
  }


  mpi_matrix z_spec ( test_nbin, 1, gang_grid );

  mpi_noise_weighted_spec ( AT, img_inv, mask, img_data, z_spec );

  elem_matrix_local loc_z;

  globloc.Attach( El::GLOBAL_TO_LOCAL, z_spec );
  if ( grank == 0 ) {
    loc_z.Resize ( z_spec.Height(), 1 );
    local_matrix_zero ( loc_z );
    globloc.Axpy ( 1.0, loc_z, 0, 0 );
  }
  globloc.Detach();

  if ( grank == 0 ) {

    vector_double serial_z;

    noise_weighted_spec ( serial_AT, serial_img_inv, mask, serial_img_data, serial_z );

    for ( size_t i = 0; i < test_nbin; ++i ) {
      if ( fabs( serial_z[i] ) > std::numeric_limits < double > :: epsilon() ) {
        double rel = fabs ( ( loc_z.Get(i,0) - serial_z[i] ) / serial_z[i] );
        if ( rel > std::numeric_limits < float > :: epsilon() ) {
          cerr << "FAIL on z [ " << i << " ], " << loc_z.Get(i,0) << " != " << serial_z[i] << endl;
          exit(1);
        }
      }
    }

  }

  if ( myp == 0 ) {
    cout << "  (PASSED)" << endl;
  }


  if ( myp == 0 ) {
    cout << "Testing high-level, chunked extraction..." << endl;
  }
  
  if ( myp == 0 ) {
    string outfile = datadir + "/mpi_extract_image_input.fits.out";
    vector_double udata;
    vector_double uinv;
    elem_to_ublas( img_data, udata );
    elem_to_ublas( img_inv, uinv );
    image_fits::write ( outfile, img->n_rows(), udata, uinv );
  }

  mpi_matrix Rtruth ( nbin, 1, grid );
  mpi_matrix f ( nbin, 1, grid );
  mpi_matrix Rf ( nbin, 1, grid );
  mpi_matrix err ( nbin, 1, grid );
  mpi_matrix Rdiag ( nbin, Rband, grid );

  // do extraction

  bool lambda_mask = true;

  string prefix = "  extract:  ";

  map < string, double > timing;

  mpi_extract_slices ( gang_slice, gauss_psf, img_data, img_inv, truth, Rband, Rdiag, Rf, f, err, Rtruth, timing, lambda_mask, prefix );

  elem_matrix_local loc_truth;
  elem_matrix_local loc_Rtruth;
  elem_matrix_local loc_Rf;
  elem_matrix_local loc_f;
  elem_matrix_local loc_err;

  globloc.Attach( El::GLOBAL_TO_LOCAL, Rtruth );
  if ( myp == 0 ) {
    loc_Rtruth.Resize ( Rtruth.Height(), 1 );
    local_matrix_zero ( loc_Rtruth );
    globloc.Axpy ( 1.0, loc_Rtruth, 0, 0 );
  }
  globloc.Detach();

  globloc.Attach( El::GLOBAL_TO_LOCAL, truth );
  if ( myp == 0 ) {
    loc_truth.Resize ( truth.Height(), 1 );
    local_matrix_zero ( loc_truth );
    globloc.Axpy ( 1.0, loc_truth, 0, 0 );
  }
  globloc.Detach();

  globloc.Attach( El::GLOBAL_TO_LOCAL, Rf );
  if ( myp == 0 ) {
    loc_Rf.Resize ( Rf.Height(), 1 );
    local_matrix_zero ( loc_Rf );
    globloc.Axpy ( 1.0, loc_Rf, 0, 0 );
  }
  globloc.Detach();

  globloc.Attach( El::GLOBAL_TO_LOCAL, f );
  if ( myp == 0 ) {
    loc_f.Resize ( f.Height(), 1 );
    local_matrix_zero ( loc_f );
    globloc.Axpy ( 1.0, loc_f, 0, 0 );
  }
  globloc.Detach();

  globloc.Attach( El::GLOBAL_TO_LOCAL, err );
  if ( myp == 0 ) {
    loc_err.Resize ( err.Height(), 1 );
    local_matrix_zero ( loc_err );
    globloc.Axpy ( 1.0, loc_err, 0, 0 );
  }
  globloc.Detach();

  string outfile;

  if ( myp == 0 ) {

    vector_double ubuf;
    vector_double errbuf;

    elem_to_ublas ( loc_err, errbuf );

    vector_double fake_err ( errbuf );
    fake_err.clear();

    elem_to_ublas ( loc_truth, ubuf );
    outfile = datadir + "/mpi_extract_spec_truth.fits.out";
    spec_fits::write ( outfile, ubuf, fake_err, lambda );

    elem_to_ublas ( loc_Rtruth, ubuf );
    outfile = datadir + "/mpi_extract_spec_Rtruth.fits.out";
    spec_fits::write ( outfile, ubuf, fake_err, lambda );

    elem_to_ublas ( loc_Rf, ubuf );
    outfile = datadir + "/mpi_extract_spec_Rf.fits.out";
    spec_fits::write ( outfile, ubuf, errbuf, lambda );

    elem_to_ublas ( loc_f, ubuf );
    outfile = datadir + "/mpi_extract_spec_f.fits.out";
    spec_fits::write ( outfile, ubuf, fake_err, lambda );

    double chisq_reduced = 0.0;
    double val;
    for ( size_t i = 0; i < loc_Rf.Height(); ++i ) {
      if ( loc_err.Get(i,0) > std::numeric_limits < double > :: epsilon() ) {
        val = ( loc_Rf.Get ( i, 0 ) - loc_Rtruth.Get ( i, 0 ) ) / loc_err.Get ( i, 0 );
        val *= val;
        chisq_reduced += val;
      }
    }

    chisq_reduced /= (double)( loc_Rf.Height() - 1 );

    cout << prefix << "Reduced Chi square = " << chisq_reduced << endl;

    cout << "  (PASSED)" << endl;
  }

  // write out block diagonal resolution

  string outtxt = datadir + "/mpi_extract_spec_res.txt.out";
  ofstream fout;
  
  outfile = datadir + "/mpi_extract_spec_res.fits.out";

  fitsfile * fp;
  int ret;
  int status = 0;
  
  long naxes[3];
  naxes[0] = nlambda;
  naxes[1] = Rband;
  naxes[2] = nspec;

  int fitstype = fits::ftype < double > :: datatype();
  long fpixel[2];

  elem_matrix_local loc_Rdiag;

  size_t write_spec_chunk = 2;
  size_t write_spec_offset = 0;

  // write buffer
  size_t write_spec_nbuf = write_spec_chunk * Rband * nlambda;
  vector_double buffer;

  if ( myp == 0 ) {
    fout.open ( outtxt.c_str(), ios::out );
    fout.precision(3);

    fits::create ( fp, outfile );

    ret = fits_create_img ( fp, fits::ftype< double >::bitpix(), 3, naxes, &status );
    fits::check ( status );

    fits::key_write ( fp, "EXTNAME", string("RESOLUTION"), "" );

    buffer.resize(write_spec_nbuf);
  }

  while ( write_spec_offset < nspec ) {
    if ( write_spec_offset + write_spec_chunk > nspec ) {
      write_spec_chunk = nspec - write_spec_offset;
    }

    if ( myp == 0 ) {
      buffer.clear();
    }

    // we have to copy the data into the buffer one spec
    // at a time, due to memory layout.
    for ( size_t k = 0; k < write_spec_chunk; ++k ) {

      globloc.Attach( El::GLOBAL_TO_LOCAL, Rdiag );
      if ( myp == 0 ) {
        loc_Rdiag.Resize ( nlambda, Rband );
        local_matrix_zero ( loc_Rdiag );
        globloc.Axpy ( 1.0, loc_Rdiag, (write_spec_offset+k)*nlambda, 0 );
      }
      globloc.Detach();

      if ( myp == 0 ) {
        fout << "spec " << k << endl;

        size_t boff = k * Rband * nlambda;

        for ( size_t j = 0; j < Rband; ++j ) {
          for ( size_t i = 0; i < nlambda; ++i ) {
          
            size_t outdiag = (Rband - 1) - j;
            int64_t outlambda = -1;
            
            int64_t loff = (int64_t)j - (int64_t)Rwidth;

            if ( loff >= 0 ) {
              outlambda = (int64_t)i + loff;
              if ( outlambda >= (int64_t)nlambda ) {
                outlambda = -1;
              }
            } else {
              outlambda = (int64_t)i + loff;
            }
            if ( outlambda >= 0 ) {
              double val = loc_Rdiag.Get(i, j);
              if ( fabs(val) < 1.0e-100 ) {
                val = 0.0;
              }
              buffer[boff + outdiag * nlambda + outlambda] = val;
            }
          }
        }

        for ( size_t j = 0; j < Rband; ++j ) {
          for ( size_t i = 0; i < nlambda; ++i ) {        
            fout << buffer[boff+j*nlambda+i] << " ";
          }
          fout << endl;
        }
      }
    }

    if ( myp == 0 ) {
      fpixel[0] = 1;
      fpixel[1] = 1;
      fpixel[2] = write_spec_offset + 1;
      long npix = (long)( write_spec_chunk * Rband * nlambda );

      ret = fits_write_pix ( fp, fitstype, fpixel, npix, &(buffer[0]), &status );
      fits::check ( status );
    }

    write_spec_offset += write_spec_chunk;
  }

  if ( myp == 0 ) {
    fits::close( fp );
    fout.close();
  }


  return;
}



