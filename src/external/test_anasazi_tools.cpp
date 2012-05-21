
#include <test_anasazi.hpp>

using namespace Anasazi;


test_params::test_params ( ) {

  np = 1;
  myp = 0;

  #ifdef HAVE_MPI
  MPI_Comm_size ( MPI_COMM_WORLD, &np );
  MPI_Comm_rank ( MPI_COMM_WORLD, &myp );
  #endif

  // testcase:  20 traces of 100 flux bins each.  One bin per pixel in wavelength, 
  // 4 pixels in between traces, 2 pixels on either side of frame.

  n_trace = 2;
  n_flux_trace = 10;
  n_flux = n_trace * n_flux_trace;

  n_pix_margin = 2;
  n_pix_gap = 4;

  n_pix_lambda = n_flux_trace;
  n_pix_trace = 2 * n_pix_margin + (n_trace - 1) * n_pix_gap + n_trace;
  n_pix = n_pix_lambda * n_pix_trace;

  std::cout << "Using " << n_flux << " flux bins (" << n_trace << " x " << n_flux_trace << ")" << std::endl;
  std::cout << "Using " << n_pix << " pixels (" << n_pix_trace << " x " << n_pix_lambda << ")" << std::endl;

  // amplitude of "spikes"
  peak_amp = 2000.0;

  // spacing of spikes (10 % of samples)
  peak_space = (int)( n_flux_trace / 10 );

  // magnitude of continuum
  background = 20.0;

  // response fwhm
  psf_fwhm = 1.2;

  // response correlation length in pixels
  psf_corr = (int)( 5.0 * (double)psf_fwhm );
  std::cout << "Using corr " << psf_corr << std::endl;

}


void block_dist ( int n, int myp, int np, std::vector < int > & my_elems ) {

  int myn = (int) ( n / np );
  
  int leftover = n % np;
  
  int offset;

  if ( myp < leftover ) {
    ++myn;
    offset = myn * myp;
  } else {
    offset = ( (myn + 1) * leftover ) + ( myn * (myp - leftover) );
  }

  my_elems.resize ( myn );

  for ( int i = 0; i < myn; ++i ) {
    my_elems[i] = offset + i;
  }

  return;
}


void gauss_sample ( std::vector < double > & vals, std::vector < double > & xrel, std::vector < double > & yrel, double amp, double maj, double min, double ang ) {

  vals.resize( xrel.size() );
  
  amp /= maj * min * TWOPI;
  
  size_t nvals = xrel.size();
  
  double cang = cos ( ang );
  double sang = sin ( ang );
  
  double invmaj = 1.0 / maj;
  double invmin = 1.0 / min;
  
  std::vector < double > buf ( nvals );
  
  double xt, yt;
  
  size_t i;
  
  #ifdef _OPENMP
  #pragma omp parallel for default(none) private(i, xt, yt) shared(nvals, xrel, yrel, buf, cang, sang, invmaj, invmin) schedule(static)
  #endif
  for ( i = 0; i < nvals; ++i ) {
    xt = xrel[i] * cang + yrel[i] * sang;
    yt = - xrel[i] * sang + yrel[i] * cang;
    buf[i] = - 0.5 * ( xt * xt * invmaj * invmaj + yt * yt * invmin * invmin );
    // vectorize this...
    vals[i] = exp( buf[i] );
  }
  
  #ifdef _OPENMP
  #pragma omp parallel for default(none) private(i) shared(nvals, vals, amp) schedule(static)
  #endif
  for ( i = 0; i < nvals; ++i ) {
    vals[i] *= amp;
  }
  
  return;
}


Teuchos::RCP < Epetra_CrsMatrix > build_design_matrix ( test_params & par, Epetra_Map & fluxmap, Epetra_Map & pixmap, double & timing ) {

  double start_time;
  double stop_time;

  #ifdef HAVE_MPI
  MPI_Barrier( MPI_COMM_WORLD );
  start_time = MPI_Wtime();
  #endif

  int my_n_flux = fluxmap.NumMyElements();
  std::vector < int > my_flux ( my_n_flux );
  fluxmap.MyGlobalElements( &(my_flux[0]) );

  int firstX;
  int firstY;
  int lastX;
  int lastY;

  int info = 0;

  int pix_elem;
  int xcenter;
  int ycenter;

  std::vector < int > :: const_iterator fit;

  // Create transpose design matrix

  Teuchos::RCP < Epetra_CrsMatrix > AT = Teuchos::rcp( new Epetra_CrsMatrix ( Copy, fluxmap, 0 ) );

  // Populate design matrix

  for ( fit = my_flux.begin(); fit != my_flux.end(); ++fit ) {
    int trace = (int)( (*fit) / par.n_flux_trace );
    int trace_bin = (*fit) - trace * par.n_flux_trace;

    // center pixel
    xcenter = par.n_pix_margin + ( trace * ( par.n_pix_gap + 1 ) );
    ycenter = trace_bin;

    if ( xcenter > par.psf_corr ) {
      firstX = xcenter - par.psf_corr;
    } else {
      firstX = 0;
    }
    if ( ycenter > par.psf_corr ) {
      firstY = ycenter - par.psf_corr;
    } else {
      firstY = 0;
    }
    if ( xcenter + par.psf_corr < par.n_pix_trace ) {
      lastX = xcenter + par.psf_corr;
    } else {
      lastX = par.n_pix_trace - 1;
    }
    if ( ycenter + par.psf_corr < par.n_pix_lambda ) {
      lastY = ycenter + par.psf_corr;
    } else {
      lastY = par.n_pix_lambda - 1;
    }

    std::vector < double > xrel;
    std::vector < double > yrel;
    std::vector < double > vals;
    std::vector < int > pixabs;

    for ( int j = firstY; j <= lastY; ++j ) {
      for ( int k = firstX; k <= lastX; ++k ) {
        pix_elem = ( j * par.n_pix_trace ) + k;
        pixabs.push_back ( pix_elem );
        xrel.push_back ( static_cast < double > ( k - xcenter ) );
        yrel.push_back ( static_cast < double > ( j - ycenter ) );
      }
    }

    gauss_sample ( vals, xrel, yrel, 1.0, par.psf_fwhm, par.psf_fwhm, 0.0 );

    info = AT->InsertGlobalValues ( (*fit), (int)pixabs.size(), &(vals[0]), &(pixabs[0]) );

  }

  // Finish up
  
  info = AT->FillComplete( pixmap, fluxmap );

  if ( info != 0 ) {
    std::cerr << "FillComplete error = " << info << std::endl;
    #ifdef HAVE_MPI
    MPI_Abort(MPI_COMM_WORLD, 1);
    #endif
    return Teuchos::null;
  }

  //A->SetTracebackMode(1); // Shutdown Epetra Warning tracebacks

  #ifdef HAVE_MPI
  MPI_Barrier( MPI_COMM_WORLD );
  stop_time = MPI_Wtime();
  timing += stop_time - start_time;
  #endif

  std::cout << "Final dimensions A^T = " << AT->NumGlobalRows() << " rows x " << AT->NumGlobalCols() << " columns" << std::endl;

  info = EpetraExt::RowMatrixToMatrixMarketFile ( "A.out", (*AT), NULL, NULL, true );
  assert( info==0 );

  return AT;
}


Teuchos::RCP < Epetra_Vector > simulate_signal ( test_params & par, Teuchos::RCP < Epetra_CrsMatrix > AT, double & timing ) {

  double start_time;
  double stop_time;

  int info;

  #ifdef HAVE_MPI
  MPI_Barrier( MPI_COMM_WORLD );
  start_time = MPI_Wtime();
  #endif

  int my_n_flux = (AT->RangeMap()).NumMyElements();
  std::vector < int > my_flux ( my_n_flux );
  (AT->RangeMap()).MyGlobalElements( &(my_flux[0]) );

  Teuchos::RCP < Epetra_Vector > truth = Teuchos::rcp( new Epetra_Vector ( AT->RangeMap() ) );

  int half_peak = par.peak_space >> 1;

  std::vector < int > :: const_iterator fit;

  for ( fit = my_flux.begin(); fit != my_flux.end(); ++fit ) {
    int trace = (int)( (*fit) / par.n_flux_trace );
    int trace_bin = (*fit) - trace * par.n_flux_trace;
    int iflux_elem = (*fit);
    if ( ( ( trace_bin + half_peak ) % par.peak_space ) == 0 ) {
      //std::cerr << "peak at trace " << trace << ", bin " << trace_bin << std::endl;
      truth->ReplaceGlobalValues ( 1, &(par.peak_amp), &iflux_elem );
    } else {
      truth->ReplaceGlobalValues ( 1, &(par.background), &iflux_elem );
    }
  }

  #ifdef HAVE_MPI
  MPI_Barrier( MPI_COMM_WORLD );
  stop_time = MPI_Wtime();
  timing += stop_time - start_time;
  #endif

  info = EpetraExt::VectorToMatrixMarketFile ( "truth.out", (*truth), NULL, NULL, true );
  assert( info==0 );

  return truth;
}


Teuchos::RCP < Epetra_CrsMatrix > noise_covariance ( test_params & par, Epetra_Map const & pixmap, Epetra_Vector & signal, Epetra_Vector & noise, double & timing ) {

  double start_time;
  double stop_time;

  int info;

  #ifdef HAVE_MPI
  MPI_Barrier( MPI_COMM_WORLD );
  start_time = MPI_Wtime();
  #endif

  int my_n_pix = pixmap.NumMyElements();
  std::vector < int > my_pix ( my_n_pix );
  pixmap.MyGlobalElements( &(my_pix[0]) );

  // create (diagonal) inverse pixel noise covariance

  Teuchos::RCP < Epetra_CrsMatrix > invpixcov = Teuchos::rcp( new Epetra_CrsMatrix ( Copy, pixmap, 1 ) );

  // add white noise and populate inverse pixel covariance

  typedef boost::ecuyer1988 base_generator_type;
  base_generator_type generator( (unsigned) par.myp );
  
  double rms;
  double invval;

  for ( int i = 0; i < my_n_pix; ++i ) {
    rms = sqrt( 16.0 + signal[i] );

    invval = 1.0 / ( rms * rms );

    //std::cerr << "proc " << par.myp << " inserting invpixcov diagonal " << i << std::endl;

    invpixcov->InsertGlobalValues ( my_pix[i], 1, &invval, &(my_pix[i]) );
    
    boost::normal_distribution < double > dist ( 0.0, rms );
    
    boost::variate_generator < base_generator_type&, boost::normal_distribution < double > > gauss ( generator, dist );

    noise[i] = gauss();
  }


  #ifdef HAVE_MPI
  MPI_Barrier( MPI_COMM_WORLD );
  #endif

  info = invpixcov->FillComplete( pixmap, pixmap );
  
  if ( info != 0 ) {
    std::cerr << "FillComplete invpixcov error = " << info << std::endl;
    #ifdef HAVE_MPI
    MPI_Abort(MPI_COMM_WORLD, 1);
    #endif
    return Teuchos::null;
  }

  /*
  for ( int p = 0; p < par.np; ++p ) {
    if ( par.myp == p ) {
      for ( int r = 0; r < par.n_pix; ++r ) {
        int num;
        std::vector < double > vals ( par.n_pix );
        std::vector < int > indx ( par.n_pix );
        invpixcov->ExtractGlobalRowCopy (r, par.n_pix, num, &(vals[0]), &(indx[0]) ); 
        std::cerr << "Proc " << par.myp << " (" << r << ") = ";
        for ( int c = 0; c < num; ++c ) {
           std::cerr << "(" << indx[c] << ":" << vals[c] << ") ";
        }
        std::cerr << std::endl << std::endl;
      }
    }
    #ifdef HAVE_MPI
    MPI_Barrier( MPI_COMM_WORLD );
    #endif
  }
  */

  #ifdef HAVE_MPI
  MPI_Barrier( MPI_COMM_WORLD );
  stop_time = MPI_Wtime();
  timing += stop_time - start_time;
  #endif

  return invpixcov;
}


void write_frame ( test_params & par, Epetra_Vector & frame, std::string const & filename ) {

  Epetra_Map outpixmap ( par.n_pix_lambda, 0, frame.Comm() );

  Teuchos::RCP < Epetra_MultiVector > outpix = Teuchos::rcp( new Epetra_MultiVector ( outpixmap, par.n_pix_trace ) );

  outpix->PutScalar ( 0.0 );

  int my_n_pix = frame.Map().NumMyElements();
  std::vector < int > my_pix ( my_n_pix );
  frame.Map().MyGlobalElements( &(my_pix[0]) );

  for ( int i = 0; i < my_n_pix; ++i ) {
    int bin = (int)( my_pix[i] / par.n_pix_trace );
    int trace = my_pix[i] - ( bin * par.n_pix_trace );
    outpix->ReplaceGlobalValue ( bin, trace, frame[ i ] );
  }

  int info = EpetraExt::MultiVectorToMatrixMarketFile ( filename.c_str(), (*outpix), NULL, NULL, true );
  assert( info==0 );

  return;
}





