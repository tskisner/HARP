
#include <test_anasazi.hpp>

using namespace Anasazi;


test_params::test_params ( ) {

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

  std::cerr << "Using " << n_flux << " flux bins (" << n_trace << " x " << n_flux_trace << ")" << std::endl;
  std::cerr << "Using " << n_pix << " pixels (" << n_pix_trace << " x " << n_pix_lambda << ")" << std::endl;

  // amplitude of "spikes"
  peak_amp = 2000.0;

  // spacing of spikes (10 % of samples)
  peak_space = (int)( n_flux_trace / 10 );

  // magnitude of continuum
  background = 20.0;

  // response fwhm
  psf_fwhm = 1.4;

  // response correlation length in pixels
  psf_corr = (int)( 6.0 * (double)psf_fwhm );
  std::cerr << "Using corr " << psf_corr << std::endl;

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


Teuchos::RCP < Epetra_CrsMatrix > build_design_matrix ( test_params & par, Epetra_Map & pixmap, Epetra_Map & fluxmap, double & timing ) {

  double start_time;
  double stop_time;

  #ifdef HAVE_MPI
  MPI_Barrier( MPI_COMM_WORLD );
  start_time = MPI_Wtime();
  #endif

  int my_n_flux = fluxmap.NumMyElements();
  std::vector < int > my_flux ( my_n_flux );
  fluxmap.MyGlobalElements( &(my_flux[0]) );

  std::vector < int > my_pix_nz ( par.n_pix );

  int firstX;
  int firstY;
  int lastX;
  int lastY;

  int info = 0;

  int pix_elem;
  int xcenter;
  int ycenter;

  // compute my non-zeros

  std::vector < int > :: const_iterator fit;

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

    //std::cerr << "trace " << trace << " bin " << trace_bin << " center(" << xcenter << "," << ycenter << ") " << "X(" << firstX << "," << lastX << ") " << "Y(" << firstY << "," << lastY << ")" << endl;

    for ( int j = firstY; j <= lastY; ++j ) {
      for ( int k = firstX; k <= lastX; ++k ) {
        pix_elem = ( j * par.n_pix_trace ) + k;
        my_pix_nz[ pix_elem ] += 1;
      }
    }

  }

  int total_nz = 0;
  for ( int i = 0; i < par.n_pix; ++i ) {
    //std::cerr << "Proc " << myp << ", pix " << i << " has " << my_pix_nz[i] << " local non-zeros" << std::endl;
    total_nz += my_pix_nz[i];
  }

  // Create communication graph for design matrix

  Teuchos::RCP < Epetra_CrsGraph > graph = Teuchos::rcp( new Epetra_CrsGraph ( Copy, pixmap, fluxmap, &(my_pix_nz[0]), true ) );

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

    int iflux_elem = (*fit);

    for ( int j = firstY; j <= lastY; ++j ) {
      for ( int k = firstX; k <= lastX; ++k ) {
        pix_elem = ( j * par.n_pix_trace ) + k;
        graph->InsertGlobalIndices ( pix_elem, 1, &iflux_elem );
      }
    }

  }

  graph->FillComplete ( fluxmap, pixmap );
  //graph.PrintGraphData ( std::cerr );

  // Create design matrix

  Teuchos::RCP < Epetra_CrsMatrix > A = Teuchos::rcp( new Epetra_CrsMatrix ( Copy, (*graph) ) );

  std::cerr << "Start dimensions A = " << A->NumGlobalRows() << " rows x " << A->NumGlobalCols() << " columns" << std::endl;


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

    int iflux_elem = (*fit);

    gauss_sample ( vals, xrel, yrel, 1.0, par.psf_fwhm, par.psf_fwhm, 0.0 );
    for ( int j = 0; j < (int)pixabs.size(); ++j ) {
      //std::cerr << "replace A(" << pixabs[j] << "," << iflux_elem << ")" << std::endl;
      info = A->ReplaceGlobalValues ( pixabs[j], 1, &(vals[j]), &iflux_elem );
      if ( info != 0 ) {
        std::cerr << "replace error = " << info << std::endl;
        #ifdef HAVE_MPI
        MPI_Abort(MPI_COMM_WORLD, 1);
        #endif
        return Teuchos::null;
      }
    }

  }

  // Finish up
  
  //Epetra_CrsGraph const & Agraph = A->Graph();
  //Agraph.PrintGraphData( std::cerr );

  info = A->FillComplete( pixmap, fluxmap, false );
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

  std::cerr << "Final dimensions A = " << A->NumGlobalRows() << " rows x " << A->NumGlobalCols() << " columns" << std::endl;

  info = EpetraExt::RowMatrixToMatrixMarketFile ( "A.out", (*A), NULL, NULL, true );
  assert( info==0 );

  return A;
}


Teuchos::RCP < Epetra_Vector > simulate_signal ( test_params & par, Teuchos::RCP < Epetra_CrsMatrix > A, double & timing ) {

  double start_time;
  double stop_time;

  int info;

  #ifdef HAVE_MPI
  MPI_Barrier( MPI_COMM_WORLD );
  start_time = MPI_Wtime();
  #endif

  int my_n_flux = (A->DomainMap()).NumMyElements();
  std::vector < int > my_flux ( my_n_flux );
  (A->DomainMap()).MyGlobalElements( &(my_flux[0]) );

  Teuchos::RCP < Epetra_Vector > truth = Teuchos::rcp( new Epetra_Vector ( A->DomainMap() ) );

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


Teuchos::RCP < Epetra_CrsMatrix > noise_covariance ( test_params & par, Epetra_Map const & map, Epetra_Vector & signal, Epetra_Vector & noise, double & timing ) {

  double start_time;
  double stop_time;

  #ifdef HAVE_MPI
  MPI_Barrier( MPI_COMM_WORLD );
  start_time = MPI_Wtime();
  #endif

  // create (diagonal) inverse pixel noise covariance

  Teuchos::RCP < Epetra_CrsMatrix > invpixcov = Teuchos::rcp( new Epetra_CrsMatrix ( Copy, map, 1 ) );

  // add white noise and populate inverse pixel covariance

  typedef boost::ecuyer1988 base_generator_type;
  base_generator_type generator(42u);
  
  double rms;
  double invval;

  int info;

  for ( int i = 0; i < par.n_pix; ++i ) {
    rms = sqrt( 16.0 + signal[i] );

    invval = 1.0 / ( rms * rms );

    invpixcov->InsertGlobalValues ( i, 1, &invval, &i );
    
    boost::normal_distribution < double > dist ( 0.0, rms );
    
    boost::variate_generator < base_generator_type&, boost::normal_distribution < double > > gauss ( generator, dist );
    
    noise[i] = gauss();
  }

  info = invpixcov->FillComplete();
  if ( info != 0 ) {
    std::cerr << "FillComplete invpixcov error = " << info << std::endl;
    #ifdef HAVE_MPI
    MPI_Abort(MPI_COMM_WORLD, 1);
    #endif
    return Teuchos::null;
  }

  #ifdef HAVE_MPI
  MPI_Barrier( MPI_COMM_WORLD );
  stop_time = MPI_Wtime();
  timing += stop_time - start_time;
  #endif

  return invpixcov;
}



