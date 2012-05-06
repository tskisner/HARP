#include <vector>
#include <set>
#include <sstream>

#include <boost/random.hpp>

#include "AnasaziConfigDefs.hpp"
#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziSimpleLOBPCGSolMgr.hpp"
#include "AnasaziBasicOutputManager.hpp"
#include "AnasaziEpetraAdapter.hpp"
#include "Epetra_CrsMatrix.h"
#include "EpetraExt_VectorOut.h"
#include "EpetraExt_MultiVectorOut.h"
#include "EpetraExt_RowMatrixOut.h"
#include "EpetraExt_MatrixMatrix.h"
#include "Teuchos_CommandLineProcessor.hpp"

#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#include <mpi.h>
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"

using namespace Anasazi;


// 2*PI
static double const TWOPI = 6.28318530717958647693;


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







int main(int argc, char *argv[]) {

  int np = 1;
  int myp = 1;
  double start_time;
  double stop_time;
  double time_build_A;
  double time_build_rhs;
  double time_build_invC;
  double time_eigen;

#ifdef HAVE_MPI
  // Initialize MPI
  //
  MPI_Init ( &argc, &argv );
  MPI_Comm_size ( MPI_COMM_WORLD, &np );
  MPI_Comm_rank ( MPI_COMM_WORLD, &myp );
#endif

  // Create an Epetra communicator
  //
#ifdef HAVE_MPI
  Epetra_MpiComm comm ( MPI_COMM_WORLD );
#else
  Epetra_SerialComm comm;
#endif

  // Create an Anasazi output manager
  
  BasicOutputManager < double > printer;
  printer.stream ( Errors ) << Anasazi_Version() << endl << endl;

  // Get the sorting string from the command line
  
  std::string which("LM");
  Teuchos::CommandLineProcessor cmdp ( false, true );
  cmdp.setOption ( "sort", &which, "Targetted eigenvalues (SM or LM)." );
  if ( cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL ) {
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return -1;
  }

  int info = 0;


  // testcase:  20 traces of 100 flux bins each.  One bin per pixel in wavelength, 
  // 4 pixels in between traces, 2 pixels on either side of frame.

  int n_trace = 20;
  int n_flux_trace = 100;
  int n_flux = n_trace * n_flux_trace;

  int n_pix_margin = 2;
  int n_pix_gap = 4;

  int n_pix_lambda = n_flux_trace;
  int n_pix_trace = 2 * n_pix_margin + (n_trace - 1) * n_pix_gap + n_trace;
  int n_pix = n_pix_lambda * n_pix_trace;

  std::cerr << "Using " << n_flux << " flux bins (" << n_trace << " x " << n_flux_trace << ")" << std::endl;
  std::cerr << "Using " << n_pix << " pixels (" << n_pix_trace << " x " << n_pix_lambda << ")" << std::endl;

  // amplitude of "spikes"
  double peak_amp = 2000.0;

  // spacing of spikes (10 % of samples)
  int peak_space = (int)( n_flux_trace / 10 );

  // magnitude of continuum
  double background = 20.0;

  // response fwhm
  double psf_fwhm = 1.4;

  // response correlation length in pixels
  int psf_corr = (int)( 6.0 * (double)psf_fwhm );
  std::cerr << "Using corr " << psf_corr << std::endl;


  // We distribute the data by flux bin.  In this case, the design matrix
  // on each process contains all pixels, though many of these rows have
  // no non-zeros.

  #ifdef HAVE_MPI
  MPI_Barrier( MPI_COMM_WORLD );
  start_time = MPI_Wtime();
  #endif

  // construct a row map which includes all pixels

  Epetra_Map pixmap ( n_pix, n_pix, 0, comm );

  // construct a column map which distributes the flux bins uniformly

  Epetra_Map fluxmap ( n_flux, 0, comm );

  // determine locally assigned flux bins

  int my_n_flux = fluxmap.NumMyElements();
  std::vector < int > my_flux ( my_n_flux );
  fluxmap.MyGlobalElements( &(my_flux[0]) );

  std::cerr << "Proc " << myp << " has " << my_n_flux << " flux bins" << std::endl;

  // compute number of non-zeros for my flux bins

  std::vector < int > my_pix_nz ( n_pix );

  int firstX;
  int firstY;
  int lastX;
  int lastY;

  int pix_elem;
  int xcenter;
  int ycenter;

  std::vector < int > :: const_iterator fit;

  for ( fit = my_flux.begin(); fit != my_flux.end(); ++fit ) {
    int trace = (int)( (*fit) / n_flux_trace );
    int trace_bin = (*fit) - trace * n_flux_trace;

    // center pixel
    xcenter = n_pix_margin + ( trace * ( n_pix_gap + 1 ) );
    ycenter = trace_bin;

    if ( xcenter > psf_corr ) {
      firstX = xcenter - psf_corr;
    } else {
      firstX = 0;
    }
    if ( ycenter > psf_corr ) {
      firstY = ycenter - psf_corr;
    } else {
      firstY = 0;
    }
    if ( xcenter + psf_corr < n_pix_trace ) {
      lastX = xcenter + psf_corr;
    } else {
      lastX = n_pix_trace - 1;
    }
    if ( ycenter + psf_corr < n_pix_lambda ) {
      lastY = ycenter + psf_corr;
    } else {
      lastY = n_pix_lambda - 1;
    }

    //std::cerr << "trace " << trace << " bin " << trace_bin << " center(" << xcenter << "," << ycenter << ") " << "X(" << firstX << "," << lastX << ") " << "Y(" << firstY << "," << lastY << ")" << endl;

    for ( int j = firstY; j <= lastY; ++j ) {
      for ( int k = firstX; k <= lastX; ++k ) {
        pix_elem = ( j * n_pix_trace ) + k;
        my_pix_nz[ pix_elem ] += 1;
      }
    }

  }

  int total_nz = 0;
  for ( int i = 0; i < n_pix; ++i ) {
    //std::cerr << "Proc " << myp << ", pix " << i << " has " << my_pix_nz[i] << " local non-zeros" << std::endl;
    total_nz += my_pix_nz[i];
  }
  std::cerr << "Proc " << myp << " has total of non-zeros = " << total_nz << std::endl;


  // Create communication graph for design matrix

  Teuchos::RCP < Epetra_CrsGraph > graph = Teuchos::rcp( new Epetra_CrsGraph ( Copy, pixmap, fluxmap, &(my_pix_nz[0]), true ) );

  for ( fit = my_flux.begin(); fit != my_flux.end(); ++fit ) {
    int trace = (int)( (*fit) / n_flux_trace );
    int trace_bin = (*fit) - trace * n_flux_trace;

    // center pixel
    xcenter = n_pix_margin + ( trace * ( n_pix_gap + 1 ) );
    ycenter = trace_bin;

    if ( xcenter > psf_corr ) {
      firstX = xcenter - psf_corr;
    } else {
      firstX = 0;
    }
    if ( ycenter > psf_corr ) {
      firstY = ycenter - psf_corr;
    } else {
      firstY = 0;
    }
    if ( xcenter + psf_corr < n_pix_trace ) {
      lastX = xcenter + psf_corr;
    } else {
      lastX = n_pix_trace - 1;
    }
    if ( ycenter + psf_corr < n_pix_lambda ) {
      lastY = ycenter + psf_corr;
    } else {
      lastY = n_pix_lambda - 1;
    }

    int iflux_elem = (*fit);

    for ( int j = firstY; j <= lastY; ++j ) {
      for ( int k = firstX; k <= lastX; ++k ) {
        pix_elem = ( j * n_pix_trace ) + k;
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
    int trace = (int)( (*fit) / n_flux_trace );
    int trace_bin = (*fit) - trace * n_flux_trace;

    // center pixel
    xcenter = n_pix_margin + ( trace * ( n_pix_gap + 1 ) );
    ycenter = trace_bin;

    if ( xcenter > psf_corr ) {
      firstX = xcenter - psf_corr;
    } else {
      firstX = 0;
    }
    if ( ycenter > psf_corr ) {
      firstY = ycenter - psf_corr;
    } else {
      firstY = 0;
    }
    if ( xcenter + psf_corr < n_pix_trace ) {
      lastX = xcenter + psf_corr;
    } else {
      lastX = n_pix_trace - 1;
    }
    if ( ycenter + psf_corr < n_pix_lambda ) {
      lastY = ycenter + psf_corr;
    } else {
      lastY = n_pix_lambda - 1;
    }

    std::vector < double > xrel;
    std::vector < double > yrel;
    std::vector < double > vals;
    std::vector < int > pixabs;

    for ( int j = firstY; j <= lastY; ++j ) {
      for ( int k = firstX; k <= lastX; ++k ) {
        pix_elem = ( j * n_pix_trace ) + k;
        pixabs.push_back ( pix_elem );
        xrel.push_back ( static_cast < double > ( k - xcenter ) );
        yrel.push_back ( static_cast < double > ( j - ycenter ) );
      }
    }

    int iflux_elem = (*fit);

    gauss_sample ( vals, xrel, yrel, 1.0, psf_fwhm, psf_fwhm, 0.0 );
    for ( int j = 0; j < (int)pixabs.size(); ++j ) {
      //std::cerr << "replace A(" << pixabs[j] << "," << iflux_elem << ")" << std::endl;
      info = A->ReplaceGlobalValues ( pixabs[j], 1, &(vals[j]), &iflux_elem );
      if ( info != 0 ) {
        std::cerr << "replace error = " << info << std::endl;
        #ifdef HAVE_MPI
        MPI_Finalize();
        #endif
        return -1;
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
    MPI_Finalize();
    #endif
    return -1;
  }

  //A->SetTracebackMode(1); // Shutdown Epetra Warning tracebacks

  #ifdef HAVE_MPI
  MPI_Barrier( MPI_COMM_WORLD );
  stop_time = MPI_Wtime();
  time_build_A = stop_time - start_time;
  #endif

  std::cerr << "Final dimensions A = " << A->NumGlobalRows() << " rows x " << A->NumGlobalCols() << " columns" << std::endl;


  info = EpetraExt::RowMatrixToMatrixMarketFile ( "A.out", (*A), NULL, NULL, true );
  assert( info==0 );

  // create input signal spectra

  Teuchos::RCP < Epetra_Vector > truth = Teuchos::rcp( new Epetra_Vector ( A->DomainMap() ) );

  int half_peak = peak_space >> 1;

  for ( fit = my_flux.begin(); fit != my_flux.end(); ++fit ) {
    int trace = (int)( (*fit) / n_flux_trace );
    int trace_bin = (*fit) - trace * n_flux_trace;
    int iflux_elem = (*fit);
    if ( ( ( trace_bin + half_peak ) % peak_space ) == 0 ) {
      //std::cerr << "peak at trace " << trace << ", bin " << trace_bin << std::endl;
      truth->ReplaceGlobalValues ( 1, &peak_amp, &iflux_elem );
    } else {
      truth->ReplaceGlobalValues ( 1, &background, &iflux_elem );
    }
  }

  info = EpetraExt::VectorToMatrixMarketFile ( "truth.out", (*truth), NULL, NULL, true );
  assert( info==0 );

  // project to get signal image

  Teuchos::RCP < Epetra_Vector > signal = Teuchos::rcp( new Epetra_Vector ( A->RangeMap() ) );

  A->Multiply ( false, (*truth), (*signal) );


  // write out signal

  Epetra_Map outpixmap ( n_pix_lambda, n_pix_lambda, 0, comm );
  Teuchos::RCP < Epetra_MultiVector > outpix = Teuchos::rcp( new Epetra_MultiVector ( outpixmap, n_pix_trace ) );

  for ( int j = 0; j < n_pix_lambda; ++j ) {
    for ( int k = 0; k < n_pix_trace; ++k ) {
      pix_elem = ( j * n_pix_trace ) + k;
      outpix->ReplaceGlobalValue ( j, k, (*signal)[ pix_elem ] );
    }
  }


  info = EpetraExt::MultiVectorToMatrixMarketFile ( "signal.out", (*outpix), NULL, NULL, true );
  assert( info==0 );

  // create (diagonal) inverse pixel noise covariance

  Teuchos::RCP < Epetra_CrsMatrix > invpixcov = Teuchos::rcp( new Epetra_CrsMatrix ( Copy, A->RangeMap(), 1 ) );

  // add white noise and populate inverse pixel covariance

  Teuchos::RCP < Epetra_Vector > noise = Teuchos::rcp( new Epetra_Vector ( A->RangeMap() ) );

  Teuchos::RCP < Epetra_Vector > measured = Teuchos::rcp( new Epetra_Vector ( A->RangeMap() ) );

  typedef boost::ecuyer1988 base_generator_type;
  base_generator_type generator(42u);
  
  double rms;
  double invval;

  for ( int i = 0; i < n_pix; ++i ) {
    rms = sqrt( 16.0 + (*signal)[i] );

    invval = 1.0 / ( rms * rms );

    invpixcov->InsertGlobalValues ( i, 1, &invval, &i );
    
    boost::normal_distribution < double > dist ( 0.0, rms );
    
    boost::variate_generator < base_generator_type&, boost::normal_distribution < double > > gauss ( generator, dist );
    
    (*noise)[i] = gauss();

    (*measured)[i] = (*noise)[i] + (*signal)[i];
  }

  info = invpixcov->FillComplete();
  if ( info != 0 ) {
    std::cerr << "FillComplete invpixcov error = " << info << std::endl;
    #ifdef HAVE_MPI
    MPI_Finalize();
    #endif
    return -1;
  }

  // write out noise

  for ( int j = 0; j < n_pix_lambda; ++j ) {
    for ( int k = 0; k < n_pix_trace; ++k ) {
      pix_elem = ( j * n_pix_trace ) + k;
      outpix->ReplaceGlobalValue ( j, k, (*noise)[ pix_elem ] );
    }
  }

  info = EpetraExt::MultiVectorToMatrixMarketFile ( "noise.out", (*outpix), NULL, NULL, true );
  assert( info==0 );

  // write out S + N

  for ( int j = 0; j < n_pix_lambda; ++j ) {
    for ( int k = 0; k < n_pix_trace; ++k ) {
      pix_elem = ( j * n_pix_trace ) + k;
      outpix->ReplaceGlobalValue ( j, k, (*measured)[ pix_elem ] );
    }
  }


  info = EpetraExt::MultiVectorToMatrixMarketFile ( "measured.out", (*outpix), NULL, NULL, true );
  assert( info==0 );

  // Free up memory

  noise = Teuchos::null;
  signal = Teuchos::null;
  truth = Teuchos::null;
  outpix = Teuchos::null;

  // Compute RHS of GLS equation

  #ifdef HAVE_MPI
  MPI_Barrier( MPI_COMM_WORLD );
  start_time = MPI_Wtime();
  #endif

  Teuchos::RCP < Epetra_Vector > rhs = Teuchos::rcp( new Epetra_Vector ( A->DomainMap() ) );

  Teuchos::RCP < Epetra_Vector > pixtemp = Teuchos::rcp( new Epetra_Vector ( A->RangeMap() ) );

  invpixcov->Multiply ( false, (*measured), (*pixtemp) );

  A->Multiply ( true, (*pixtemp), (*rhs) );

  #ifdef HAVE_MPI
  MPI_Barrier( MPI_COMM_WORLD );
  stop_time = MPI_Wtime();
  time_build_rhs = stop_time - start_time;
  #endif

  // Free up memory

  measured = Teuchos::null;
  pixtemp = Teuchos::null;

  // Compute LHS ( inverse flux bin covariance )

  #ifdef HAVE_MPI
  MPI_Barrier( MPI_COMM_WORLD );
  start_time = MPI_Wtime();
  #endif

  // This is a temporary matrix for the N^-1 x A product...
  Teuchos::RCP < Epetra_CrsMatrix > NA = Teuchos::rcp( new Epetra_CrsMatrix ( Copy, (*graph) ) );

  // First Multiply
  info = EpetraExt::MatrixMatrix::Multiply ( (*invpixcov), false, (*A), false, (*NA), true );
  if ( info != 0 ) {
    std::cerr << "NA MatrixMatrix error = " << info << std::endl;
    #ifdef HAVE_MPI
    MPI_Finalize();
    #endif
    return -1;
  }

  // Create matrix for final inverse flux covariance

  // Create communication graph with an estimate of the number of non-zeros, based on
  // the pixel-space correlation distance and spacing of the flux bins.

  int approx_nz = ( 2 * psf_corr ) * ( 2 * (int)( psf_corr / ( n_pix_gap + 1 ) ) );

  //Teuchos::RCP < Epetra_CrsGraph > Cgraph = Teuchos::rcp( new Epetra_CrsGraph ( Copy, fluxmap, fluxmap, 0, false ) );

  // create the matrix
  Teuchos::RCP < Epetra_CrsMatrix > invC = Teuchos::rcp( new Epetra_CrsMatrix ( Copy, fluxmap, fluxmap, approx_nz, false ) );

  // second multiply
  info = EpetraExt::MatrixMatrix::Multiply ( (*A), true, (*NA), false, (*invC), true );
  if ( info != 0 ) {
    std::cerr << "A^TNA MatrixMatrix error = " << info << std::endl;
    #ifdef HAVE_MPI
    MPI_Finalize();
    #endif
    return -1;
  }

  #ifdef HAVE_MPI
  MPI_Barrier( MPI_COMM_WORLD );
  stop_time = MPI_Wtime();
  time_build_invC = stop_time - start_time;
  #endif

  // Free up memory

  NA = Teuchos::null;
  A = Teuchos::null;
  invpixcov = Teuchos::null;

  // Now we finally have the RHS (Z) vector, and the inverse spectral covariance...

  info = EpetraExt::RowMatrixToMatrixMarketFile ( "invC.out", (*invC), NULL, NULL, true );
  assert( info==0 );

  
  //************************************
  // Call the LOBPCG solver manager
  //***********************************
  //
  //  Variables used for the LOBPCG Method
  //

  #ifdef HAVE_MPI
  MPI_Barrier( MPI_COMM_WORLD );
  start_time = MPI_Wtime();
  #endif


  const int    nev       = 100;
  const int    blockSize = 10;
  const int    maxIters  = 1000;
  const double tol       = 1.0e-8;

  typedef Epetra_MultiVector MV;
  typedef Epetra_Operator OP;
  typedef MultiVecTraits < double, Epetra_MultiVector > MVT;

  // Create an Epetra_MultiVector for an initial vector to start the solver.
  // Note:  This needs to have the same number of columns as the blocksize.

  Teuchos::RCP < Epetra_MultiVector > ivec = Teuchos::rcp ( new Epetra_MultiVector ( fluxmap, blockSize ) );
  ivec->Random();

  // Create the eigenproblem.
  
  Teuchos::RCP<BasicEigenproblem < double, MV, OP > > MyProblem =
    Teuchos::rcp( new BasicEigenproblem < double, MV, OP > ( invC, ivec ) );

  // Inform the eigenproblem that the operator A is symmetric
  
  MyProblem->setHermitian(true);

  // Set the number of eigenvalues requested
  
  MyProblem->setNEV( nev );

  // Inform the eigenproblem that you are finishing passing it information

  bool boolret = MyProblem->setProblem();
  if (boolret != true) {
    printer.print(Errors,"Anasazi::BasicEigenproblem::setProblem() returned an error.\n");
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return -1;
  }

  // Create parameter list to pass into the solver manager

  Teuchos::ParameterList MyPL;
  MyPL.set( "Which", which );
  MyPL.set( "Block Size", blockSize );
  MyPL.set( "Maximum Iterations", maxIters );
  MyPL.set( "Convergence Tolerance", tol );
  
  // Create the solver manager
  SimpleLOBPCGSolMgr < double, MV, OP > MySolverMan ( MyProblem, MyPL );

  // Solve the problem

  ReturnType returnCode = MySolverMan.solve();

  // Get the eigenvalues and eigenvectors from the eigenproblem

  Eigensolution < double, MV > sol = MyProblem->getSolution();
  std::vector < Value < double > > evals = sol.Evals;
  Teuchos::RCP < MV > evecs = sol.Evecs;

  #ifdef HAVE_MPI
  MPI_Barrier( MPI_COMM_WORLD );
  stop_time = MPI_Wtime();
  time_eigen = stop_time - start_time;
  #endif

  // Compute residuals.
  
  std::vector < double > normR ( sol.numVecs );
  if ( sol.numVecs > 0 ) {
    Teuchos::SerialDenseMatrix < int,double > T ( sol.numVecs, sol.numVecs );
    Epetra_MultiVector tempAevec ( fluxmap, sol.numVecs );
    T.putScalar ( 0.0 ); 
    for ( int i = 0; i < sol.numVecs; i++ ) {
      T ( i, i ) = evals[i].realpart;
    }
    invC->Apply ( *evecs, tempAevec );
    MVT::MvTimesMatAddMv ( -1.0, *evecs, T, 1.0, tempAevec );
    MVT::MvNorm ( tempAevec, normR );
  }

  // Print the results

  std::ostringstream os;
  os.setf ( std::ios_base::right, std::ios_base::adjustfield );
  os << "Solver manager returned " << ( returnCode == Converged ? "converged." : "unconverged." ) << std::endl;
  os << std::endl;
  os << "------------------------------------------------------" << endl;
  os << std::setw(16) << "Eigenvalue"
    << std::setw(18) << "Direct Residual"
    << std::endl;
  os << "------------------------------------------------------" << endl;
  for ( int i = 0; i < sol.numVecs; i++ ) {
    os << std::setw(16) << evals[i].realpart
      << std::setw(18) << ( normR[i] / evals[i].realpart )
      << std::endl;
  }
  os << "------------------------------------------------------" << endl;
  printer.print ( Errors, os.str() );
  

  #ifdef HAVE_MPI
  if ( myp == 0 ) {
    std::cout << std::endl;
    std::cout << "Timing information:" << std::endl;
    std::cout << "  Building A matrix:     " << time_build_A << " seconds" << std::endl;
    std::cout << "  Building RHS vector:   " << time_build_rhs << " seconds" << std::endl;
    std::cout << "  Building C^-1 matrix:  " << time_build_invC << " seconds" << std::endl;
    std::cout << "  Compute " << nev << "/" << n_flux << " eigenpairs:  " << time_eigen << " seconds" << std::endl;
    std::cout << std::endl;
  }
  #endif

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  return 0;
}

