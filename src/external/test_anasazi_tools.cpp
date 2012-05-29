
#include <test_anasazi.hpp>

#include <cmath>

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

  n_bundle = 5;
  bundle_size = 20;

  n_trace = n_bundle * bundle_size;
  n_flux_trace = 4000; //20;
  n_flux = n_trace * n_flux_trace;

  n_pix_margin = 5;
  n_pix_gap = 7;
  n_pix_lambda = n_flux_trace;

  n_pix_bundle = 2 * n_pix_margin + (bundle_size - 1) * n_pix_gap + bundle_size;
  n_pix_trace = n_pix_bundle * n_bundle;
  n_pix = n_pix_lambda * n_pix_trace;

  std::cout << "Using " << n_flux << " flux bins (" << n_trace << " x " << n_flux_trace << ")" << std::endl;
  std::cout << "Using " << n_pix << " pixels (" << n_pix_trace << " x " << n_pix_lambda << ")" << std::endl;

  // magnitude of continuum
  background = 20.0;

  // amplitude of atmosphere "spikes"
  peak_amp = 100.0 * background;

  // spacing of spikes
  peak_space = (int)( n_flux_trace / 50 );

  // object size
  peak_obj = 3.0 * background;

  // response fwhm
  psf_fwhm = 1.5;

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
  #pragma omp parallel for default(none) private(i, xt, yt) shared(nvals, xrel, yrel, buf, cang, sang, invmaj, invmin, vals) schedule(static)
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

  int last = -1;
  cerr << endl;

  for ( fit = my_flux.begin(); fit != my_flux.end(); ++fit ) {
    int trace = (int)( (*fit) / par.n_flux_trace );
    int trace_bin = (*fit) - trace * par.n_flux_trace;

    int bundle = (int)( trace / par.bundle_size );
    int bundle_trace = trace - bundle * par.bundle_size;

    if ( trace != last ) {
      fprintf ( stderr, "\rTrace %d        ", trace );
      last = trace;
    }

    // center pixel
    xcenter = (bundle * par.n_pix_bundle) + par.n_pix_margin + ( bundle_trace * ( par.n_pix_gap + 1 ) );
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

  cerr << endl;

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

  //info = EpetraExt::RowMatrixToMatrixMarketFile ( "A.out", (*AT), NULL, NULL, true );
  //assert( info==0 );

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
  int half_obj = half_peak / 3;
  int obj_space = par.peak_space / 3;

  std::vector < int > :: const_iterator fit;

  for ( fit = my_flux.begin(); fit != my_flux.end(); ++fit ) {
    int trace = (int)( (*fit) / par.n_flux_trace );
    int trace_bin = (*fit) - trace * par.n_flux_trace;
    int iflux_elem = (*fit);

    double val = par.background * sin ( 3.0 * (double)trace_bin / ( par.peak_space * 2.0 * M_PI ) ) * sin ( 7.0 * (double)trace_bin / ( par.peak_space * 2.0 * M_PI ) ) * sin ( 11.0 * (double)trace_bin / ( par.peak_space * 2.0 * M_PI ) );

    val += 2.0 * par.background;

    if ( ( ( trace_bin + half_peak ) % par.peak_space ) == 0 ) {
      //std::cerr << "peak at trace " << trace << ", bin " << trace_bin << std::endl;
      val += par.peak_amp;
    }

    if ( ( ( trace_bin + (trace % 25) + half_obj ) % obj_space ) == 0 ) {
      val += par.peak_obj;
    }
    truth->ReplaceGlobalValues ( 1, &(val), &iflux_elem );
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


void eigenvals_lobpcg ( test_params & par, Teuchos::RCP < gls_operator > mat, std::vector < double > & vals, double & timing ) {
//void eigenvals ( test_params & par, Teuchos::RCP < Epetra_CrsMatrix > mat, std::vector < double > & vals, double & timing ) {

  double start_time;
  double stop_time;

  #ifdef HAVE_MPI
  MPI_Barrier( MPI_COMM_WORLD );
  start_time = MPI_Wtime();
  #endif

  Anasazi::BasicOutputManager < double > printer;

  int nev = 50;
  int blockSize = 10;
  int maxIters = 1000;
  double tol = 1.0e-12;

  typedef Epetra_MultiVector MV;
  typedef Epetra_Operator OP;
  typedef MultiVecTraits < double, Epetra_MultiVector > MVT;

  vals.resize( nev );

  // Create an Epetra_MultiVector for an initial vector to start the solver.
  // Note:  This needs to have the same number of columns as the blocksize.

  Teuchos::RCP < Epetra_MultiVector > ivec = Teuchos::rcp ( new Epetra_MultiVector ( mat->OperatorDomainMap(), blockSize ) );
  ivec->Random();

  // Create the eigenproblem.
  
  Teuchos::RCP< Anasazi::BasicEigenproblem < double, MV, OP > > MyProblem = Teuchos::rcp( new Anasazi::BasicEigenproblem < double, MV, OP > ( mat, ivec ) );

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
    return;
  }

  // Create parameter list to pass into the solver manager

  Teuchos::ParameterList MyPL;
  MyPL.set( "Which", "LM" );
  MyPL.set( "Block Size", blockSize );
  MyPL.set( "Maximum Iterations", maxIters );
  MyPL.set( "Convergence Tolerance", tol );
  
  // Create the solver manager
  Anasazi::SimpleLOBPCGSolMgr < double, MV, OP > MySolverMan ( MyProblem, MyPL );

  // Solve the problem

  Anasazi::ReturnType returnCode = MySolverMan.solve();

  // Get the eigenvalues from the eigenproblem

  Anasazi::Eigensolution < double, MV > sol = MyProblem->getSolution();
  std::vector < Anasazi::Value < double > > evals = sol.Evals;

  // Print the results

  std::ostringstream os;
  os.setf ( std::ios_base::right, std::ios_base::adjustfield );
  os << "LOBPCG Solver manager returned " << ( returnCode == Anasazi::Converged ? "converged." : "unconverged." ) << std::endl;
  os << std::endl;
  os << "------------------------------------------------------" << endl;
  os << std::setw(16) << "Eigenvalue" << std::endl;
  os << "------------------------------------------------------" << endl;
  for ( int i = 0; i < sol.numVecs; i++ ) {
    os << std::setw(16) << evals[i].realpart << std::endl;
  }
  os << "------------------------------------------------------" << endl;
  printer.print ( Anasazi::Errors, os.str() );

  std::cout << endl << endl;

  #ifdef HAVE_MPI
  MPI_Barrier( MPI_COMM_WORLD );
  stop_time = MPI_Wtime();
  timing = stop_time - start_time;
  #endif


  return;
}


void eigenvals_lapack ( test_params & par, Teuchos::RCP < Epetra_CrsMatrix > mat, std::vector < double > & vals, double & timing ) {

  double start_time;
  double stop_time;

  #ifdef HAVE_MPI
  MPI_Barrier( MPI_COMM_WORLD );
  start_time = MPI_Wtime();
  #endif

  int info = 0;

  Teuchos::LAPACK<int, double> lapack;

  Teuchos::SerialDenseMatrix<int, double> s_invC ( par.n_flux, par.n_flux );

  double rowview[ par.n_flux ];
  int rowindx[ par.n_flux ];

  for ( int i = 0; i < par.n_flux; ++i ) {
    int n;
    info = (*mat).ExtractGlobalRowCopy (i, par.n_flux, n, rowview, rowindx);
    for ( int j = 0; j < n; ++j ) {
      s_invC( i, rowindx[j] ) = rowview[ rowindx[j] ];
    }
  }

  int lwork = 3 * par.n_flux - 1;
  double work[ lwork ];

  vals.resize ( par.n_flux );

  lapack.SYEV( 'N', 'L', par.n_flux, s_invC.values(), par.n_flux, &(vals[0]), work, lwork, &info );

  std::cout << std::endl << "SYEV:" << std::endl;
  for ( int i = 0; i < par.n_flux; ++i ) {
    std::cout << "  " << vals[i] << std::endl;
  }
  std::cout << std::endl;


  #ifdef HAVE_MPI
  MPI_Barrier( MPI_COMM_WORLD );
  stop_time = MPI_Wtime();
  timing = stop_time - start_time;
  #endif

  return;
}


/*
int eigenvals_full ( test_params & par, Teuchos::RCP < Epetra_CrsMatrix > mat, std::vector < double > & vals, double & timing ) {
//void eigenvals ( test_params & par, Teuchos::RCP < Epetra_CrsMatrix > mat, std::vector < double > & vals, double & timing ) {

  double start_time;
  double stop_time;

  #ifdef HAVE_MPI
  MPI_Barrier( MPI_COMM_WORLD );
  start_time = MPI_Wtime();
  #endif

  // =============================================================== //
  // B E G I N N I N G   O F   I F P A C K   C O N S T R U C T I O N //
  // =============================================================== //

  Teuchos::ParameterList List;

  // Allocate an IFPACK factory.  The object contains no data, only
  // the Create() method for creating preconditioners.
  Ifpack Factory;

  // Create the preconditioner.  For the list of PrecType values that
  // Create() accepts, please check the IFPACK documentation.
  string PrecType = "ILU"; // incomplete LU
  int OverlapLevel = 1; // must be >= 0. If Comm.NumProc() == 1,
                        // it is ignored.

  Teuchos::RCP<Ifpack_Preconditioner> Prec = 
    Teuchos::rcp (Factory.Create (PrecType, &*mat, OverlapLevel));
  TEST_FOR_EXCEPTION(Prec == Teuchos::null, std::runtime_error,
                     "IFPACK failed to create a preconditioner of type \"" 
                     << PrecType << "\" with overlap level " 
                     << OverlapLevel << ".");

  // Specify parameters for ILU.  ILU is local to each MPI process.
  List.set("fact: drop tolerance", 1e-10);
  List.set("fact: level-of-fill", 1);

  // IFPACK uses overlapping Schwarz domain decomposition over all
  // participating processes to combine the results of ILU on each
  // process.  IFPACK's Schwarz method can use any of the following
  // combine modes to combine overlapping results:
  //
  // "Add", "Zero", "Insert", "InsertAdd", "Average", "AbsMax"
  //
  // The Epetra_CombineMode.h header file defines their meaning.
  List.set("schwarz: combine mode", "Add");
  // Set the parameters.
  IFPACK_CHK_ERR(Prec->SetParameters(List));

  // Initialize the preconditioner. At this point the matrix must have
  // been FillComplete()'d, but actual values are ignored.
  IFPACK_CHK_ERR(Prec->Initialize());

  // Build the preconditioner, by looking at the values of the matrix.
  IFPACK_CHK_ERR(Prec->Compute());

  // Create the Belos preconditioned operator from the Ifpack preconditioner.
  // NOTE:  This is necessary because Belos expects an operator to apply the
  //        preconditioner with Apply() NOT ApplyInverse().
  Teuchos::RCP<Belos::EpetraPrecOp> belosPrec = Teuchos::rcp (new Belos::EpetraPrecOp (Prec));

  // =================================================== //
  // E N D   O F   I F P A C K   C O N S T R U C T I O N //
  // =================================================== //


  // ******************************************************
  // Set up Belos Block GMRES operator for inner iteration
  // ******************************************************
  //
  int belosblockSize = 2; // block size used by linear solver and eigensolver [ not required to be the same ]
  int maxits = mat->NumGlobalRows(); // maximum number of iterations to run
  //
  // Create the Belos::LinearProblem
  //
  Teuchos::RCP<Belos::LinearProblem<double,Epetra_MultiVector,Epetra_Operator> > My_LP = Teuchos::rcp( new Belos::LinearProblem<double,Epetra_MultiVector,Epetra_Operator>() );

  My_LP->setOperator( mat );
  My_LP->setRightPrec( belosPrec );

  //
  // Create the ParameterList for the Belos Operator
  // 
  Teuchos::RCP<Teuchos::ParameterList> My_List = Teuchos::rcp( new Teuchos::ParameterList() );
  My_List->set( "Solver", "BlockCG" );
  My_List->set( "Maximum Iterations", maxits );
  My_List->set( "Block Size", belosblockSize );
  My_List->set( "Convergence Tolerance", 1e-13 );
  //
  // Create the Belos::EpetraOperator
  //
  Teuchos::RCP<Belos::EpetraOperator> BelosOp = Teuchos::rcp( new Belos::EpetraOperator( My_LP, My_List ));


  // ************************************
  // Start the block Arnoldi iteration
  // ************************************
  //
  //  Variables used for the Block Arnoldi Method
  //
  double tol = 1.0e-12;
  int nev = 10;
  int blockSize = 5;  
  int numBlocks = 3*nev / blockSize;
  int maxRestarts = 20;
  //int step = 5;
  std::string which = "LM";
  int verbosity = Anasazi::Errors + Anasazi::Warnings + Anasazi::FinalSummary;
  //
  // Create parameter list to pass into solver
  //
  Teuchos::ParameterList MyPL;
  MyPL.set( "Verbosity", verbosity );
  MyPL.set( "Which", which );
  //MyPL.set( "Extra NEV Blocks", Overlap );
  MyPL.set( "Block Size", blockSize );
  MyPL.set( "Num Blocks", numBlocks );
  MyPL.set( "Maximum Restarts", maxRestarts );
  MyPL.set( "Convergence Tolerance", tol );
  //MyPL.set( "Step Size", step );
  
  typedef Epetra_MultiVector MV;
  typedef Epetra_Operator OP;
  typedef Anasazi::MultiVecTraits<double, MV> MVT;
  typedef Anasazi::OperatorTraits<double, MV, OP> OPT;
  
  // Create an Epetra_MultiVector for an initial vector to start the solver.
  // Note:  This needs to have the same number of columns as the blocksize.
  Teuchos::RCP<Epetra_MultiVector> ivec = Teuchos::rcp( new Epetra_MultiVector( mat->DomainMap(), blockSize) );
  MVT::MvRandom( *ivec );
  
  // Call the ctor that calls the petra ctor for a matrix
  //Teuchos::RCP<Anasazi::EpetraGenOp> Aop = Teuchos::rcp( new Anasazi::EpetraGenOp(precOperator, M) );
  
  //Teuchos::RCP<Anasazi::BasicEigenproblem<double,MV,OP> > MyProblem = Teuchos::rcp( new Anasazi::BasicEigenproblem<double,MV,OP>(Aop, M, ivec) );

  Teuchos::RCP<Anasazi::BasicEigenproblem<double,MV,OP> > MyProblem = Teuchos::rcp( new Anasazi::BasicEigenproblem<double,MV,OP>(BelosOp, ivec) );

  //Teuchos::RCP<Anasazi::BasicEigenproblem<double,MV,OP> > MyProblem = Teuchos::rcp( new Anasazi::BasicEigenproblem<double,MV,OP>(invC, ivec) );
  
  // Inform the eigenproblem that the matrix pencil (K,M) is symmetric
  MyProblem->setHermitian(true);
  
  // Set the number of eigenvalues requested 
  MyProblem->setNEV( nev );

  // Inform the eigenproblem that you are finished passing it information
  bool boolret = MyProblem->setProblem();
  if (boolret != true) {
    if ( par.myp == 0 ) {
      cout << "Anasazi::BasicEigenproblem::setProblem() returned with error." << endl;
    }
#ifdef HAVE_MPI
    MPI_Finalize() ;
#endif
    return -1;
  }

  // Initialize the Block Arnoldi solver
  Anasazi::BlockKrylovSchurSolMgr<double, MV, OP> MySolverMgr(MyProblem, MyPL);

  // Solve the problem to the specified tolerances or length
  Anasazi::ReturnType returnCode = MySolverMgr.solve();
  if (returnCode != Anasazi::Converged && par.myp==0) {
    cout << "Anasazi::EigensolverMgr::solve() returned unconverged." << endl;
  }

  // Get the eigenvalues and eigenvectors from the eigenproblem
  Anasazi::Eigensolution<double,MV> sol = MyProblem->getSolution();
  std::vector<Anasazi::Value<double> > evals = sol.Evals;
  vals.resize ( evals.size() );
  Teuchos::RCP<MV> evecs = sol.Evecs;
  int numev = sol.numVecs;

  if (numev > 0) {

    Teuchos::SerialDenseMatrix<int,double> dmatr(numev,numev);
    Epetra_MultiVector tempvec(mat->Map(), MVT::GetNumberVecs( *evecs ));
    OPT::Apply( *mat, *evecs, tempvec );
    MVT::MvTransMv( 1.0, tempvec, *evecs, dmatr );
    
    if ( par.myp == 0 ) {
      double compeval = 0.0;
      cout.setf(std::ios_base::right, std::ios_base::adjustfield);
      cout<<"Actual Eigenvalues (obtained by Rayleigh quotient) : "<<endl;
      cout<<"------------------------------------------------------"<<endl;
      cout<<std::setw(16)<<"Real Part"
        <<std::setw(16)<<"Rayleigh Error"<<endl;
      cout<<"------------------------------------------------------"<<endl;
      for ( int i=0; i<numev; i++) {
        compeval = dmatr(i,i);
        cout<<std::setw(16)<<compeval
          <<std::setw(16)<<Teuchos::ScalarTraits<double>::magnitude(compeval-1.0/evals[i].realpart)
          <<endl;
      }
      cout<<"------------------------------------------------------"<<endl;
    }
    
  }

  #ifdef HAVE_MPI
  MPI_Barrier( MPI_COMM_WORLD );
  stop_time = MPI_Wtime();
  timing = stop_time - start_time;
  #endif


  return 0;
}

*/




