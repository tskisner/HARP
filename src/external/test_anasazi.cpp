#include <test_anasazi.hpp>

#include <ietl/interface/epetra.h>
#include <ietl/interface/lapack.h>
#include <ietl/vectorspace.h>
#include <ietl/iteration.h>
#include <ietl/lanczos.h>

using namespace Anasazi;


int main(int argc, char *argv[]) {

  double start_time;
  double stop_time;
  double time_build_A;
  double time_build_rhs;
  double time_build_invC;
  double time_lobpcg;
  double time_syev;
  double time_ietl;

#ifdef HAVE_MPI
  // Initialize MPI
  //
  MPI_Init ( &argc, &argv );
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

  test_params par;


  Epetra_Map pixmap ( par.n_pix, 0, comm );

  int my_n_pix = pixmap.NumMyElements();
  std::vector < int > my_pix ( my_n_pix );
  pixmap.MyGlobalElements( &(my_pix[0]) );

  Epetra_Map fluxmap ( par.n_flux, 0, comm );

  int my_n_flux = fluxmap.NumMyElements();
  std::vector < int > my_flux ( my_n_flux );
  fluxmap.MyGlobalElements( &(my_flux[0]) );


  // build transpose of the design matrix

  Teuchos::RCP < Epetra_CrsMatrix > AT = build_design_matrix ( par, fluxmap, pixmap, time_build_A );

  // create input signal spectra

  Teuchos::RCP < Epetra_Vector > truth = simulate_signal ( par, AT, time_build_rhs );

  // project to get signal image

  Teuchos::RCP < Epetra_Vector > pixel_signal = Teuchos::rcp( new Epetra_Vector ( pixmap ) );

  AT->Multiply ( true, (*truth), (*pixel_signal) );


  // construct noise vector

  Teuchos::RCP < Epetra_Vector > pixel_noise = Teuchos::rcp( new Epetra_Vector ( pixmap ) );

  // populate noise and pixel noise covariance

  Teuchos::RCP < Epetra_CrsMatrix > invpixcov = noise_covariance ( par, pixmap, (*pixel_signal), (*pixel_noise), time_build_invC );

  Teuchos::RCP < Epetra_Vector > pixel_data = Teuchos::rcp( new Epetra_Vector ( (*pixel_signal) ) );

  for ( int i = 0; i < my_n_pix; ++i ) {
    pixel_data->SumIntoGlobalValues ( 1, &((*pixel_noise)[ i ]), &(my_pix[i]) );
  }


  // Free up memory

  pixel_noise = Teuchos::null;
  pixel_signal = Teuchos::null;
  truth = Teuchos::null;

  // Compute RHS of GLS equation

  #ifdef HAVE_MPI
  MPI_Barrier( MPI_COMM_WORLD );
  start_time = MPI_Wtime();
  #endif

  Teuchos::RCP < Epetra_Vector > rhs = Teuchos::rcp( new Epetra_Vector ( fluxmap ) );

  Teuchos::RCP < Epetra_Vector > pixtemp = Teuchos::rcp( new Epetra_Vector ( pixmap ) );

  invpixcov->Multiply ( false, (*pixel_data), (*pixtemp) );

  AT->Multiply ( false, (*pixtemp), (*rhs) );

  #ifdef HAVE_MPI
  MPI_Barrier( MPI_COMM_WORLD );
  stop_time = MPI_Wtime();
  time_build_rhs = stop_time - start_time;
  #endif

  // Free up memory

  pixel_data = Teuchos::null;
  pixtemp = Teuchos::null;


  // Compute LHS ( inverse flux bin covariance )

  #ifdef HAVE_MPI
  MPI_Barrier( MPI_COMM_WORLD );
  start_time = MPI_Wtime();
  #endif

  // This is a temporary matrix for the N^-1 x A product...

  Teuchos::RCP < Epetra_CrsMatrix > NA = Teuchos::rcp( new Epetra_CrsMatrix ( Copy, pixmap, 0 ) );

  // First Multiply
  info = EpetraExt::MatrixMatrix::Multiply ( (*invpixcov), false, (*AT), true, (*NA) );
  if ( info != 0 ) {
    std::cerr << "NA MatrixMatrix error = " << info << std::endl;
    #ifdef HAVE_MPI
    MPI_Finalize();
    #endif
    return -1;
  }

  // Create matrix for final inverse flux covariance

  Teuchos::RCP < Epetra_CrsMatrix > invC = Teuchos::rcp( new Epetra_CrsMatrix ( Copy, fluxmap, 0 ) );
  
  // second multiply
  info = EpetraExt::MatrixMatrix::Multiply ( (*AT), false, (*NA), false, (*invC) );
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
  AT = Teuchos::null;
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

  
  const int    nev       = par.n_flux;
  const int    blockSize = par.n_flux;
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
  os << "LOBPCG Solver manager returned " << ( returnCode == Converged ? "converged." : "unconverged." ) << std::endl;
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

  std::cout << endl << endl;

  #ifdef HAVE_MPI
  MPI_Barrier( MPI_COMM_WORLD );
  stop_time = MPI_Wtime();
  time_lobpcg = stop_time - start_time;
  #endif



  // Dense LAPACK version...

  /*

  #ifdef HAVE_MPI
  MPI_Barrier( MPI_COMM_WORLD );
  start_time = MPI_Wtime();
  #endif

  if ( par.myp == 0 ) {
    Teuchos::LAPACK<int, double> lapack;

    Teuchos::SerialDenseMatrix<int, double> s_invC ( par.n_flux, par.n_flux );

    double rowview[ par.n_flux ];
    int rowindx[ par.n_flux ];

    for ( int i = 0; i < par.n_flux; ++i ) {
      int n;
      info = (*invC).ExtractGlobalRowCopy (i, par.n_flux, n, rowview, rowindx);
      for ( int j = 0; j < n; ++j ) {
        s_invC( i, rowindx[j] ) = rowview[ rowindx[j] ];
      }
    }

    int lwork = 3 * par.n_flux - 1;
    double work[ lwork ];
    double s_eigen[ par.n_flux ];

    lapack.SYEV( 'N', 'L', par.n_flux, s_invC.values(), par.n_flux, s_eigen, work, lwork, &info );

    std::cout << std::endl << "SYEV:" << std::endl;
    for ( int i = 0; i < par.n_flux; ++i ) {
      std::cout << "  " << s_eigen[i] << std::endl;
    }
    std::cout << std::endl;

    // simple IETL case

    typedef ietl::vectorspace < ietl::lapack_vector > Vecspace;
    typedef boost::lagged_fibonacci607 Gen;

    Vecspace vec( par.n_flux );
    Gen mygen;

    ietl::lapack_matrix invL ( par.n_flux, par.n_flux );
    for ( int i = 0; i < par.n_flux; ++i ) {
      int n;
      info = (*invC).ExtractGlobalRowCopy (i, par.n_flux, n, rowview, rowindx);
      for ( int j = 0; j < n; ++j ) {
        invL( i, rowindx[j] ) = rowview[ rowindx[j] ];
      }
    }

    ietl::lanczos < ietl::lapack_matrix, Vecspace > lanczos ( invL, vec );


    int max_iter = 2 * par.n_flux;  
    std::cout << "Computation of eigenvalues with fixed size of T-matrix\n\n";
    std::cout << "-----------------------------------\n\n";

    std::vector < double > eigen;
    std::vector < double > err;
    std::vector < int > multiplicity;
    //ietl::lanczos_iteration_nhighest < double > iter ( max_iter, n_flux );
    ietl::fixed_lanczos_iteration < double > iter ( max_iter );  
    try {
      lanczos.calculate_eigenvalues ( iter, mygen );
      eigen = lanczos.eigenvalues();
      err = lanczos.errors();
      multiplicity = lanczos.multiplicities();
    }
    catch (std::runtime_error & e) {
      std::cout << e.what() << std::endl;
    } 

    // Printing eigenvalues with error & multiplicities:  
    std::cout << "#         eigenvalue         error         multiplicity\n";  
    std::cout.precision(10);
    for (size_t i=0;i<eigen.size();++i)
      std::cout << i << "\t" << eigen[i] << "\t" << err[i] << "\t" 
          << multiplicity[i] << "\n";


  }

  #ifdef HAVE_MPI
  MPI_Barrier( MPI_COMM_WORLD );
  stop_time = MPI_Wtime();
  time_syev = stop_time - start_time;
  #endif

  */


  // USE IETL to do the same...

  #ifdef HAVE_MPI
  MPI_Barrier( MPI_COMM_WORLD );
  start_time = MPI_Wtime();
  #endif  

  ietl::ept_matrix invCwrap ( invC );

  typedef ietl::vectorspace < ietl::ept_vector > Vecspace;
  typedef boost::lagged_fibonacci607 Gen;

  Vecspace vec( par.n_flux );
  Gen mygen;

  ietl::lanczos < ietl::ept_matrix, Vecspace > lanczos ( invCwrap, vec );


  int max_iter = 10 * par.n_flux;  
  if ( par.myp == 0 ) {
    std::cout << "Computation of eigenvalues with fixed size of T-matrix\n\n";
    std::cout << "-----------------------------------\n\n";
  }

  std::vector < double > eigen;
  std::vector < double > err;
  std::vector < int > multiplicity;
  //ietl::lanczos_iteration_nhighest < double > iter ( max_iter, n_flux );
  ietl::fixed_lanczos_iteration < double > iter ( max_iter );  
  try {
    lanczos.calculate_eigenvalues ( iter, mygen );
    eigen = lanczos.eigenvalues();
    err = lanczos.errors();
    multiplicity = lanczos.multiplicities();
  }
  catch (std::runtime_error & e) {
    std::cout << e.what() << std::endl;
  } 

  // Printing eigenvalues with error & multiplicities:  
  if ( par.myp == 0 ) {
    std::cout << "#         eigenvalue         error         multiplicity\n";  
    std::cout.precision(10);
    for (size_t i=0;i<eigen.size();++i)
      std::cout << i << "\t" << eigen[i] << "\t" << err[i] << "\t" 
          << multiplicity[i] << "\n";
  }

  #ifdef HAVE_MPI
  MPI_Barrier( MPI_COMM_WORLD );
  stop_time = MPI_Wtime();
  time_ietl = stop_time - start_time;
  #endif



  // timing dump

  #ifdef HAVE_MPI
  if ( par.myp == 0 ) {
    std::cout << std::endl;
    std::cout << "Timing information:" << std::endl;
    std::cout << "  Building A matrix:     " << time_build_A << " seconds" << std::endl;
    std::cout << "  Building RHS vector:   " << time_build_rhs << " seconds" << std::endl;
    std::cout << "  Building C^-1 matrix:  " << time_build_invC << " seconds" << std::endl;
    std::cout << "  LOBPCG Compute eigenpairs:  " << time_lobpcg << " seconds" << std::endl;
    std::cout << "  SYEV Compute eigenpairs:  " << time_syev << " seconds" << std::endl;
    std::cout << "  IETL Compute eigenpairs:  " << time_ietl << " seconds" << std::endl;
    std::cout << std::endl;
  }
  #endif

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  return 0;
}

