#include <test_anasazi.hpp>

#include "AnasaziBlockKrylovSchurSolMgr.hpp"
#include "Epetra_LinearProblem.h"
#include "Epetra_InvOperator.h"

// Include header for AztecOO solver and solver interface for Epetra_Operator
#include "AztecOO.h"
#include "AztecOO_Operator.h"

// Include header for Belos solver and solver interface for Epetra_Operator
#include "BelosEpetraOperator.h"
#include "BelosEpetraAdapter.hpp"

// Ifpack
#include "Ifpack.h"
#include "Ifpack_AdditiveSchwarz.h"

using namespace Anasazi;


int main(int argc, char *argv[]) {

  int np = 1;
  int myp = 1;
  double start_time;
  double stop_time;
  double time_build_A;
  double time_build_rhs;
  double time_build_invC;
  double time_solve;

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

  int info = 0;

  test_params par;

  // We distribute the data by flux bin.  In this case, the design matrix
  // on each process contains all pixels, though many of these rows have
  // no non-zeros.

  
  // construct a row map which includes all pixels

  Epetra_Map pixmap ( par.n_pix, par.n_pix, 0, comm );

  // construct a column map which distributes the flux bins uniformly

  Epetra_Map fluxmap ( par.n_flux, 0, comm );

  // build design matrix

  Teuchos::RCP < Epetra_CrsMatrix > A = build_design_matrix ( par, pixmap, fluxmap, time_build_A );

  // create input signal spectra

  Teuchos::RCP < Epetra_Vector > truth = simulate_signal ( par, A, time_build_rhs );

  // project to get signal image

  Teuchos::RCP < Epetra_Vector > pixel_signal = Teuchos::rcp( new Epetra_Vector ( A->RangeMap() ) );

  A->Multiply ( false, (*truth), (*pixel_signal) );

  Epetra_Map outpixmap ( par.n_pix_lambda, par.n_pix_lambda, 0, comm );
  Teuchos::RCP < Epetra_MultiVector > outpix = Teuchos::rcp( new Epetra_MultiVector ( outpixmap, par.n_pix_trace ) );

  int pix_elem;
  for ( int j = 0; j < par.n_pix_lambda; ++j ) {
    for ( int k = 0; k < par.n_pix_trace; ++k ) {
      pix_elem = ( j * par.n_pix_trace ) + k;
      outpix->ReplaceGlobalValue ( j, k, (*pixel_signal)[ pix_elem ] );
    }
  }

  info = EpetraExt::MultiVectorToMatrixMarketFile ( "signal.out", (*outpix), NULL, NULL, true );
  assert( info==0 );

  // construct noise vector

  Teuchos::RCP < Epetra_Vector > pixel_noise = Teuchos::rcp( new Epetra_Vector ( A->RangeMap() ) );

  // populate noise and pixel noise covariance

  Teuchos::RCP < Epetra_CrsMatrix > invpixcov = noise_covariance ( par, A->RangeMap(), (*pixel_signal), (*pixel_noise), time_build_invC );

  for ( int j = 0; j < par.n_pix_lambda; ++j ) {
    for ( int k = 0; k < par.n_pix_trace; ++k ) {
      pix_elem = ( j * par.n_pix_trace ) + k;
      outpix->ReplaceGlobalValue ( j, k, (*pixel_noise)[ pix_elem ] );
    }
  }

  info = EpetraExt::MultiVectorToMatrixMarketFile ( "noise.out", (*outpix), NULL, NULL, true );
  assert( info==0 );

  Teuchos::RCP < Epetra_Vector > pixel_data = Teuchos::rcp( new Epetra_Vector ( (*pixel_signal) ) );

  for ( int j = 0; j < par.n_pix; ++j ) {
    pixel_data->SumIntoGlobalValues ( 1, &((*pixel_noise)[j]), &j );
  }

  // write out S + N

  for ( int j = 0; j < par.n_pix_lambda; ++j ) {
    for ( int k = 0; k < par.n_pix_trace; ++k ) {
      pix_elem = ( j * par.n_pix_trace ) + k;
      outpix->ReplaceGlobalValue ( j, k, (*pixel_data)[ pix_elem ] );
    }
  }

  info = EpetraExt::MultiVectorToMatrixMarketFile ( "data.out", (*outpix), NULL, NULL, true );
  assert( info==0 );

  // Free up memory

  pixel_noise = Teuchos::null;
  pixel_signal = Teuchos::null;
  truth = Teuchos::null;
  outpix = Teuchos::null;

  // Compute RHS of GLS equation

  #ifdef HAVE_MPI
  MPI_Barrier( MPI_COMM_WORLD );
  start_time = MPI_Wtime();
  #endif

  Teuchos::RCP < Epetra_Vector > rhs = Teuchos::rcp( new Epetra_Vector ( A->DomainMap() ) );

  Teuchos::RCP < Epetra_Vector > pixtemp = Teuchos::rcp( new Epetra_Vector ( A->RangeMap() ) );

  invpixcov->Multiply ( false, (*pixel_data), (*pixtemp) );

  A->Multiply ( true, (*pixtemp), (*rhs) );

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
  Teuchos::RCP < Epetra_CrsMatrix > NA = Teuchos::rcp( new Epetra_CrsMatrix ( Copy, (*A).Graph() ) );

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

  int approx_nz = ( 2 * par.psf_corr ) * ( 2 * (int)( par.psf_corr / ( par.n_pix_gap + 1 ) ) );

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



  // Use Block Krylov iteration with PCG inner iteration using Ifpack preconditioner
  // the mass matrix in this case is the identity.

  // mass matrix ( identity )
  /*
  Teuchos::RCP < Epetra_CrsMatrix > M = Teuchos::rcp( new Epetra_CrsMatrix ( Copy, fluxmap, fluxmap, par.n_flux, false ) );
  double one = 1.0;
  for ( int i = 0; i < par.n_flux; ++i ) {
    M->InsertGlobalValues ( i, 1, &one, &i );
  }
  info = M->FillComplete();
  if ( info != 0 ) {
    std::cerr << "FillComplete M error = " << info << std::endl;
    #ifdef HAVE_MPI
    MPI_Finalize();
    #endif
    return -1;
  }
  */

  // shifted matrix (A - sigma * I)

  Teuchos::RCP < Epetra_CrsMatrix > invCshift = Teuchos::rcp( new Epetra_CrsMatrix ( *invC ) );

  Teuchos::RCP < Epetra_Vector > diag = Teuchos::rcp( new Epetra_Vector ( fluxmap ) );

  invC->ExtractDiagonalCopy ( *diag );

  double shift = 0.0;
  double diagtemp;

  for ( int i = 0; i < par.n_flux; ++i ) {
    diagtemp = (*diag)[i] - shift;
    invCshift->ReplaceGlobalValues ( i, 1, &diagtemp, &i );
  }
  info = invCshift->FillComplete();
  if ( info != 0 ) {
    std::cerr << "FillComplete invCshift error = " << info << std::endl;
    #ifdef HAVE_MPI
    MPI_Finalize();
    #endif
    return -1;
  }



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
    Teuchos::rcp (Factory.Create (PrecType, &*invCshift, OverlapLevel));
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
  int belosblockSize = 3; // block size used by linear solver and eigensolver [ not required to be the same ]
  int maxits = invCshift->NumGlobalRows(); // maximum number of iterations to run
  //
  // Create the Belos::LinearProblem
  //
  Teuchos::RCP<Belos::LinearProblem<double,Epetra_MultiVector,Epetra_Operator> > My_LP = Teuchos::rcp( new Belos::LinearProblem<double,Epetra_MultiVector,Epetra_Operator>() );

  My_LP->setOperator( invCshift );
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
  int nev = 3;
  int blockSize = 3;  
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
  Teuchos::RCP<Epetra_MultiVector> ivec = Teuchos::rcp( new Epetra_MultiVector(invCshift->Map(), blockSize) );
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
    if (myp == 0) {
      cout << "Anasazi::BasicEigenproblem::setProblem() returned with error." << endl;
    }
#ifdef HAVE_MPI
    MPI_Finalize() ;
#endif
    return -1;
  }

  #ifdef HAVE_MPI
  MPI_Barrier( MPI_COMM_WORLD );
  start_time = MPI_Wtime();
  #endif

  // Initialize the Block Arnoldi solver
  Anasazi::BlockKrylovSchurSolMgr<double, MV, OP> MySolverMgr(MyProblem, MyPL);

  // Solve the problem to the specified tolerances or length
  Anasazi::ReturnType returnCode = MySolverMgr.solve();
  if (returnCode != Anasazi::Converged && myp==0) {
    cout << "Anasazi::EigensolverMgr::solve() returned unconverged." << endl;
  }

  // Get the eigenvalues and eigenvectors from the eigenproblem
  Anasazi::Eigensolution<double,MV> sol = MyProblem->getSolution();
  std::vector<Anasazi::Value<double> > evals = sol.Evals;
  Teuchos::RCP<MV> evecs = sol.Evecs;
  int numev = sol.numVecs;

  if (numev > 0) {

    Teuchos::SerialDenseMatrix<int,double> dmatr(numev,numev);
    Epetra_MultiVector tempvec(invCshift->Map(), MVT::GetNumberVecs( *evecs ));
    OPT::Apply( *invCshift, *evecs, tempvec );
    MVT::MvTransMv( 1.0, tempvec, *evecs, dmatr );
    
    if (myp==0) {
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
  time_solve = stop_time - start_time;
  #endif



  // timing dump

  #ifdef HAVE_MPI
  if ( myp == 0 ) {
    std::cout << std::endl;
    std::cout << "Timing information:" << std::endl;
    std::cout << "  Building A matrix:     " << time_build_A << " seconds" << std::endl;
    std::cout << "  Building RHS vector:   " << time_build_rhs << " seconds" << std::endl;
    std::cout << "  Building C^-1 matrix:  " << time_build_invC << " seconds" << std::endl;
    std::cout << "  Compute eigenpairs:  " << time_solve << " seconds" << std::endl;
    std::cout << std::endl;
  }
  #endif

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  return 0;
}

