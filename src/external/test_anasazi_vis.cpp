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

  write_frame ( par, (*pixel_signal), "signal.out" );
  write_frame ( par, (*pixel_noise), "noise.out" );
  write_frame ( par, (*pixel_data), "data.out" );

  #ifdef HAVE_MPI
    MPI_Finalize();
  #endif
  return 0;
}

