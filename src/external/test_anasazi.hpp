#ifndef TEST_ANASAZI_HPP
#define TEST_ANASAZI_HPP

#include <vector>
#include <set>
#include <sstream>

#include <boost/random.hpp>

#include "Teuchos_LAPACK.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"

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


// 2*PI
static double const TWOPI = 6.28318530717958647693;


class test_params {

  public :

    test_params ( );
    ~test_params ( ) { }

    int np;
    int myp;

    int n_trace;
    int n_flux_trace;
    int n_flux;

    int n_pix_margin;
    int n_pix_gap;

    int n_pix_lambda;
    int n_pix_trace;
    int n_pix;

    double peak_amp;

    int peak_space;

    double background;

    double psf_fwhm;

    int psf_corr;

};


void block_dist ( int n, int myp, int np, std::vector < int > & my_elems );


void gauss_sample ( std::vector < double > & vals, std::vector < double > & xrel, std::vector < double > & yrel, double amp, double maj, double min, double ang );


Teuchos::RCP < Epetra_CrsMatrix > build_design_matrix ( test_params & par, Epetra_Map & fluxmap, Epetra_Map & pixmap, double & timing );


Teuchos::RCP < Epetra_Vector > simulate_signal ( test_params & par, Teuchos::RCP < Epetra_CrsMatrix > AT, double & timing );


Teuchos::RCP < Epetra_CrsMatrix > noise_covariance ( test_params & par, Epetra_Map const & pixmap, Epetra_Vector & signal, Epetra_Vector & noise, double & timing );


void write_frame ( test_params & par, Epetra_Vector & frame, std::string const & filename );


#endif


