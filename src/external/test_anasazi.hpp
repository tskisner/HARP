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
#include "Epetra_Operator.h"
#include "EpetraExt_VectorOut.h"
#include "EpetraExt_MultiVectorOut.h"
#include "EpetraExt_RowMatrixOut.h"
#include "EpetraExt_MatrixMatrix.h"
#include "Teuchos_CommandLineProcessor.hpp"

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

// Pliris
#include "Pliris.h"

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

    int n_bundle;
    int bundle_size;

    int n_trace;
    int n_flux_trace;
    int n_flux;

    int n_pix_margin;
    int n_pix_gap;

    int n_pix_bundle;
    int n_pix_lambda;
    int n_pix_trace;
    int n_pix;

    double peak_amp;
    double peak_obj;

    int peak_space;

    double background;

    double psf_fwhm;

    int psf_corr;

};


class gls_operator : public Epetra_Operator {

  public :

    gls_operator ( Teuchos::RCP < Epetra_CrsMatrix > AT, Teuchos::RCP < Epetra_CrsMatrix > invN ) {
      AT_ = AT;
      invN_ = invN;
      strcpy ( label_, "GLS" );
    }

    ~gls_operator () { }

    gls_operator ( gls_operator const & orig ) {
      AT_ = orig.AT_;
      invN_ = orig.invN_;
      strcpy ( label_, orig.label_ );
    }

    int SetUseTranspose ( bool UseTranspose ) { return -1; }

    int Apply ( const Epetra_MultiVector& X, Epetra_MultiVector& Y ) const {
      Teuchos::RCP < Epetra_MultiVector > temp1 = Teuchos::rcp( new Epetra_MultiVector ( AT_->DomainMap(), X.NumVectors() ) );
      Teuchos::RCP < Epetra_MultiVector > temp2 = Teuchos::rcp( new Epetra_MultiVector ( AT_->DomainMap(), X.NumVectors() ) );
      int info = AT_->Multiply ( true, X, (*temp1) );
      info = invN_->Multiply ( false, (*temp1), (*temp2) );
      info = AT_->Multiply ( false, (*temp2), Y );
      return info;
    }

    int ApplyInverse ( const Epetra_MultiVector& X, Epetra_MultiVector& Y ) const {
      return -1;
    }

    double NormInf() const {
      return 0.0;
    }

    const char * Label() const {
      return label_;
    }

    bool UseTranspose() const {
      return false;
    }

    bool HasNormInf() const {
      return false;
    }

    const Epetra_Comm & Comm() const {
      return AT_->Comm();
    }

    //! Returns the Epetra_Map object associated with the domain of this operator.
    const Epetra_Map & OperatorDomainMap() const {
      return AT_->RangeMap();
    }

    //! Returns the Epetra_Map object associated with the range of this operator.
    const Epetra_Map & OperatorRangeMap() const {
      return AT_->RangeMap();
    }

  private :

    Teuchos::RCP < Epetra_CrsMatrix > AT_;
    Teuchos::RCP < Epetra_CrsMatrix > invN_;
    char label_[20];

};




void block_dist ( int n, int myp, int np, std::vector < int > & my_elems );


void gauss_sample ( std::vector < double > & vals, std::vector < double > & xrel, std::vector < double > & yrel, double amp, double maj, double min, double ang );


Teuchos::RCP < Epetra_CrsMatrix > build_design_matrix ( test_params & par, Epetra_Map & fluxmap, Epetra_Map & pixmap, double & timing );


Teuchos::RCP < Epetra_Vector > simulate_signal ( test_params & par, Teuchos::RCP < Epetra_CrsMatrix > AT, double & timing );


Teuchos::RCP < Epetra_CrsMatrix > noise_covariance ( test_params & par, Epetra_Map const & pixmap, Epetra_Vector & signal, Epetra_Vector & noise, double & timing );


void eigenvals_lobpcg ( test_params & par, Teuchos::RCP < gls_operator > mat, std::vector < double > & vals, double & timing );

void eigenvals_lapack ( test_params & par, Teuchos::RCP < Epetra_CrsMatrix > mat, std::vector < double > & vals, double & timing );

int eigenvals_full ( test_params & par, Teuchos::RCP < Epetra_CrsMatrix > mat, std::vector < double > & vals, double & timing );


void write_frame ( test_params & par, Epetra_Vector & frame, std::string const & filename );


#endif


