// @COPYRIGHT@

#include <harp_math_internal.hpp>


using namespace std;
using namespace harp;


void harp::check_column_major ( boost::numeric::ublas::column_major_tag ) {
  return;
}


void harp::check_column_major ( boost::numeric::ublas::row_major_tag ) {
  HARP_THROW( "Matrices must be column-major, in order to maintain compatibility with LAPACK" );
  return;
}


void harp::eigen_decompose ( matrix_double const & invcov, vector_double & D, matrix_double & W ) {

  D.resize ( invcov.size1() );
  W.resize ( invcov.size1(), invcov.size2() );
  matrix_double temp ( invcov );

  int nfound;

  boost::numeric::ublas::vector < int > support ( 2 * invcov.size1() );

  boost::numeric::bindings::lapack::syevr ( 'V', 'A', boost::numeric::bindings::lower ( temp ), 0.0, 0.0, 0, 0, 0.0, nfound, D, W, support );

  return;
}


void harp::eigen_compose ( eigen_op op, vector_double const & D, matrix_double const & W, matrix_double & out ) {

  double threshold = 1.0e-12;

  out.resize ( W.size1(), W.size2() );

  vector_double scaled ( D );

  // scale eigenvalues and apply threshold

  double max = 0.0;
  double min;
  double val;
  double invmin;

  for ( size_t i = 0; i < scaled.size(); ++i ) {
    val = scaled[i];
    if ( val > max ) {
      max = val;
    }
  }

  switch ( op ) {
    case EIG_SQRT:
      min = sqrt(max) * threshold;
      for ( size_t i = 0; i < scaled.size(); ++i ) {
        val = sqrt ( scaled[i] );
        val = ( val < min ) ? min : val;
        scaled[i] = val;
      }
      break;
    case EIG_INVSQRT:
      min = sqrt(max) * threshold;
      invmin = 1.0 / min;
      for ( size_t i = 0; i < scaled.size(); ++i ) {
        val = sqrt ( scaled[i] );
        val = ( val < min ) ? invmin : ( 1.0 / val );
        scaled[i] = val;
      }
      break;
    default:
      min = max * threshold;
      for ( size_t i = 0; i < scaled.size(); ++i ) {
        val = scaled[i];
        val = ( val < min ) ? min : val;
        scaled[i] = val;
      }
      break;
  }

  // Compute temp = W * op(D)

  matrix_double temp ( W );

  for ( size_t i = 0; i < temp.size2(); ++i ) {
    for ( size_t j = 0; j < temp.size1(); ++j ) {
      temp( j, i ) *= scaled[ i ];
    }
  }

  // compute out = ( W * op(D) ) * W^T

  boost::numeric::bindings::blas::gemm ( 1.0, temp, boost::numeric::bindings::trans ( W ), 0.0, out );

  return;
}


void harp::column_norm ( matrix_double const & mat, vector_double & S ) {

  S.resize( mat.size1() );
  S.clear();

  for ( size_t i = 0; i < mat.size2(); ++i ) {
    for ( size_t j = 0; j < mat.size1(); ++j ) {
      S[ j ] += mat( j, i );
    }
  }

  // Invert

  for ( size_t i = 0; i < S.size(); ++i ) {
    S[i] = 1.0 / S[i];
  }

  return;
}


void harp::apply_norm ( vector_double const & S, matrix_double & mat ) {

  for ( size_t i = 0; i < mat.size2(); ++i ) {
    for ( size_t j = 0; j < mat.size1(); ++j ) {
      mat( j, i ) *= S[j];
    }
  }

  return;
}


void harp::apply_inverse_norm ( vector_double const & S, matrix_double & mat ) {

  for ( size_t i = 0; i < mat.size2(); ++i ) {
    for ( size_t j = 0; j < mat.size1(); ++j ) {
      mat( j, i ) /= S[j];
    }
  }

  return;
}


void harp::apply_norm ( vector_double const & S, vector_double & vec ) {

  for ( size_t i = 0; i < vec.size(); ++i ) {
    vec[i] *= S[i];
  }

  return;
}


void harp::apply_inverse_norm ( vector_double const & S, vector_double & vec ) {

  for ( size_t i = 0; i < vec.size(); ++i ) {
    vec[i] /= S[i];
  }

  return;
}


void harp::norm ( vector_double const & D, matrix_double const & W, vector_double & S ) {

  matrix_double temp ( W );
  temp.clear();
 
  eigen_compose ( EIG_SQRT, D, W, temp );

  column_norm ( temp, S );
  
  return;
}





