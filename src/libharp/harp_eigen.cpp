// @COPYRIGHT@

#include <harp_internal.hpp>


using namespace std;
using namespace harp;


extern "C" void dsyevd_ ( char *, char *, int *, double *, int *, double *, double *, int *, int *, int *, int * );

int harp::lapack_ev ( int dim, double * mat, double * eigen ) {
  char jobz = 'V';
  char uplo = 'L';
  int lwork;
  int liwork;
  double * work;
  int * iwork;
  int info;

  lwork = 1 + 6 * dim + 2 * dim * dim;
  work = moat::double_alloc ( lwork );

  liwork = 3 + 5 * dim;
  iwork = moat::int_alloc ( liwork );

  dsyevd_ ( &jobz, &uplo, &dim, mat, &dim, eigen, work, &lwork, iwork, &liwork, &info );

  free ( iwork );
  free ( work );

  return info;
}


