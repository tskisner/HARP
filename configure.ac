dnl
dnl This is the Autoconf file for the HARP
dnl
dnl +------------------------
dnl | Initialize package info
dnl +------------------------
AC_INIT(Makefile.am)
AC_CONFIG_HEADERS(config.h)
AM_INIT_AUTOMAKE(harp, 0.0.1)
AC_CONFIG_MACRO_DIR([m4])
dnl +-------------------------
dnl | If the library source code has changed at all since the last update, 
dnl | then increment revision (c:r:a becomes c:r+1:a).  If any interfaces 
dnl | have been added, removed, or changed since the last update, increment
dnl | current, and set revision to 0.  If any interfaces have been added 
dnl | since the last public release, then increment age.  If any interfaces 
dnl | have been removed since the last public release, then set age to 0.
dnl +-------------------------
HARP_LT_VERSION="0:0:0"
AC_SUBST(HARP_LT_VERSION)
dnl +----------------------------------------
dnl | Compiler features and BLAS / LAPACK
dnl +----------------------------------------
dnl Set default optimization to '-O3' instead of the default '-g -O2'
if test x"${CXXFLAGS}" = x; then
   CXXFLAGS="-O3"
fi
AC_CANONICAL_HOST
AC_PROG_INSTALL
AC_PROG_LIBTOOL
dnl +-------------------------
dnl | Option to perform full static linking of executables, not just
dnl | building of static libraries.
dnl +-------------------------
harp_static=no
AC_ARG_ENABLE([all-static], [AC_HELP_STRING([--enable-all-static], [force static linking of executables.])])
if test x"$enable_all_static" = x"yes"; then
  if test x"$enable_static" != x"yes"; then
    AC_MSG_ERROR([Cannot perform static linking of executables if static libraries are disabled!])
  fi
  harp_static=yes
fi
AM_CONDITIONAL([USE_AM_ALLSTATIC], [test "$harp_static" = "yes"])
harp_fortran=yes
AC_ARG_ENABLE(fortran, [AC_HELP_STRING([--disable-fortran], [Disable features that depend on fortran])], [harp_fortran=$enable_fortran])
if test x"$harp_fortran" != xno; then
  AC_PROG_FC
dnl  AC_FC_LIBRARY_LDFLAGS
dnl  AC_FC_WRAPPERS
dnl  ACX_FC_PP_DEFINE
fi
AC_PROG_CXX
AC_LANG([C++])
AC_CHECK_HEADERS([iostream iostream.h])
AC_CHECK_HEADERS([iomanip iomanip.h])
AC_CHECK_HEADERS([cmath math.h])
AC_CHECK_HEADERS([cstdlib stdlib.h])
harp_openmp=yes
ax_have_openmp=no
AC_ARG_ENABLE(openmp, [AC_HELP_STRING([--disable-openmp], [Disable OpenMP])], [harp_openmp=$enable_openmp])
if test x"$harp_openmp" != xno; then
  AX_OPENMP([ax_have_openmp=yes;AC_DEFINE(HAVE_OPENMP,1,[Define if OpenMP is enabled])])
fi
AM_CONDITIONAL([HAVE_AM_OPENMP], [test $ax_have_openmp = yes])
dnl +------------------------------------------------
dnl | BLAS / LAPACK
dnl +------------------------------------------------
AX_CHECK_BLAS(,AC_MSG_ERROR([Could not find a BLAS library!]))
AX_CHECK_LAPACK(,AC_MSG_ERROR([Could not find a LAPACK library!]))
dnl +------------------------------------------------
dnl | Check for CFITSIO
dnl +------------------------------------------------
ACX_CFITSIO([], [AC_MSG_ERROR([Could not find the CFITSIO library!])])
dnl +------------------------------------------------
dnl | Check for Python
dnl +------------------------------------------------
harp_python=yes
AC_ARG_ENABLE(python, [AC_HELP_STRING([--disable-python], [Disable python bindings])], [harp_python=$enable_python])
if test x"$harp_python" != xno; then
  AM_PATH_PYTHON
  ACX_PYTHON_DEV
  ACX_NUMPY
fi
dnl +------------------------------------------------
dnl | Check for MPI
dnl +------------------------------------------------
have_mpi=no
ACX_MPI([have_mpi=yes], [AC_MSG_WARN([Could not find the MPI C++ compiler, MPI tools will not be built])])
dnl +------------------------------------------------
dnl | Check for BOOST
dnl +------------------------------------------------
BOOST_REQUIRE([1.53], [AC_MSG_ERROR([Could not find BOOST library >= 1.53!])])
BOOST_SERIALIZATION([], [AC_MSG_ERROR([Could not find BOOST serialization library])])
BOOST_PROGRAM_OPTIONS([], [AC_MSG_ERROR([Could not find BOOST program options library])])
if test x"$harp_python" != xno; then
  BOOST_PYTHON([], [AC_MSG_WARN([Could not find boost::python, bindings will not be built])])
fi
if test x"${MPICXX}" != x; then
  BOOST_MPI([], [AC_MSG_WARN([Could not find boost::mpi, MPI tools will not be built])])
fi
dnl +------------------------------------------------
dnl | Check for Elemental and any dependencies
dnl +------------------------------------------------
harp_elemental=no
if test x"$harp_fortran" != xno; then
  ACX_ELEMENTAL([harp_elemental=yes], [AC_MSG_WARN([Could not find Elemental library, MPI tools will not be built])])
fi
dnl +------------------------------------------------------------------------
dnl | Set automake variables which control building of python and MPI tools
dnl +------------------------------------------------------------------------
AM_CONDITIONAL([HAVE_AM_PYTHON], [test "x$PYTHON_VERSION" != x && test "x$BOOST_PYTHON_LIBS" != x])
AM_CONDITIONAL([HAVE_AM_MPI], [test "x$MPICXX" != x && test "x$BOOST_MPI_LIBS" != x && test "x$ELEMENTAL" != x])
dnl +------------------------------------------------
dnl | Set outputs
dnl +------------------------------------------------
AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([src/Makefile])
AC_CONFIG_FILES([src/libharp/math/Makefile])
AC_CONFIG_FILES([src/libharp/Makefile])
AC_CONFIG_FILES([src/tests/Makefile])
AC_CONFIG_FILES([src/apps/Makefile])
dnl AC_CONFIG_FILES([src/python/Makefile])
dnl AC_CONFIG_FILES([src/tests-python/Makefile])
AC_CONFIG_FILES([docs/Makefile])
dnl +------------------------------------------------
dnl | Generate Makefiles
dnl +------------------------------------------------
AC_OUTPUT
dnl +------------------------------------------------
dnl | Print out detected build options
dnl +------------------------------------------------
AC_MSG_NOTICE([                                           ])
AC_MSG_NOTICE([=========== Build Configuration ===========])
AC_MSG_NOTICE([                                           ])
AC_MSG_NOTICE([  C++ Compiler              : $CXX])
if test x"${have_mpi}" = xyes; then
  AC_MSG_NOTICE([  MPI C++ Compiler          : $MPICXX])
else
  AC_MSG_NOTICE([  MPI C++ Compiler          : Not Found])
fi
AC_MSG_NOTICE([  C++ Compile flags         : $CXXFLAGS])
if test x"${ax_have_openmp}" = xyes; then
  AC_MSG_NOTICE([  OpenMP                    : Detected ($OPENMP_CXXFLAGS)])
else
  AC_MSG_NOTICE([  OpenMP                    : Disabled])
fi
AC_MSG_NOTICE([  CFITSIO Compile           : $CFITSIO_CPPFLAGS])
AC_MSG_NOTICE([  CFITSIO Link              : $CFITSIO])
AC_MSG_NOTICE([  CBLAS Link                : $CBLAS_LIBS])
AC_MSG_NOTICE([  BLAS Link                 : $BLAS_LIBS])
AC_MSG_NOTICE([  LAPACK Link               : $LAPACK_LIBS])
if test x"${PYTHON_VERSION}" != x; then
  AC_MSG_NOTICE([  Python                    : $PYTHON])
else
  AC_MSG_NOTICE([  Python                    : Disabled])
fi
AC_MSG_NOTICE([  Boost Compile             : $BOOST_CPPFLAGS])
AC_MSG_NOTICE([  Boost Program Options     : $BOOST_PROGRAM_OPTIONS_LDFLAGS $BOOST_PROGRAM_OPTIONS_LIBS])
AC_MSG_NOTICE([  Boost Serialization       : $BOOST_SERIALIZATION_LDFLAGS $BOOST_SERIALIZATION_LIBS])
if test x"$BOOST_MPI_LIBS" != x; then
  AC_MSG_NOTICE([  Boost MPI                 : $BOOST_MPI_LDFLAGS $BOOST_MPI_LIBS])
else
  AC_MSG_NOTICE([  Boost MPI                 : Disabled])
fi
if test x"${acx_elemental_ok}" = xyes; then
  AC_MSG_NOTICE([  Elemental Compile         : $ELEMENTAL_CPPFLAGS])
  AC_MSG_NOTICE([  Elemental Link            : $ELEMENTAL])
else
  AC_MSG_NOTICE([  Elemental                 : Disabled])
fi
AC_MSG_NOTICE([                                           ])
