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
dnl | Check OS.  Basically on Darwin we have to jump through hoops
dnl | to build loadable modules...
dnl +----------------------------------------
AC_CANONICAL_HOST
HARP_PLUGIN_COMPILE=""
HARP_PLUGIN_LINK=""
case $host_os in
  darwin* )
    HARP_PLUGIN_COMPILE="-bundle"
    HARP_PLUGIN_LINK="-flat_namespace -undefined suppress"
    HARP_PLUGIN_EXT="dylib"
    ;;
  linux*)
    HARP_PLUGIN_COMPILE="-fPIC"
    HARP_PLUGIN_LINK="-shared"
    HARP_PLUGIN_EXT="so"
    ;;
  *BSD*)
    HARP_PLUGIN_COMPILE="-fPIC"
    HARP_PLUGIN_LINK="-shared"
    HARP_PLUGIN_EXT="so"
    ;;
  *)
    HARP_PLUGIN_COMPILE="-fPIC"
    HARP_PLUGIN_LINK="-shared"
    HARP_PLUGIN_EXT="so"
    AC_MSG_WARN([Your platform was not detected- using Linux defaults for Plugin compilation])
    ;;
esac
dnl +----------------------------------------
dnl | Compiler features
dnl +----------------------------------------
dnl Set default optimization to '-O3' instead of the default '-g -O2'
if test x"${CXXFLAGS}" = x; then
   CXXFLAGS="-O3"
fi
AC_CANONICAL_HOST
AC_PROG_INSTALL
AC_PROG_LIBTOOL
LT_SYS_MODULE_EXT
LT_FUNC_DLSYM_USCORE
AC_PROG_CXX
AC_PROG_CC
AC_LANG([C++])
AC_CHECK_HEADERS([iostream])
AC_CHECK_HEADERS([iomanip])
AC_CHECK_HEADERS([cmath])
AC_CHECK_HEADERS([cstdlib])
harp_openmp=yes
ax_have_openmp=no
AC_ARG_ENABLE(openmp, [AC_HELP_STRING([--disable-openmp], [Disable OpenMP])], [harp_openmp=$enable_openmp])
if test x"$harp_openmp" != xno; then
  AX_OPENMP([ax_have_openmp=yes;AC_DEFINE(HAVE_OPENMP,1,[Define if OpenMP is enabled])])
fi
AM_CONDITIONAL([HAVE_AM_OPENMP], [test $ax_have_openmp = yes])
AC_SUBST([OPENMP_CXXFLAGS])
dnl +-------------------------
dnl | Should we enable expensive debugging checks in Boost?
dnl +-------------------------
harp_debug=no
AC_ARG_ENABLE([debug], [AC_HELP_STRING([--enable-debug], [Enable expensive range checking in Boost])])
if test x"$enable_debug" = x"yes"; then
  harp_debug=yes
fi
AM_CONDITIONAL([ENABLE_AM_DEBUG], [test "$harp_debug" = "yes"])
dnl +-------------------------
dnl | Option to perform full static linking of executables, not just
dnl | building of static libraries.
dnl +-------------------------
harp_static=no
harp_plugins=no
AC_ARG_ENABLE([all-static], [AC_HELP_STRING([--enable-all-static], [force static linking of executables.])])
if test x"$enable_all_static" = x"yes"; then
  if test x"$enable_static" != x"yes"; then
    AC_MSG_ERROR([Cannot perform static linking of executables if static libraries are disabled!])
  fi
  harp_static=yes
fi
AM_CONDITIONAL([USE_AM_ALLSTATIC], [test "$harp_static" = "yes"])
if test x"$enable_shared" = x"yes"; then
  if test x"$harp_static" = xno; then
    harp_plugins=yes
    AC_DEFINE([USE_PLUGINS], 1, [Enable support for external plugins])
  else
    AC_MSG_WARN([Executables are linked statically (--enable-all-static), so external Plugins are disabled])
    HARP_PLUGIN_COMPILE=""
    HARP_PLUGIN_LINK="-static"
    HARP_PLUGIN_EXT="a"
  fi
else
  AC_MSG_WARN([Shared libraries disabled (with --disable-shared), so external Plugins are disabled])
  HARP_PLUGIN_COMPILE=""
  HARP_PLUGIN_LINK="-static"
  HARP_PLUGIN_EXT="a"
fi
AC_SUBST(HARP_PLUGIN_COMPILE)
AC_SUBST(HARP_PLUGIN_LINK)
AC_SUBST(HARP_PLUGIN_EXT)
dnl
dnl =========================================
dnl Do these checks with the C compiler
dnl
AC_LANG_PUSH([C])
dnl +------------------------------------------------
dnl | Check for CFITSIO
dnl +------------------------------------------------
ACX_CFITSIO([], [AC_MSG_ERROR([Could not find the CFITSIO library!])])
dnl +------------------------------------------------
dnl | BLAS / LAPACK
dnl +------------------------------------------------
AX_CHECK_BLAS(,AC_MSG_ERROR([Could not find a BLAS library!]))
AX_CHECK_LAPACK(,AC_MSG_ERROR([Could not find a LAPACK library!]))
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
AC_LANG_POP([C])
dnl 
dnl End C language checks
dnl
dnl ==========================================
dnl
dnl +------------------------------------------------
dnl | Check for MPI
dnl +------------------------------------------------
harp_mpi=yes
AC_ARG_ENABLE([mpi], [AC_HELP_STRING([--disable-mpi], [Disable all checks for MPI])], [harp_mpi=$enable_mpi])
if test x"$harp_mpi" != xno; then
  ACX_MPI([harp_mpi=yes], [harp_mpi=no;AC_MSG_WARN([Could not find the MPI C++ compiler, MPI tools will not be built])])
fi
dnl +------------------------------------------------
dnl | Check for BOOST
dnl +------------------------------------------------
BOOST_REQUIRE([1.53], [AC_MSG_ERROR([Could not find BOOST library >= 1.53!])])
BOOST_REGEX([], [AC_MSG_ERROR([Could not find BOOST regex library])])
BOOST_SERIALIZATION([], [AC_MSG_ERROR([Could not find BOOST serialization library])])
BOOST_PROGRAM_OPTIONS([], [AC_MSG_ERROR([Could not find BOOST program options library])])
if test x"$harp_python" != xno; then
  BOOST_PYTHON([], [AC_MSG_WARN([Could not find boost::python, bindings will not be built])])
fi
if test x"$harp_mpi" != xno; then
  BOOST_MPI([harp_mpi=yes], [harp_mpi=no;AC_MSG_WARN([Could not find boost::mpi, MPI tools will not be built])])
fi
dnl +------------------------------------------------
dnl | Check for Elemental and any dependencies
dnl +------------------------------------------------
if test x"$harp_mpi" != xno; then
  ACX_ELEMENTAL([harp_mpi=yes], [harp_mpi=no;AC_MSG_WARN([Could not find Elemental library, MPI tools will not be built])])
fi
dnl +------------------------------------------------------------------------
dnl | Set automake variables which control building of python and MPI tools
dnl +------------------------------------------------------------------------
AM_CONDITIONAL([HAVE_AM_PYTHON], [test "x$PYTHON_VERSION" != x && test "x$BOOST_PYTHON_LIBS" != x])
AM_CONDITIONAL([HAVE_AM_MPI], [test "x$harp_mpi" = xyes])
dnl +------------------------------------------------
dnl | Set outputs
dnl +------------------------------------------------
AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([src/Makefile])
AC_CONFIG_FILES([src/libharp/math/Makefile])
AC_CONFIG_FILES([src/libharp/data/Makefile])
AC_CONFIG_FILES([src/libharp/plugins/Makefile])
AC_CONFIG_FILES([src/libharp/extract/Makefile])
AC_CONFIG_FILES([src/libharp/Makefile])
AC_CONFIG_FILES([src/tests/Makefile])
AC_CONFIG_FILES([src/apps/Makefile])
dnl AC_CONFIG_FILES([src/python/Makefile])
dnl AC_CONFIG_FILES([src/tests-python/Makefile])
AC_CONFIG_FILES([src/libharp-mpi/math/Makefile])
AC_CONFIG_FILES([src/libharp-mpi/plugins/Makefile])
AC_CONFIG_FILES([src/libharp-mpi/Makefile])
AC_CONFIG_FILES([src/tests-mpi/Makefile])
AC_CONFIG_FILES([src/apps-mpi/Makefile])
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
if test x"${harp_mpi}" = xyes; then
  AC_MSG_NOTICE([  MPI C++ Compiler          : $MPICXX])
else
  AC_MSG_NOTICE([  MPI C++ Compiler          : Not Found])
fi
if test x"${harp_debug}" = xyes; then
  AC_MSG_NOTICE([  C++ Compile flags         : $CXXFLAGS -DNDEBUG])
else
  AC_MSG_NOTICE([  C++ Compile flags         : $CXXFLAGS])
fi
if test x"${ax_have_openmp}" = xyes; then
  AC_MSG_NOTICE([  OpenMP                    : Detected ($OPENMP_CXXFLAGS)])
else
  AC_MSG_NOTICE([  OpenMP                    : Disabled])
fi
AC_MSG_NOTICE([  CFITSIO Compile           : $CFITSIO_CPPFLAGS])
AC_MSG_NOTICE([  CFITSIO Link              : $CFITSIO])
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
if test x"${harp_plugins}" = xyes; then
  AC_MSG_NOTICE([  External Plugins          : ENABLED])
else
  AC_MSG_NOTICE([  External Plugins          : DISABLED])
fi
AC_MSG_NOTICE([                                           ])
