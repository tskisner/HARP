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
dnl +---------------------
dnl | Compiler features
dnl +---------------------
dnl Set default optimization to '-O3' instead of the default '-g -O2'
if test x"${CXXFLAGS}" = x; then
   CXXFLAGS="-O3"
fi
AC_CANONICAL_HOST
AC_PROG_INSTALL
AC_PROG_LIBTOOL
AC_PROG_CXX
AC_LANG([C++])
AC_CHECK_HEADERS([iostream iostream.h])
AC_CHECK_HEADERS([iomanip iomanip.h])
AC_CHECK_HEADERS([cmath math.h])
AC_CHECK_HEADERS([cstdlib stdlib.h])
dnl +------------------------------------------------
dnl | Check for OpenMP
dnl +------------------------------------------------
harp_openmp=yes
ax_have_openmp=no
AC_ARG_ENABLE(openmp, [AC_HELP_STRING([--disable-openmp], [Disable OpenMP])], [harp_openmp=$enable_openmp])
if test x"$harp_openmp" != xno; then
  AX_OPENMP([ax_have_openmp=yes;AC_DEFINE(HAVE_OPENMP,1,[Define if OpenMP is enabled])])
fi
AM_CONDITIONAL([HAVE_AM_OPENMP], [test $ax_have_openmp = yes])
dnl +------------------------------------------------
dnl | Check for MPI
dnl +------------------------------------------------
ACX_MPI([], [AC_MSG_ERROR([Could not find the MPI C++ compiler])])
dnl +------------------------------------------------
dnl | Check for Python
dnl +------------------------------------------------
harp_python=yes
AC_ARG_ENABLE(python, [AC_HELP_STRING([--disable-python], [Disable python bindings])], [harp_python=$enable_python])
if test x"$harp_python" != xno; then
  AM_PATH_PYTHON
  AX_PYTHON_DEVEL
  ACX_NUMPY
fi
AM_CONDITIONAL([HAVE_AM_PYTHON], [test "x$PYTHON_VERSION" != x])
dnl +------------------------------------------------
dnl | Check for BOOST
dnl +------------------------------------------------
AX_BOOST_BASE([], [], [AC_MSG_ERROR([Could not find the BOOOST library])])
AX_BOOST_PROGRAM_OPTIONS
dnl +------------------------------------------------
dnl | Check for Clique and any dependencies
dnl +------------------------------------------------
ACX_CLIQUE([], [AC_MSG_ERROR([Could not find the Clique library!])])
dnl +------------------------------------------------
dnl | Check for CFITSIO
dnl +------------------------------------------------
ACX_CFITSIO([], [AC_MSG_ERROR([Could not find the CFITSIO library!])])
dnl +------------------------------------------------
dnl | Check for Doxygen
dnl +------------------------------------------------
DOXYGEN=""
AC_CHECK_PROGS([DOXYGEN], [doxygen])
if test x"$DOXYGEN" = x; then
   AC_MSG_WARN([Doxygen not found - continuing without Doxygen support])
fi
AM_CONDITIONAL([HAVE_AM_DOXYGEN], [test x"$DOXYGEN" != x])
dnl +------------------------------------------------
dnl | Format config
dnl +------------------------------------------------
format_list=""
m4_include([instruments.m4]) 
dnl +------------------------------------------------
dnl | Set outputs
dnl +------------------------------------------------
AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([src/Makefile])
AC_CONFIG_FILES([src/libharp/Makefile])
AC_CONFIG_FILES([src/tests/Makefile])
AC_CONFIG_FILES([src/apps/Makefile])
dnl AC_CONFIG_FILES([src/python/Makefile])
dnl AC_CONFIG_FILES([src/tests-python/Makefile])
dnl AC_CONFIG_FILES([docs/Makefile])
dnl AM_COND_IF([HAVE_AM_DOXYGEN], [AC_CONFIG_FILES([docs/toast.dox])])
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
AC_MSG_NOTICE([  MPICXX Compiler           : $MPICXX])
AC_MSG_NOTICE([  C++ Compile flags         : $CXXFLAGS])
if test x"${ax_have_openmp}" = xyes; then
  AC_MSG_NOTICE([  OpenMP                    : Detected ($OPENMP_CXXFLAGS)])
else
  AC_MSG_NOTICE([  OpenMP                    : Disabled])
fi
AC_MSG_NOTICE([  Boost Compile             : $BOOST_CPPFLAGS])
AC_MSG_NOTICE([  Boost Link Flags          : $BOOST_LDFLAGS])
AC_MSG_NOTICE([  Boost Program Options     : $BOOST_PROGRAM_OPTIONS_LIB])
AC_MSG_NOTICE([  Clique Compile         : $CLIQUE_CPPFLAGS])
AC_MSG_NOTICE([  Clique Link            : $CLIQUE])
AC_MSG_NOTICE([  CFITSIO Compile           : $CFITSIO_CPPFLAGS])
AC_MSG_NOTICE([  CFITSIO Link              : $CFITSIO])
if test x"${PYTHON_VERSION}" != x; then
  AC_MSG_NOTICE([  Python                    : $PYTHON])
else
  AC_MSG_NOTICE([  Python                    : Disabled])
fi
AC_MSG_NOTICE([                                           ])
AC_MSG_NOTICE([============== Format Status ==============])
AC_MSG_NOTICE([                                           ])
AC_MSG_NOTICE([  Formats Enabled:                         ])
AC_MSG_NOTICE([  $format_list ])
AC_MSG_NOTICE([                                           ])
AC_MSG_NOTICE([===========================================])
AC_MSG_NOTICE([                                           ])
