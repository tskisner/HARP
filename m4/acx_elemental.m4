#
# SYNOPSIS
#
#   ACX_ELEMENTAL([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
#
# DESCRIPTION
#
#   This macro looks for a version of the Elemental library.  The ELEMENTAL_CPPFLAGS
#   and ELEMENTAL output variables hold the compile and link flags.
#
#   To link an application with Elemental, you should link with:
#
#   	$ELEMENTAL
#
#   The user may use:
# 
#       --with-elemental-cpp=<flags> --with-elemental-libs=<flags> 
#
#   to manually specify the Elemental include and linking flags.
#
#   ACTION-IF-FOUND is a list of shell commands to run if an Elemental library is
#   found, and ACTION-IF-NOT-FOUND is a list of commands to run it if it is
#   not found. If ACTION-IF-FOUND is not specified, the default action will
#   define HAVE_ELEMENTAL and set output variables above.
#
#   This macro requires autoconf 2.50 or later.
#
# LAST MODIFICATION
#
#   2012-10-01
#
# COPYING
#
#   Copyright (c) 2012 Theodore Kisner <tskisner@lbl.gov>
#
#   All rights reserved.
#
#   Redistribution and use in source and binary forms, with or without modification,
#   are permitted provided that the following conditions are met:
#
#   o  Redistributions of source code must retain the above copyright notice, 
#      this list of conditions and the following disclaimer.
#
#   o  Redistributions in binary form must reproduce the above copyright notice, 
#      this list of conditions and the following disclaimer in the documentation
#      and/or other materials provided with the distribution.
#
#   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND 
#   ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED 
#   WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. 
#   IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, 
#   INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, 
#   BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
#   DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
#   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE 
#   OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
#   OF THE POSSIBILITY OF SUCH DAMAGE.
#

AC_DEFUN([ACX_ELEMENTAL], [
AC_PREREQ(2.50)
AC_REQUIRE([ACX_MPI])

acx_elemental_ok=no
acx_elemental_default="-lelemental -lplcg -lpmrrr -lcmake-dummy-lib -llapack -lblas"

ELEMENTAL_CPPFLAGS=""
ELEMENTAL=""

AC_ARG_WITH(elemental-cpp, [AC_HELP_STRING([--with-elemental-cpp=<flags>], [use Elemental preprocessing flags <flags>.  Set to "no" to disable.])])

AC_ARG_WITH(elemental-libs, [AC_HELP_STRING([--with-elemental-libs=<flags>], [use Elemental linking flags <flags>.  Set to "no" to disable.])])

if test x"$with_elemental_cpp" != x; then
   if test x"$with_elemental_cpp" != xno; then
      ELEMENTAL_CPPFLAGS="$with_elemental_cpp"
   else
      acx_elemental_ok=disable
   fi
fi

if test x"$with_elemental_libs" != x; then
   if test x"$with_elemental_libs" != xno; then
      ELEMENTAL="$with_elemental_libs"
   else
      acx_elemental_ok=disable
   fi
fi

if test $acx_elemental_ok = disable; then
   echo "**** Elemental explicitly disabled by configure."
else

   # Save environment

   acx_elemental_save_CXX="$CXX"
   acx_elemental_save_CPP="$CPP"
   acx_elemental_save_CPPFLAGS="$CPPFLAGS"
   acx_elemental_save_LIBS="$LIBS"

   # Test MPI compile and linking.  We need the eigensolver capabilities,
   # so check that the library has been built with eigen routines.

   CXX="$MPICXX"
   CXXCPP="$MPICXX -E"
   CPPFLAGS="$CPPFLAGS $ELEMENTAL_CPPFLAGS"
   LIBS="$ELEMENTAL $acx_elemental_save_LIBS -lm $OPENMP_CXXFLAGS"

   AC_CHECK_HEADERS([elemental.hpp])

   AC_MSG_CHECKING([for elem::SortEig in user specified location])
   AC_LINK_IFELSE([AC_LANG_PROGRAM([
     [#include <elemental.hpp>]
   ],[
      using namespace std;
      using namespace elem;
      typedef double R;
      typedef Complex<R> C;
      int m = 10;
      int n = 10;
      DistMatrix<C> A;
      Uniform( m, n, A );
      DistMatrix<C> U, V;
      DistMatrix<R,VR,STAR> s;
      U = A;
      SVD( U, s, V );
   ])],[acx_elemental_ok=yes;AC_DEFINE(HAVE_ELEMENTAL,1,[Define if you have the Elemental library.])])

   AC_MSG_RESULT($acx_elemental_ok)

   if test $acx_elemental_ok = no; then
      ELEMENTAL="$acx_elemental_default"
      LIBS="$acx_elemental_default $acx_elemental_save_LIBS -lm $OPENMP_CXXFLAGS"

      AC_MSG_CHECKING([for elem::SortEig in default location])
      AC_LINK_IFELSE([AC_LANG_PROGRAM([
        [#include <elemental.hpp>]
      ],[
         using namespace std;
         using namespace elem;
         typedef double R;
         typedef Complex<R> C;
         int m = 10;
         int n = 10;
         DistMatrix<C> A;
         Uniform( m, n, A );
         DistMatrix<C> U, V;
         DistMatrix<R,VR,STAR> s;
         U = A;
         SVD( U, s, V );
      ])],[acx_elemental_ok=yes;AC_DEFINE(HAVE_ELEMENTAL,1,[Define if you have the Elemental library.])])

      AC_MSG_RESULT($acx_elemental_ok)
   fi

   if test $acx_elemental_ok = no; then
      ELEMENTAL=""
   fi

   # Restore environment

   CC="$acx_elemental_save_CC"
   CPP="$acx_elemental_save_CPP"
   LIBS="$acx_elemental_save_LIBS"
   CPPFLAGS="$acx_elemental_save_CPPFLAGS"

fi

# Define exported variables

AC_SUBST(ELEMENTAL_CPPFLAGS)
AC_SUBST(ELEMENTAL)

# Execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
   
if test x"$acx_elemental_ok" = xyes; then
   ifelse([$1],,[echo "**** Enabling support for Elemental."],[$1])
else
   ifelse([$2],,[echo "**** Elemental not found - disabling support."],[$2])
fi

])
