#
# SYNOPSIS
#
#   ACX_CLIQUE([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
#
# DESCRIPTION
#
#   This macro looks for a version of the Clique library.  The CLIQUE_CPPFLAGS
#   and CLIQUE output variables hold the compile and link flags.
#
#   To link an application with Clique, you should link with:
#
#   	$CLIQUE
#
#   The user may use:
# 
#       --with-clique=<path> 
#
#   to manually specify the installation prefix of Clique.
#
#   ACTION-IF-FOUND is a list of shell commands to run if a clique library is
#   found, and ACTION-IF-NOT-FOUND is a list of commands to run it if it is
#   not found. If ACTION-IF-FOUND is not specified, the default action will
#   define HAVE_CLIQUE and set output variables above.
#
#   This macro requires autoconf 2.50 or later.
#
# LAST MODIFICATION
#
#   2013-08-14
#
# COPYING
#
#   Copyright (c) 2012-2013 Theodore Kisner <tskisner@lbl.gov>
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

AC_DEFUN([ACX_CLIQUE], [
AC_PREREQ(2.50)
AC_REQUIRE([ACX_ELEMENTAL])

acx_clique_ok=no
acx_clique_default="-lclique -lparmetis-addons -lmetis-addons -lparmetis -lmetis"

CLIQUE_CPPFLAGS=""
CLIQUE=""

AC_ARG_WITH(clique, [AC_HELP_STRING([--with-clique=<PATH>], [use the Clique installed in <PATH>.])])

if test x"$with_elemental" != x; then
   if test x"$with_elemental" != xno; then
      ELEMENTAL_CPPFLAGS="-I$with_elemental/include"
      ELEMENTAL="-L$with_elemental/lib -lelemental -lpmrrr"
   else
      acx_elemental_ok=disable
   fi
fi

AC_ARG_WITH(clique-cpp, [AC_HELP_STRING([--with-clique-cpp=<flags>], [use Clique preprocessing flags <flags>.  Set to "no" to disable.])])

AC_ARG_WITH(clique-libs, [AC_HELP_STRING([--with-clique-libs=<flags>], [use Clique linking flags <flags>.  Set to "no" to disable.])])

if test x"$with_clique_cpp" != x; then
   if test x"$with_clique_cpp" != xno; then
      CLIQUE_CPPFLAGS="$with_clique_cpp"
   else
      acx_clique_ok=disable
   fi
fi

if test x"$with_clique_libs" != x; then
   if test x"$with_clique_libs" != xno; then
      CLIQUE="$with_clique_libs"
   else
      acx_clique_ok=disable
   fi
fi

if test $acx_clique_ok = disable; then
   echo "**** Clique explicitly disabled by configure."
else

   # Save environment

   acx_clique_save_CXX="$CXX"
   acx_clique_save_CPP="$CPP"
   acx_clique_save_CPPFLAGS="$CPPFLAGS"
   acx_clique_save_LIBS="$LIBS"

   # Test MPI compile and linking.  We need the eigensolver capabilities,
   # so check that the library has been built with eigen routines.

   CXX="$MPICXX"
   CXXCPP="$MPICXX -E"
   CPPFLAGS="$CPPFLAGS $CLIQUE_CPPFLAGS"
   LIBS="$CLIQUE $acx_clique_save_LIBS -lm $OPENMP_CXXFLAGS"

   AC_CHECK_HEADERS([clique.hpp])

   AC_MSG_CHECKING([for cliq::DistSparseMatrix::StartAssembly in user specified location])
   AC_LINK_IFELSE([AC_LANG_PROGRAM([
     [#include <clique.hpp>]
   ],[
      using namespace std;
      using namespace cliq;
      DistSparseMatrix < double > test;
      test.StartAssembly();
   ])],[acx_clique_ok=yes;AC_DEFINE(HAVE_CLIQUE,1,[Define if you have the Clique library.])])

   AC_MSG_RESULT($acx_clique_ok)

   if test $acx_clique_ok = no; then
      CLIQUE="$acx_clique_default"
      LIBS="$acx_clique_default $LAPACK_LIBS $BLAS_LIBS $acx_clique_save_LIBS -lm $OPENMP_CXXFLAGS"

      AC_MSG_CHECKING([for cliq::DistSparseMatrix::StartAssembly in default location])
      AC_LINK_IFELSE([AC_LANG_PROGRAM([
        [#include <clique.hpp>]
      ],[
         using namespace std;
         using namespace cliq;
         DistSparseMatrix < double > test;
         test.StartAssembly();
      ])],[acx_clique_ok=yes;AC_DEFINE(HAVE_CLIQUE,1,[Define if you have the Clique library.])])

      AC_MSG_RESULT($acx_clique_ok)
   fi

   if test $acx_clique_ok = no; then
      CLIQUE=""
   fi

   # Restore environment

   CC="$acx_clique_save_CC"
   CPP="$acx_clique_save_CPP"
   LIBS="$acx_clique_save_LIBS"
   CPPFLAGS="$acx_clique_save_CPPFLAGS"

fi

# Define exported variables

AC_SUBST(CLIQUE_CPPFLAGS)
AC_SUBST(CLIQUE)

# Execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
   
if test x"$acx_clique_ok" = xyes; then
   ifelse([$1],,[echo "**** Enabling support for Clique."],[$1])
else
   ifelse([$2],,[echo "**** Clique not found - disabling support."],[$2])
fi

])
