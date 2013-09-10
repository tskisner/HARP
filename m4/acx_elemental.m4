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
#       --with-elemental=<PATH>
#
#   to manually specify the Elemental installation prefix.
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
#   2013-03-09
#
# COPYING
#
#   Copyright (c) 2013 Theodore Kisner <tskisner@lbl.gov>
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
AC_REQUIRE([ACX_LAPACK])
AC_REQUIRE([AX_OPENMP])

acx_elemental_ok=no
acx_elemental_default="-lelemental -lpmrrr"

ELEMENTAL_CPPFLAGS=""
ELEMENTAL=""

AC_ARG_WITH(elemental, [AC_HELP_STRING([--with-elemental=<PATH>], [use the Elemental installed in <PATH>.])])

if test x"$with_elemental" != x; then
   if test x"$with_elemental" != xno; then
      ELEMENTAL_CPPFLAGS="-I$with_elemental/include"
      ELEMENTAL="-L$with_elemental/lib -lelemental -lpmrrr"
   else
      acx_elemental_ok=disable
   fi
fi

if test $acx_elemental_ok = disable; then
   echo "**** Elemental explicitly disabled by configure."
else

   # Save environment

   acx_elemental_save_CC="$CC"
   acx_elemental_save_CXX="$CXX"
   acx_elemental_save_CPP="$CPP"
   acx_elemental_save_CPPFLAGS="$CPPFLAGS"
   acx_elemental_save_LIBS="$LIBS"

   # Test MPI compile and linking.  We need the eigensolver capabilities,
   # so check that the library has been built with eigen routines.

   CXX="$MPICXX"
   CXXCPP="$MPICXX -E"
   CPPFLAGS="$CPPFLAGS $ELEMENTAL_CPPFLAGS"
   LIBS="$ELEMENTAL $acx_elemental_save_LIBS $LAPACK_LIBS $BLAS_LIBS $LIBS $FCLIBS -lm $OPENMP_CXXFLAGS"

   AC_CHECK_HEADERS([elemental.hpp])

   AC_MSG_CHECKING([for elem::HermitianEig in user specified location])
   AC_LINK_IFELSE([AC_LANG_PROGRAM([
     [#include <elemental.hpp>]
   ],[
      using namespace std;
      using namespace elem;
      elem::DistMatrix < double, elem::MC, elem::MR > cov ( 4, 4 );
      elem::DistMatrix < double, elem::MC, elem::MR > W ( 4, 4 );
      elem::DistMatrix < double, elem::VR, elem::STAR > eigvals ( 4, 1 );
      elem::HermitianEig ( elem::LOWER, cov, eigvals, W );
   ])],[acx_elemental_ok=yes;AC_DEFINE(HAVE_ELEMENTAL,1,[Define if you have the Elemental library.])])

   AC_MSG_RESULT($acx_elemental_ok)

   if test $acx_elemental_ok = no; then
      ELEMENTAL="$ELEMENTAL -llapack-addons"
      LIBS="$ELEMENTAL $acx_elemental_save_LIBS $LAPACK_LIBS $BLAS_LIBS $LIBS $FCLIBS -lm $OPENMP_CXXFLAGS"

      AC_MSG_CHECKING([for elem::HermitianEig in user specified location with lapack addons])
      AC_LINK_IFELSE([AC_LANG_PROGRAM([
        [#include <elemental.hpp>]
      ],[
         using namespace std;
         using namespace elem;
         elem::DistMatrix < double, elem::MC, elem::MR > cov ( 4, 4 );
         elem::DistMatrix < double, elem::MC, elem::MR > W ( 4, 4 );
         elem::DistMatrix < double, elem::VR, elem::STAR > eigvals ( 4, 1 );
         elem::HermitianEig ( elem::LOWER, cov, eigvals, W );
      ])],[acx_elemental_ok=yes;AC_DEFINE(HAVE_ELEMENTAL,1,[Define if you have the Elemental library.])])

      AC_MSG_RESULT($acx_elemental_ok)
   fi

   if test $acx_elemental_ok = no; then
      ELEMENTAL="$acx_elemental_default"
      LIBS="$ELEMENTAL $acx_elemental_save_LIBS -lm $OPENMP_CXXFLAGS"

      AC_MSG_CHECKING([for elem::HermitianEig in default location])
      AC_LINK_IFELSE([AC_LANG_PROGRAM([
        [#include <elemental.hpp>]
      ],[
         using namespace std;
         using namespace elem;
         elem::DistMatrix < double, elem::MC, elem::MR > cov ( 4, 4 );
         elem::DistMatrix < double, elem::MC, elem::MR > W ( 4, 4 );
         elem::DistMatrix < double, elem::VR, elem::STAR > eigvals ( 4, 1 );
         elem::HermitianEig ( elem::LOWER, cov, eigvals, W );
      ])],[acx_elemental_ok=yes;AC_DEFINE(HAVE_ELEMENTAL,1,[Define if you have the Elemental library.])])

      AC_MSG_RESULT($acx_elemental_ok)
   fi

   if test $acx_elemental_ok = no; then
      ELEMENTAL="$acx_elemental_default -llapack-addons"
      LIBS="$ELEMENTAL $acx_elemental_save_LIBS -lm $OPENMP_CXXFLAGS"

      AC_MSG_CHECKING([for elem::HermitianEig in default location with lapack addons])
      AC_LINK_IFELSE([AC_LANG_PROGRAM([
        [#include <elemental.hpp>]
      ],[
         using namespace std;
         using namespace elem;
         elem::DistMatrix < double, elem::MC, elem::MR > cov ( 4, 4 );
         elem::DistMatrix < double, elem::MC, elem::MR > W ( 4, 4 );
         elem::DistMatrix < double, elem::VR, elem::STAR > eigvals ( 4, 1 );
         elem::HermitianEig ( elem::LOWER, cov, eigvals, W );
      ])],[acx_elemental_ok=yes;AC_DEFINE(HAVE_ELEMENTAL,1,[Define if you have the Elemental library.])])

      AC_MSG_RESULT($acx_elemental_ok)
   fi

   if test $acx_elemental_ok = no; then
      ELEMENTAL=""
   fi

   # Restore environment

   CC="$acx_elemental_save_CC"
   CXX="$acx_elemental_save_CXX"
   CPP="$acx_elemental_save_CPP"
   CPPFLAGS="$acx_elemental_save_CPPFLAGS"
   LIBS="$acx_elemental_save_LIBS"

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
