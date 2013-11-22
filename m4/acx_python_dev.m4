#
# SYNOPSIS
#
#   ACX_PYTHON_DEV([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
#
# DESCRIPTION
#
#   This macro uses the python interpreter found by AM_PATH_PYTHON to
#   get the include and link flags.  The PYTHON_CONFIG, PYTHON_CPPFLAGS, 
#   PYTHON_LDFLAGS, and PYTHON_LIBS variables are set.
#
#   The user may use:
# 
#       --with-python-config=<executable>
#
#   to manually specify the location of the python-config executable.
#
#   ACTION-IF-FOUND is a list of shell commands to run if python is found,
#   and ACTION-IF-NOT-FOUND is a list of commands to run it if it is not found.
#   If ACTION-IF-FOUND is not specified, the default action will define 
#   HAVE_PYTHON and set output variables above.
#
#   This macro requires autoconf 2.50 or later.
#
# LAST MODIFICATION
#
#   2013-11-07
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

AC_DEFUN([ACX_PYTHON_DEV], [
AC_PREREQ(2.50)
AC_REQUIRE([AM_PATH_PYTHON])

acx_python_dev_ok=no

PYTHON_CONFIG="python-config"
PYTHON_CPPFLAGS=""
PYTHON_LDFLAGS=""
PYTHON_LIBS=""

AC_ARG_WITH(python-config, [AC_HELP_STRING([--with-python-config=<executable>], [use the speficied python-config executable.])])

if test x"$with_python_config" != x; then
   if test x"$with_python_config" = xno; then
      acx_python_dev_ok=disable
      PYTHON_CONFIG=""
   else
      if test x"$with_python_config" != xyes; then
         PYTHON_CONFIG="$with_python_config"
      fi
   fi
else
   PYTHON_CONFIG="$PYTHON-config"
fi
      
if test $acx_python_dev_ok = disable; then
   AC_MSG_NOTICE([**** python explicitly disabled by configure.])
else

   # print python interpreter
   AC_MSG_NOTICE([using python interpreter $PYTHON])
   AC_MSG_NOTICE([using python-config $PYTHON_CONFIG])

   # use python-config to set variables
   PYTHON_CPPFLAGS=`$PYTHON_CONFIG --includes`
   PYTHON_LDFLAGS=`$PYTHON_CONFIG --ldflags`
   PYTHON_LIBS=`$PYTHON_CONFIG --libs`

   # save environment
   acx_python_dev_save_CPPFLAGS="$CPPFLAGS"
   acx_python_dev_save_LIBS="$LIBS"

   # set compile variables
   CPPFLAGS="$acx_python_dev_save_CPPFLAGS $PYTHON_CPPFLAGS"
   LIBS="$acx_python_dev_save_LIBS $PYTHON_LDFLAGS $PYTHON_LIBS"

   # check for header
   AC_CHECK_HEADERS([Python.h])

   # check linking
   AC_MSG_CHECKING([python linking])
   AC_LINK_IFELSE([
      AC_LANG_PROGRAM([[#include <Python.h>]],
            [[Py_Initialize();]])
      ],[acx_python_dev_ok=yes;AC_DEFINE(HAVE_PYTHON,1,[Define if you have python development libraries.])])
   AC_MSG_RESULT([$acx_python_dev_ok])

   # reset environment
   CPPFLAGS="$acx_python_dev_save_CPPFLAGS"
   LIBS="$acx_python_dev_save_LIBS"

fi

AC_SUBST(PYTHON_CONFIG)
AC_SUBST(PYTHON_CPPFLAGS)
AC_SUBST(PYTHON_LDFLAGS)
AC_SUBST(PYTHON_LIBS)

# Execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:

if test x"$acx_python_dev_ok" = xyes; then
   ifelse([$1],,[echo "**** Enabling python support."],[$1])
else
   ifelse([$2],,[echo "**** python not found - disabling support."],[$2])
fi

])
