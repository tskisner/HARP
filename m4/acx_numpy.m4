#
# SYNOPSIS
#
#   ACX_NUMPY([])
#
# DESCRIPTION
#
#   This macro checks for NumPy headers and sets the NUMPY_CPPFLAGS
#   output variable. 
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

AC_DEFUN([ACX_NUMPY],[
  AC_REQUIRE([ACX_PYTHON_DEV])

  acx_numpy_ok="no"

  if test "x$PYTHON" != "x"; then
    AC_MSG_CHECKING([for NumPy])
    cat > conftest.py << EOF
import sys
try:
  import numpy
except ImportError:
  sys.exit(1)
EOF
    if $PYTHON conftest.py > /dev/null 2>&1; then
      acx_numpy_ok="yes"
    fi
    AC_MSG_RESULT([$acx_numpy_ok])
  fi

  if test "x$acx_numpy_ok" = "xyes"; then
    AC_MSG_CHECKING([NumPy includes])
    NUMPY_CPPFLAGS=-I`$PYTHON -c "import numpy; print numpy.get_include()"`
    AC_MSG_RESULT([$NUMPY_CPPFLAGS])
    AC_SUBST([NUMPY_CPPFLAGS])

    acx_numpy_save_CPPFLAGS=$CPPFLAGS

    CPPFLAGS="$CPPFLAGS $PYTHON_CPPFLAGS $NUMPY_CPPFLAGS"
    AC_CHECK_HEADERS([numpy/arrayobject.h],,[acx_numpy_ok="no"],[
      #include<Python.h>
    ])

    CPPFLAGS=$acx_numpy_save_CPPFLAGS
  fi

  if test "x$acx_numpy_ok" = "xyes"; then
    AC_DEFINE([HAVE_NUMPY], [1], [Define if we have NumPy header files])
  fi

])


