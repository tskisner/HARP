#
# SYNOPSIS
#
#   ACX_NUMPY([])
#
# DESCRIPTION
#
#   This macro checks for NumPy and sets the $(NUMPY_CPPFLAGS) 
#   output variable. 
#
# LICENSE
#
#   Copyright (c) 2011 Ted Kisner <tskisner.public@gmail.com>
#
#   This program is free software: you can redistribute it and/or modify it
#   under the terms of the GNU General Public License as published by the
#   Free Software Foundation, either version 3 of the License, or (at your
#   option) any later version.
#
#   This program is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
#   Public License for more details.
#
#   You should have received a copy of the GNU General Public License along
#   with this program. If not, see <http://www.gnu.org/licenses/>.
#
#   As a special exception, the respective Autoconf Macro's copyright owner
#   gives unlimited permission to copy, distribute and modify the configure
#   scripts that are the output of Autoconf when processing the Macro. You
#   need not follow the terms of the GNU General Public License when using
#   or distributing such scripts, even though portions of the text of the
#   Macro appear in them. The GNU General Public License (GPL) does govern
#   all other use of the material that constitutes the Autoconf Macro.
#
#   This special exception to the GPL applies to versions of the Autoconf
#   Macro released by the Autoconf Archive. When you make and distribute a
#   modified version of the Autoconf Macro, you may extend this special
#   exception to the GPL to apply to your modified version as well.

AC_DEFUN([ACX_NUMPY],[

  AC_PATH_PROG([PYTHON],[python[$PYTHON_VERSION]])
  if test -z "$PYTHON"; then
    AC_MSG_ERROR([Cannot find python$PYTHON_VERSION in your system path])
    PYTHON_VERSION=""
  fi

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

    saved_cppflags=$CPPFLAGS
    CPPFLAGS="${CPPFLAGS} ${PYTHON_CPPFLAGS} ${NUMPY_CPPFLAGS}"
    AC_CHECK_HEADERS([numpy/arrayobject.h],,[acx_numpy_ok="no"],[
#include<Python.h>
])
    CPPFLAGS=$saved_cppflags
  fi

  if test "x$acx_numpy_ok" = "xyes"; then
    AC_DEFINE([HAVE_NUMPY], [], [ Define to to enable NumPy support in the Python bindings ])
  fi

])


