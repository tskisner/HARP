
# Initialize compile and linking flags

AM_CPPFLAGS =
AM_CXXFLAGS =
AM_LDFLAGS =
AM_LIBS = 

# Append compile and linking flags based on configure variables
include $(top_srcdir)/src/build_options.am


EXTRA_DIST = gen_python.sh

pkgharpdir = $(pyexecdir)/harp

pkgharp_PYTHON = __init__.py




#noinst_HEADERS = pytoast.hpp

#pkgtoastcoredir = $(pyexecdir)/pytoast/core

#pkgtoastcore_LTLIBRARIES = _pytoast.la
#_pytoast_la_LDFLAGS = -module
#_pytoast_la_LIBADD = $(top_builddir)/src/libtoast/libtoast.la $(AM_LIBS)
#_pytoast_la_SOURCES = pytoast.cpp

#pkgtoastcore_PYTHON = core/__init__.py

#pkgtoastdir = $(pyexecdir)/pytoast
#pkgtoast_PYTHON = __init__.py



