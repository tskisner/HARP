SUBDIRS = boost_numpy

AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = -I$(top_srcdir)/src/libharp -I$(top_srcdir)/src/libharp/math -I$(top_srcdir)/src/libharp/data -I$(top_srcdir)/src/libharp/plugin $(CFITSIO_CPPFLAGS) $(BOOST_CPPFLAGS) $(PYTHON_CPPFLAGS) $(NUMPY_CPPFLAGS)
AM_LDFLAGS = 
AM_LIBS = $(OPENMP_CXXFLAGS) $(BOOST_PYTHON_LDFLAGS) $(BOOST_PYTHON_LIBS) $(PYTHON_LDFLAGS)
AM_CXXFLAGS = $(OPENMP_CXXFLAGS)

if ENABLE_AM_DEBUG
else
  AM_CXXFLAGS += -DNDEBUG
endif

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



