AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = -I$(top_srcdir)/src/libharp -I$(top_srcdir)/src/libharp/math $(CFITSIO_CPPFLAGS) $(BOOST_CPPFLAGS) $(PYTHON_CPPFLAGS) $(NUMPY_CPPFLAGS)
AM_LIBS = $(BOOST_PYTHON_LDFLAGS) $(BOOST_PYTHON_LIBS) $(PYTHON_LDFLAGS)
AM_CXXFLAGS =

if HAVE_AM_OPENMP
  AM_CXXFLAGS += $(OPENMP_CXXFLAGS)
  AM_LIBS += $(OPENMP_CXXFLAGS)
endif

EXTRA_DIST = gen_python.sh

pkgharp_PYTHON = __init__.py
pkgharpdir = $(pyexecdir)/harp




