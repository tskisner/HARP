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

pkgharp_PYTHON = __init__.py
pkgharpdir = $(pyexecdir)/harp




