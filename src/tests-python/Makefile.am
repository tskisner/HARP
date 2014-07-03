
AM_CPPFLAGS = -I$(top_srcdir)/src/libharp/math -I$(top_srcdir)/src/libharp/data -I$(top_srcdir)/src/libharp/extract -I$(top_srcdir)/src/libharp/plugins -I$(top_srcdir)/src/python -I$(top_srcdir)/src/python/boost_numpy -I$(top_srcdir)/src/libharp -I$(top_srcdir)/src/tests-python $(BOOST_CPPFLAGS) $(CFITSIO_CPPFLAGS) $(PYTHON_CPPFLAGS) $(NUMPY_CPPFLAGS)

AM_LDFLAGS =

AM_LIBS = $(top_builddir)/src/libharp/libharp.la $(BOOST_PROGRAM_OPTIONS_LDFLAGS) $(BOOST_PROGRAM_OPTIONS_LIBS) $(BOOST_REGEX_LDFLAGS) $(BOOST_REGEX_LIBS)

if HAVE_AM_PYTHON
  AM_LIBS += $(top_builddir)/src/python/boost_numpy/libboost_numpy.la $(BOOST_PYTHON_LDFLAGS) $(BOOST_PYTHON_LIBS) $(BOOST_SERIALIZATION_LDFLAGS) $(BOOST_SERIALIZATION_LIBS) $(PYTHON_LDFLAGS) $(PYTHON_LIBS)
else
  AM_LIBS += $(BOOST_SERIALIZATION_LDFLAGS) $(BOOST_SERIALIZATION_LIBS)
endif

AM_LIBS += $(CFITSIO) $(LAPACK_LIBS) $(BLAS_LIBS) $(FCLIBS) $(OPENMP_CXXFLAGS) -lm

AM_CXXFLAGS = $(OPENMP_CXXFLAGS)

if ENABLE_AM_DEBUG
else
  AM_CXXFLAGS += -DNDEBUG
endif


bin_PROGRAMS = harp_test_python

TESTS = harp_test_python

harp_test_python_SOURCES = harp_test_python.cpp
harp_test_python_LDADD = $(AM_LIBS)

MAINTAINERCLEANFILES =

