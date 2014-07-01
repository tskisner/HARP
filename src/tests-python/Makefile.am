
AM_CPPFLAGS = $(BOOST_CPPFLAGS) $(PYTHON_CPPFLAGS) $(NUMPY_CPPFLAGS) -I$(top_srcdir)/src/python/boost_numpy
AM_LDFLAGS = 
AM_LIBS = $(top_builddir)/src/python/boost_numpy/libboost_numpy.la $(BOOST_PYTHON_LDFLAGS) $(BOOST_PYTHON_LIBS) $(BOOST_SERIALIZATION_LDFLAGS) $(BOOST_SERIALIZATION_LIBS) $(PYTHON_LDFLAGS) $(PYTHON_LIBS)
AM_CXXFLAGS =

if ENABLE_AM_DEBUG
else
  AM_CXXFLAGS += -DNDEBUG
endif

bin_PROGRAMS = harp_test_python

TESTS = harp_test_python

harp_test_python_SOURCES = harp_test_python.cpp
harp_test_python_LDADD = $(AM_LIBS)

MAINTAINERCLEANFILES =

