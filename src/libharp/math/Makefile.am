
AM_CPPFLAGS = -I$(top_srcdir)/src/libharp/math $(CFITSIO_CPPFLAGS) $(BOOST_CPPFLAGS) $(OPENMP_CXXFLAGS)
AM_LDFLAGS = -L$(top_builddir)/src/libharp/math
AM_LIBS =

noinst_LTLIBRARIES = libharpmath.la

libharpmath_la_includedir = $(includedir)/harp
libharpmath_la_include_HEADERS = harp/math.hpp harp/linalg.hpp
libharpmath_la_SOURCES = harp_math.cpp
libharp_la_LIBADD = $(AM_LIBS)

