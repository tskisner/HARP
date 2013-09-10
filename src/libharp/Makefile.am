CXX = @MPICXX@

AM_CPPFLAGS = -I$(top_srcdir)/src/libharp -I$(top_srcdir)/src/libharp/formats $(CFITSIO_CPPFLAGS) $(BOOST_CPPFLAGS) $(OPENMP_CXXFLAGS)
AM_LDFLAGS =
AM_LIBS =

noinst_HEADERS = harp_internal.hpp
include_HEADERS = harp.hpp

lib_LTLIBRARIES = libharp.la

libharp_la_includedir = $(includedir)/harp

libharp_la_include_HEADERS = harp/fits.hpp harp/matrix.hpp harp/image.hpp harp/spec.hpp harp/psf.hpp # harp/extraction.hpp harp/sky.hpp

libharp_la_SOURCES = harp.cpp harp_matrix.cpp harp_image.cpp harp_spec.cpp harp_psf.cpp harp_fits.cpp # harp_extraction.cpp harp_sky.cpp

libharp_la_LDFLAGS = -version-info @HARP_LT_VERSION@

libharp_la_LIBADD = $(AM_LIBS)

#include formats/formats.am

MAINTAINERCLEANFILES = git-version.cpp

