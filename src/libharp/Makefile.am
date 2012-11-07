CXX = @MPICXX@

AM_CPPFLAGS = -I$(top_srcdir)/src/libharp $(CFITSIO_CPPFLAGS) $(BOOST_CPPFLAGS) $(CLIQUE_CPPFLAGS) $(OPENMP_CXXFLAGS)
AM_LDFLAGS =
AM_LIBS = $(CFITSIO) $(BOOST_LDFLAGS) $(CLIQUE) $(OPENMP_CXXFLAGS) -lm

noinst_HEADERS = harp_internal.hpp
include_HEADERS = harp.hpp

lib_LTLIBRARIES = libharp.la

libharp_la_includedir = $(includedir)/harp

libharp_la_include_HEADERS = harp/fits.hpp harp/image.hpp harp/spec.hpp harp/psf.hpp harp/extract.hpp

libharp_la_SOURCES = harp.cpp harp_image.cpp harp_spec.cpp harp_psf.cpp harp_fits.cpp harp_extract.cpp


include instruments.am

MAINTAINERCLEANFILES = git-version.cpp instruments.am instruments.hpp harp_image_formats.cpp harp_spec_formats.cpp harp_psf_formats.cpp

