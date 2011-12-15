SUBDIRS = matrix

AM_CPPFLAGS = -I$(top_srcdir)/src/libharp -I$(top_srcdir)/src/libharp/matrix $(CFITSIO_CPPFLAGS) $(BOOST_CPPFLAGS) $(MOAT_CPPFLAGS)
AM_LDFLAGS = -L$(top_builddir)/src/libharp -L$(top_builddir)/src/libharp/matrix
AM_LIBS = $(CFITSIO) $(MOAT) $(BOOST_LDFLAGS) -lm

noinst_HEADERS = harp_internal.hpp
include_HEADERS = harp.hpp

lib_LTLIBRARIES = libharp.la

libharp_la_includedir = $(includedir)/harp

libharp_la_include_HEADERS = harp/fits.hpp harp/image.hpp harp/spec.hpp harp/psf.hpp harp/metadata.hpp harp/extract.hpp

libharp_la_SOURCES = harp_image.cpp harp_spec.cpp harp_psf.cpp harp_fits.cpp

libharp_la_LIBADD = $(top_builddir)/src/libharp/matrix/libharpmatrix.la $(AM_LIBS)


include instruments.am

MAINTAINERCLEANFILES = git-version.cpp instruments.am instruments.hpp harp_image_formats.cpp harp_spec_formats.cpp harp_psf_formats.cpp

