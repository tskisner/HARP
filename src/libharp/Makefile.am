AM_CPPFLAGS = -I$(top_srcdir)/src/libharp $(CFITSIO_CPPFLAGS) $(BOOST_CPPFLAGS) $(MOAT_CPPFLAGS)
AM_LDFLAGS = -L$(top_builddir)/src/libharp
AM_LIBS = $(CFITSIO) $(MOAT) $(BOOST_LDFLAGS) -lm

noinst_HEADERS = harp_internal.hpp
include_HEADERS = harp.hpp

lib_LTLIBRARIES = libharp.la

libharp_la_includedir = $(includedir)/harp

libharp_la_include_HEADERS = harp/fits.hpp harp/image.hpp harp/spectrum.hpp harp/psf.hpp

libharp_la_SOURCES = harp_image.cpp harp_spectrum.cpp harp_psf.cpp harp_fits.cpp


include instruments.am

MAINTAINERCLEANFILES = git-version.cpp instruments.am instruments.hpp harp_image_formats.cpp harp_spectrum_formats.cpp harp_psf_formats.cpp

