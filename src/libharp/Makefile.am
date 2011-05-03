AM_CPPFLAGS = -I$(top_srcdir)/src/libharp $(CFITSIO_CPPFLAGS) $(BOOST_CPPFLAGS) $(MOAT_CPPFLAGS)
AM_LDFLAGS = -L$(top_builddir)/src/libharp
AM_LIBS = $(CFITSIO) $(BOOST_LDFLAGS) $(MOAT) -lm

include_HEADERS = harp.hpp

lib_LTLIBRARIES = libharp.la

libharp_la_includedir = $(includedir)/harp

libharp_la_include_HEADERS = harp/image.hpp harp/spectra.hpp harp/psf.hpp

libharp_la_SOURCES = harp_image.cpp harp_spectra.cpp harp_psf.cpp


include instruments.am

MAINTAINERCLEANFILES = git-version.cpp instruments.am instruments.hpp harp_image_formats.cpp harp_spectra_formats.cpp harp_psf_formats.cpp

