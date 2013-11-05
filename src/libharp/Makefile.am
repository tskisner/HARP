
SUBDIRS = math

AM_CPPFLAGS = -I$(top_srcdir)/src/libharp/math -I$(top_srcdir)/src/libharp $(CFITSIO_CPPFLAGS) $(BOOST_CPPFLAGS) $(OPENMP_CXXFLAGS)
AM_LDFLAGS =
AM_LIBS =

noinst_HEADERS = harp_internal.hpp harp_serialization.hpp
include_HEADERS = harp.hpp

lib_LTLIBRARIES = libharp.la

libharp_la_includedir = $(includedir)/harp

libharp_la_include_HEADERS = harp/fits.hpp harp/image.hpp harp/spec.hpp harp/psf.hpp harp/extraction.hpp harp/sky.hpp

libharp_la_SOURCES = harp.cpp harp_image.cpp harp_spec.cpp harp_psf.cpp harp_fits.cpp # harp_extraction.cpp harp_sky.cpp

libharp_la_LDFLAGS = -version-info @HARP_LT_VERSION@

libharp_la_LIBADD = $(top_builddir)/src/libharp/math/libharpmath.la

# Add new formats here

libharp_la_include_HEADERS += \
  harp/formats_spec.hpp \
  harp/formats_image.hpp \
  harp/formats_psf.hpp

libharp_la_SOURCES += \
  harp_spec_specter.cpp \
  harp_spec_sim.cpp \
  harp_image_fits.cpp \
  harp_image_sim.cpp \
  harp_psf_gauss.cpp \
  harp_psf_gauss_sim.cpp \
  harp_psf_gh.cpp

MAINTAINERCLEANFILES = git-version.cpp

