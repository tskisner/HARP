
SUBDIRS = math

AM_CPPFLAGS = -I$(top_srcdir)/src/libharp/math -I$(top_srcdir)/src/libharp/data $(CFITSIO_CPPFLAGS) $(BOOST_CPPFLAGS) $(OPENMP_CXXFLAGS)
AM_LDFLAGS =
AM_LIBS = $(OPENMP_CXXFLAGS)

AM_CXXFLAGS = $(OPENMP_CXXFLAGS)

if ENABLE_AM_DEBUG
else
  AM_CXXFLAGS += -DNDEBUG
endif

noinst_HEADERS = harp_data_internal.hpp

noinst_LTLIBRARIES = libharpdata.la

libharpdata_la_includedir = $(includedir)/harp

libharpdata_la_include_HEADERS = \
	harp/fits.hpp \
	harp/target.hpp \
	harp/image.hpp \
	harp/spec.hpp \
	harp/psf.hpp \
	harp/data.hpp

libharpdata_la_SOURCES = \
	harp_image.cpp \
	harp_spec.cpp \
	harp_psf.cpp \
	harp_fits.cpp \
	harp_plugin_registry.cpp \
	harp_data.cpp

libharp_la_LIBADD = $(top_builddir)/src/libharp/math/libharpmath.la

MAINTAINERCLEANFILES = git-version.cpp

