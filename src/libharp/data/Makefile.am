
# Initialize compile and linking flags

AM_CPPFLAGS =
AM_CXXFLAGS =
AM_LDFLAGS =
AM_LIBS = 

# Append compile and linking flags based on configure variables
include $(top_srcdir)/src/build_options.am


noinst_HEADERS = harp_data_internal.hpp

noinst_LTLIBRARIES = libharpdata.la

libharpdata_la_includedir = $(includedir)/harp

libharpdata_la_include_HEADERS = \
	harp/fits.hpp \
	harp/targets.hpp \
	harp/image.hpp \
	harp/spec.hpp \
	harp/psf.hpp \
	harp/data.hpp \
	harp/metadata.hpp

libharpdata_la_SOURCES = \
	harp_image.cpp \
	harp_spec.cpp \
	harp_psf.cpp \
	harp_fits.cpp \
	harp_targets.cpp \
	harp_plugin_registry.cpp \
	harp_data.cpp \
	harp_metadata.cpp

libharpdata_la_LIBADD = 

MAINTAINERCLEANFILES = git-version.cpp

