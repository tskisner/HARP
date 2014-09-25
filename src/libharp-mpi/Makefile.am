
SUBDIRS = math

# Initialize compile and linking flags

AM_CPPFLAGS =
AM_CXXFLAGS =
AM_LDFLAGS =
AM_LIBS = 

# Append compile and linking flags based on configure variables
include $(top_srcdir)/src/build_options_mpi.am


noinst_HEADERS = harp_mpi_internal.hpp
include_HEADERS = harp_mpi.hpp

lib_LTLIBRARIES = libharp-mpi.la

libharp_mpi_la_includedir = $(includedir)/harp

libharp_mpi_la_include_HEADERS = \
	harp/mpi_fits.hpp \
	harp/mpi_image.hpp \
	harp/mpi_spec.hpp \
	harp/mpi_psf.hpp \
	harp/mpi_targets.hpp \
	harp/mpi_spec_slice.hpp \
	harp/mpi_extraction.hpp

libharp_mpi_la_SOURCES = \
	harp_mpi.cpp \
	harp_mpi_image.cpp \
	harp_mpi_spec.cpp \
	harp_mpi_psf.cpp \
	harp_mpi_targets.cpp \
	harp_mpi_spec_slice.cpp \
	harp_mpi_extraction.cpp

libharp_mpi_la_LDFLAGS = -version-info @HARP_LT_VERSION@

libharp_mpi_la_LIBADD = \
	$(top_builddir)/src/libharp/data/libharpdata.la \
	$(top_builddir)/src/libharp/math/libharpmath.la \
	$(top_builddir)/src/libharp/extract/libharpextract.la \
	$(top_builddir)/src/libharp/plugins/libharpplugins.la \
	$(top_builddir)/src/libharp-mpi/math/libharpmath-mpi.la

