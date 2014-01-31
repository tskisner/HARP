CXX = @MPICXX@

SUBDIRS = math

AM_CPPFLAGS = -I$(top_srcdir)/src/libharp/math -I$(top_srcdir)/src/libharp $(CFITSIO_CPPFLAGS) $(BOOST_CPPFLAGS) -I$(top_srcdir)/src/libharp-mpi/math -I$(top_srcdir)/src/libharp-mpi $(ELEMENTAL_CPPFLAGS)
AM_LDFLAGS =
AM_LIBS = $(OPENMP_CXXFLAGS)

AM_CXXFLAGS = $(OPENMP_CXXFLAGS)

if ENABLE_AM_DEBUG
else
  AM_CXXFLAGS += -DNDEBUG
endif

noinst_HEADERS = harp_mpi_internal.hpp
include_HEADERS = harp_mpi.hpp

lib_LTLIBRARIES = libharp-mpi.la

libharp_mpi_la_includedir = $(includedir)/harp

libharp_mpi_la_include_HEADERS = \
	harp/mpi_fits.hpp \
	harp/mpi_image.hpp \
	harp/mpi_spec.hpp \
	harp/mpi_psf.hpp \
	harp/mpi_spec_slice.hpp \
	harp/mpi_extraction.hpp

libharp_mpi_la_SOURCES = \
	harp_mpi.cpp \
	harp_mpi_image.cpp \
	harp_mpi_spec.cpp \
	harp_mpi_psf.cpp
#
#	harp_mpi_fits.cpp \
#	harp_mpi_spec_slice.cpp
# harp_extraction.cpp harp_sky.cpp

libharp_mpi_la_LDFLAGS = -version-info @HARP_LT_VERSION@

libharp_mpi_la_LIBADD = $(top_builddir)/src/libharp-mpi/math/libharpmath-mpi.la

