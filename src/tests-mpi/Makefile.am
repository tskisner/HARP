
# Initialize compile and linking flags

AM_CPPFLAGS =
AM_CXXFLAGS =
AM_LDFLAGS =
AM_LIBS =

# Append compile and linking flags based on configure variables
include $(top_srcdir)/src/build_options_mpi.am


noinst_HEADERS = harp_mpi_test.hpp

bin_PROGRAMS = harp_mpi_test

TESTS = harp_mpi_test

harp_mpi_test_SOURCES = \
	harp_mpi_test.cpp \
	harp_mpi_test_linalg.cpp \
	harp_mpi_test_specslice.cpp \
	harp_mpi_test_spec.cpp \
	harp_mpi_test_extract.cpp

harp_mpi_test_LDADD = \
	$(top_builddir)/src/libharp-mpi/libharp-mpi.la
	$(AM_LIBS)


MAINTAINERCLEANFILES =

