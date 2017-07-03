
# Initialize compile and linking flags

AM_CPPFLAGS =
AM_CXXFLAGS =
AM_LDFLAGS =
AM_LIBS =

# Append compile and linking flags based on configure variables
include $(top_srcdir)/src/build_options_mpi.am


bin_SCRIPTS =

bin_PROGRAMS = harp_mpi_extract harp_mpi_desi

BUILT_SOURCES =

CLEANFILES =

harp_mpi_extract_SOURCES = harp_mpi_extract.cpp
harp_mpi_extract_LDADD = $(top_builddir)/src/libharp-mpi/libharp-mpi.la

harp_mpi_desi_SOURCES = harp_mpi_desi.cpp
harp_mpi_desi_LDADD = $(top_builddir)/src/libharp-mpi/libharp-mpi.la

