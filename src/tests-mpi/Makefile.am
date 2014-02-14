CXX = @MPICXX@

AM_CPPFLAGS = -I$(top_srcdir)/src/libharp-mpi -I$(top_srcdir)/src/libharp-mpi/math -I$(top_srcdir)/src/libharp/math -I$(top_srcdir)/src/libharp/data -I$(top_srcdir)/src/libharp/plugin -I$(top_srcdir)/src/libharp -I$(top_srcdir)/src/tests-mpi $(BOOST_CPPFLAGS) $(CFITSIO_CPPFLAGS) $(ELEMENTAL_CPPFLAGS)
AM_LDFLAGS =
AM_LIBS = $(top_builddir)/src/libharp-mpi/libharp-mpi.la $(top_builddir)/src/libharp/libharp.la $(ELEMENTAL) $(BOOST_MPI_LDFLAGS) $(BOOST_MPI_LIBS) $(BOOST_SERIALIZATION_LDFLAGS) $(BOOST_PROGRAM_OPTIONS_LDFLAGS) $(BOOST_REGEX_LDFLAGS) $(BOOST_REGEX_LIBS) $(BOOST_SERIALIZATION_LIBS) $(BOOST_PROGRAM_OPTIONS_LIBS) $(CFITSIO) $(LAPACK_LIBS) $(BLAS_LIBS) $(OPENMP_CXXFLAGS) -lm
AM_CXXFLAGS = $(OPENMP_CXXFLAGS)

if ENABLE_AM_DEBUG
else
  AM_CXXFLAGS += -DNDEBUG
endif

if USE_AM_ALLSTATIC
  AM_LDFLAGS += -all-static
endif


noinst_HEADERS = harp_mpi_test.hpp

bin_PROGRAMS = harp_mpi_test

TESTS = harp_mpi_test

harp_mpi_test_SOURCES = \
	harp_mpi_test.cpp \
	harp_mpi_test_linalg.cpp \
	harp_mpi_test_specslice.cpp

harp_mpi_test_LDADD = $(AM_LIBS)


MAINTAINERCLEANFILES =

