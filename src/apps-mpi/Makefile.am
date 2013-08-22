CXX = @MPICXX@

AM_CPPFLAGS = -I$(top_srcdir)/src/libharp-mpi -I$(top_srcdir)/src/libharp -I$(top_srcdir)/src/apps-mpi $(CFITSIO_CPPFLAGS) $(BOOST_CPPFLAGS) $(ELEMENTAL_CPPFLAGS) $(CLIQUE_CPPFLAGS) $(OPENMP_CXXFLAGS)

AM_LDFLAGS =

AM_LIBS = $(top_builddir)/src/libharp-mpi/libharp-mpi.la $(top_builddir)/src/libharp/libharp.la $(CFITSIO) $(BOOST_LDFLAGS) $(BOOST_PROGRAM_OPTIONS_LIB) $(CLIQUE) $(ELEMENTAL) $(LAPACK_LIBS) $(BLAS_LIBS) $(FCLIBS) $(OPENMP_CXXFLAGS) -lm

bin_SCRIPTS =
bin_PROGRAMS = harp_extract
BUILT_SOURCES =
CLEANFILES =

harp_extract_SOURCES = harp_extract.cpp
harp_extract_LDADD = $(AM_LIBS)

