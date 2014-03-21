CXX = @MPICXX@

AM_CPPFLAGS = -I$(top_srcdir)/src/libharp-mpi -I$(top_srcdir)/src/libharp-mpi/math -I$(top_srcdir)/src/libharp/math -I$(top_srcdir)/src/libharp/data -I$(top_srcdir)/src/libharp/plugins -I$(top_srcdir)/src/libharp/extract -I$(top_srcdir)/src/libharp -I$(top_srcdir)/src/apps-mpi $(CFITSIO_CPPFLAGS) $(BOOST_CPPFLAGS) $(ELEMENTAL_CPPFLAGS) $(CLIQUE_CPPFLAGS)
AM_LDFLAGS =
AM_LIBS = $(top_builddir)/src/libharp-mpi/libharp-mpi.la $(CFITSIO) $(BOOST_SERIALIZATION_LDFLAGS) $(BOOST_PROGRAM_OPTIONS_LDFLAGS) $(BOOST_REGEX_LDFLAGS) $(BOOST_REGEX_LIBS) $(BOOST_SERIALIZATION_LIBS) $(BOOST_PROGRAM_OPTIONS_LIB) $(CLIQUE) $(ELEMENTAL) $(LAPACK_LIBS) $(BLAS_LIBS) $(FCLIBS) $(OPENMP_CXXFLAGS) -lm
AM_CXXFLAGS = $(OPENMP_CXXFLAGS)

if ENABLE_AM_DEBUG
else
  AM_CXXFLAGS += -DNDEBUG
endif

if USE_AM_ALLSTATIC
  AM_LDFLAGS += -all-static
endif

bin_SCRIPTS =
bin_PROGRAMS = harp_extract
BUILT_SOURCES =
CLEANFILES =

harp_extract_SOURCES = harp_extract.cpp
harp_extract_LDADD = $(AM_LIBS)

