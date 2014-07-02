CXX = @MPICXX@

AM_CPPFLAGS = -I$(top_srcdir)/src/libharp-mpi -I$(top_srcdir)/src/libharp-mpi/math -I$(top_srcdir)/src/libharp/math -I$(top_srcdir)/src/libharp/data -I$(top_srcdir)/src/libharp/extract -I$(top_srcdir)/src/libharp/plugins -I$(top_srcdir)/src/libharp -I$(top_srcdir)/src/apps-mpi $(BOOST_CPPFLAGS) $(CFITSIO_CPPFLAGS) $(ELEMENTAL_CPPFLAGS)

AM_LDFLAGS =

AM_LIBS = $(top_builddir)/src/libharp-mpi/libharp-mpi.la $(ELEMENTAL) $(BOOST_MPI_LDFLAGS) $(BOOST_MPI_LIBS) $(BOOST_PROGRAM_OPTIONS_LDFLAGS) $(BOOST_PROGRAM_OPTIONS_LIBS) $(BOOST_REGEX_LDFLAGS) $(BOOST_REGEX_LIBS)

if HAVE_AM_PYTHON
  AM_LIBS += $(top_builddir)/src/python/boost_numpy/libboost_numpy.la $(BOOST_PYTHON_LDFLAGS) $(BOOST_PYTHON_LIBS) $(BOOST_SERIALIZATION_LDFLAGS) $(BOOST_SERIALIZATION_LIBS) $(PYTHON_LDFLAGS) $(PYTHON_LIBS)
else
  AM_LIBS += $(BOOST_SERIALIZATION_LDFLAGS) $(BOOST_SERIALIZATION_LIBS)
endif

AM_LIBS += $(CFITSIO) $(LAPACK_LIBS) $(BLAS_LIBS) $(FCLIBS) $(OPENMP_CXXFLAGS) -lm

AM_CXXFLAGS = $(OPENMP_CXXFLAGS)

if ENABLE_AM_DEBUG
else
  AM_CXXFLAGS += -DNDEBUG
endif

if USE_AM_ALLSTATIC
  AM_LDFLAGS += -all-static
endif

bin_SCRIPTS =
bin_PROGRAMS = harp_mpi_extract
BUILT_SOURCES =
CLEANFILES =

harp_mpi_extract_SOURCES = harp_mpi_extract.cpp
harp_mpi_extract_LDADD = $(AM_LIBS)

