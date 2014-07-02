
SUBDIRS = math data extract plugins

AM_CPPFLAGS = -I$(top_srcdir)/src/libharp/math -I$(top_srcdir)/src/libharp/data -I$(top_srcdir)/src/libharp/extract -I$(top_srcdir)/src/libharp/plugins -I$(top_srcdir)/src/libharp $(CFITSIO_CPPFLAGS) $(BOOST_CPPFLAGS) $(OPENMP_CXXFLAGS)

AM_LDFLAGS =

AM_LIBS = $(OPENMP_CXXFLAGS)

AM_CXXFLAGS = $(OPENMP_CXXFLAGS)

if ENABLE_AM_DEBUG
else
  AM_CXXFLAGS += -DNDEBUG
endif

noinst_HEADERS = harp_internal.hpp
include_HEADERS = harp.hpp

noinst_SCRIPTS = gen_harpconfig.sh
bin_SCRIPTS = harpconfig


config_serial_link = -lharp $(LDFLAGS) $(BOOST_PROGRAM_OPTIONS_LDFLAGS) $(BOOST_PROGRAM_OPTIONS_LIBS) $(BOOST_REGEX_LDFLAGS) $(BOOST_REGEX_LIBS)
config_mpi_link = -lharp-mpi $(LDFLAGS) $(ELEMENTAL) $(BOOST_MPI_LDFLAGS) $(BOOST_MPI_LIBS) $(BOOST_PROGRAM_OPTIONS_LDFLAGS) $(BOOST_PROGRAM_OPTIONS_LIBS) $(BOOST_REGEX_LDFLAGS) $(BOOST_REGEX_LIBS)

if HAVE_AM_PYTHON
  config_serial_link += -lboost_numpy $(BOOST_PYTHON_LDFLAGS) $(BOOST_PYTHON_LIBS) $(BOOST_SERIALIZATION_LDFLAGS) $(BOOST_SERIALIZATION_LIBS) $(PYTHON_LDFLAGS) $(PYTHON_LIBS)
  config_mpi_link += -lboost_numpy $(BOOST_PYTHON_LDFLAGS) $(BOOST_PYTHON_LIBS) $(BOOST_SERIALIZATION_LDFLAGS) $(BOOST_SERIALIZATION_LIBS) $(PYTHON_LDFLAGS) $(PYTHON_LIBS)
else
  config_serial_link += $(BOOST_SERIALIZATION_LDFLAGS) $(BOOST_SERIALIZATION_LIBS)
  config_mpi_link += $(BOOST_SERIALIZATION_LDFLAGS) $(BOOST_SERIALIZATION_LIBS)
endif

config_serial_link += $(CFITSIO) $(LAPACK_LIBS) $(BLAS_LIBS) $(FCLIBS) $(OPENMP_CXXFLAGS) $(LIBS) -lm
config_mpi_link += $(CFITSIO) $(LAPACK_LIBS) $(BLAS_LIBS) $(FCLIBS) $(OPENMP_CXXFLAGS) $(LIBS) -lm


harpconfig : harpconfig.in
	@./gen_harpconfig.sh \
	"$(CXX)" \
	"$(CXXFLAGS) $(OPENMP_CXXFLAGS) $(HARP_PLUGIN_SHARED)" \
	"-I$(DESTDIR)$(includedir) $(CFITSIO_CPPFLAGS) $(BOOST_CPPFLAGS) $(ELEMENTAL_CPPFLAGS)" \
	"-L$(DESTDIR)$(libdir)" \
	"$(config_serial_link)" \
	"$(MPICXX)" \
	"$(config_mpi_link)" \
	"$(pyexecdir)" \
	"$(HARP_PLUGIN_COMPILE)" \
	"$(HARP_PLUGIN_LINK)" \
	"$(HARP_PLUGIN_EXT)"

lib_LTLIBRARIES = libharp.la

libharp_la_includedir = $(includedir)/harp

libharp_la_include_HEADERS =

libharp_la_SOURCES = harp.cpp

libharp_la_LDFLAGS = -version-info @HARP_LT_VERSION@

libharp_la_LIBADD = \
	$(top_builddir)/src/libharp/data/libharpdata.la \
	$(top_builddir)/src/libharp/math/libharpmath.la \
	$(top_builddir)/src/libharp/plugins/libharpplugins.la \
	$(top_builddir)/src/libharp/extract/libharpextract.la


CLEANFILES = harpconfig

EXTRA_DIST = harpconfig.in gen_harpconfig.sh

MAINTAINERCLEANFILES = 

