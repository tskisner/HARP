
SUBDIRS = math data

AM_CPPFLAGS = -I$(top_srcdir)/src/libharp/math -I$(top_srcdir)/src/libharp/data -I$(top_srcdir)/src/libharp -I$(top_srcdir)/src/plugins $(CFITSIO_CPPFLAGS) $(BOOST_CPPFLAGS) $(OPENMP_CXXFLAGS)
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

harpconfig : harpconfig.in
	@./gen_harpconfig.sh \
	"$(CXX)" \
	"$(CXXFLAGS)" \
	"-I$(DESTDIR)$(includedir) $(CFITSIO_CPPFLAGS) $(BOOST_CPPFLAGS) $(ELEMENTAL_CPPFLAGS)" \
	"-L$(DESTDIR)$(libdir)" \
	"-lharp-plugins -lharp $(LDFLAGS) $(CFITSIO) $(BOOST_SERIALIZATION_LDFLAGS) $(BOOST_PROGRAM_OPTIONS_LDFLAGS) $(BOOST_REGEX_LDFLAGS) $(BOOST_REGEX_LIBS) $(BOOST_SERIALIZATION_LIBS) $(BOOST_PROGRAM_OPTIONS_LIBS) $(LAPACK_LIBS) $(BLAS_LIBS) $(FCLIBS) $(OPENMP_CXXFLAGS) $(LIBS) -lm" \
	"$(MPICXX)" \
	"-lharp-mpi-plugins -lharp-mpi $(ELEMENTAL) $(BOOST_MPI_LDFLAGS) $(BOOST_MPI_LIBS)" \
	"$(pyexecdir)" \
	"$(HARP_PLUGIN_COMPILE)" \
	"$(HARP_PLUGIN_LINK)"

lib_LTLIBRARIES = libharp.la

libharp_la_includedir = $(includedir)/harp

libharp_la_include_HEADERS = \
	harp/spec_slice.hpp \
	harp/extraction.hpp

libharp_la_SOURCES = \
	harp.cpp \
	harp_spec_slice.cpp \
	harp_extraction.cpp 

libharp_la_LDFLAGS = -version-info @HARP_LT_VERSION@

libharp_la_LIBADD = \
	$(top_builddir)/src/libharp/data/libharpdata.la \
	$(top_builddir)/src/libharp/math/libharpmath.la


CLEANFILES = harpconfig

EXTRA_DIST = harpconfig.in gen_harpconfig.sh

MAINTAINERCLEANFILES = 

