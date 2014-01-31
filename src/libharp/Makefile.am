
SUBDIRS = math data plugin

AM_CPPFLAGS = -I$(top_srcdir)/src/libharp/math -I$(top_srcdir)/src/libharp/data -I$(top_srcdir)/src/libharp/plugin -I$(top_srcdir)/src/libharp $(CFITSIO_CPPFLAGS) $(BOOST_CPPFLAGS) $(OPENMP_CXXFLAGS)
AM_LDFLAGS =
AM_LIBS = $(OPENMP_CXXFLAGS)

AM_CXXFLAGS = $(OPENMP_CXXFLAGS)

if ENABLE_AM_DEBUG
else
  AM_CXXFLAGS += -DNDEBUG
endif

noinst_HEADERS = harp_internal.hpp harp_serialization.hpp
include_HEADERS = harp.hpp

noinst_SCRIPTS = gen_harpconfig.sh
bin_SCRIPTS = harpconfig

harpconfig : harpconfig.in
	@./gen_harpconfig.sh \
	"$(CXX)" \
	"$(CXXFLAGS)" \
	"-I$(DESTDIR)$(includedir) $(CFITSIO_CPPFLAGS) $(BOOST_CPPFLAGS)" \
	"-L$(DESTDIR)$(libdir)" \
	"-lharp $(LDFLAGS) $(CFITSIO) $(BOOST_SERIALIZATION_LDFLAGS) $(BOOST_PROGRAM_OPTIONS_LDFLAGS) $(BOOST_SERIALIZATION_LIBS) $(BOOST_PROGRAM_OPTIONS_LIBS) $(LAPACK_LIBS) $(BLAS_LIBS) $(FCLIBS) $(OPENMP_CXXFLAGS) $(LIBS) -lm" \
	"$(MPICXX)" \
	"-lharp-mpi" \
	"$(pyexecdir)" \
	"$(HARP_PLUGIN_COMPILE)" \
	"$(HARP_PLUGIN_LINK)"

lib_LTLIBRARIES = libharp.la

libharp_la_includedir = $(includedir)/harp

libharp_la_include_HEADERS = \
	harp/fits.hpp \
	harp/target.hpp \
	harp/image.hpp \
	harp/spec.hpp \
	harp/psf.hpp \
	harp/spec_slice.hpp \
	harp/extraction.hpp

libharp_la_SOURCES = \
	harp.cpp \
	harp_image.cpp \
	harp_spec.cpp \
	harp_psf.cpp \
	harp_fits.cpp \
	harp_spec_slice.cpp \
	harp_extraction.cpp 
#harp_sky.cpp

libharp_la_LDFLAGS = -version-info @HARP_LT_VERSION@

libharp_la_LIBADD = \
	$(top_builddir)/src/libharp/plugin/libharpplugin.la \
	$(top_builddir)/src/libharp/data/libharpdata.la \
	$(top_builddir)/src/libharp/math/libharpmath.la

# Add new formats here

libharp_la_include_HEADERS += \
	harp/formats_spec.hpp \
	harp/formats_image.hpp \
	harp/formats_psf.hpp

libharp_la_SOURCES += \
	harp_spec_specter.cpp \
	harp_spec_sim.cpp \
	harp_image_fits.cpp \
	harp_image_sim.cpp \
	harp_psf_gauss.cpp \
	harp_psf_gauss_sim.cpp \
	harp_psf_gh.cpp

CLEANFILES = harpconfig

EXTRA_DIST = harpconfig.in gen_harpconfig.sh

MAINTAINERCLEANFILES = 

