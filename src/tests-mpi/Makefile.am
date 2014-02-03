CXX = @MPICXX@

EXTRA_DIST = testdata

AM_CPPFLAGS = -I$(top_srcdir)/src/libharp-mpi -I$(top_srcdir)/src/libharp/math -I$(top_srcdir)/src/libharp/data -I$(top_srcdir)/src/libharp/plugin -I$(top_srcdir)/src/libharp -I$(top_srcdir)/src/tests $(BOOST_CPPFLAGS) $(CLIQUE_CPPFLAGS) $(CFITSIO_CPPFLAGS)
AM_LDFLAGS =
AM_LIBS = $(top_builddir)/src/libharp-mpi/libharp-mpi.la $(top_builddir)/src/libharp/libharp.la $(CFITSIO) $(BOOST_SERIALIZATION_LDFLAGS) $(BOOST_PROGRAM_OPTIONS_LDFLAGS) $(BOOST_REGEX_LDFLAGS) $(BOOST_REGEX_LIBS) $(BOOST_SERIALIZATION_LIBS) $(BOOST_PROGRAM_OPTIONS_LIB) $(CLIQUE) $(ELEMENTAL) $(LAPACK_LIBS) $(BLAS_LIBS) $(FCLIBS) $(OPENMP_CXXFLAGS) -lm
AM_CXXFLAGS = $(OPENMP_CXXFLAGS)

if ENABLE_AM_DEBUG
else
  AM_CXXFLAGS += -DNDEBUG
endif

if USE_AM_ALLSTATIC
  AM_LDFLAGS += -all-static
endif

noinst_HEADERS = harp_test.hpp

bin_PROGRAMS = harp_test

TESTS = harp_test

harp_test_SOURCES = \
	harp_test.cpp \
	harp_test_elemental.cpp \
	harp_test_invcov.cpp \
	harp_test_spec_specter.cpp \
	harp_test_spec_sim.cpp \
	harp_test_psf_gauss.cpp \
	harp_test_sim_extract.cpp \
	harp_test_small_extract.cpp \
	harp_test_sky_extract.cpp

harp_test_LDADD = $(AM_LIBS)

# Manual dependency on the data directory

harp_test.$(OBJEXT) : testdata/realtoy_input_spectra.fits

testdata/realtoy_input_spectra.fits :
	@BUNZIP=`which bunzip2`; if [ x"$${BUNZIP}" = x ]; then \
	echo "bunzip2 not found."; else \
	for bzfile in `ls testdata/*.bz2`; do \
	file=`echo $${bzfile} | sed -e 's/\(.*\)\.bz2/\1/'`; \
	cat $${bzfile} | $${BUNZIP} > $${file}; \
	done; \
	fi

clean-local :
	@rm -f testdata/*.out; \
	for bzfile in `ls testdata/*.bz2`; do \
	file=`echo $${bzfile} | sed -e 's/\(.*\)\.bz2/\1/'`; \
	rm -f $${file}; rm -f testdata/*.out; \
	done

MAINTAINERCLEANFILES =

