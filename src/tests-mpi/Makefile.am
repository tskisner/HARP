CXX = @MPICXX@

EXTRA_DIST = testdata

AM_CPPFLAGS = -I$(top_srcdir)/src/libharp -I$(top_srcdir)/src/libharp/formats -I$(top_srcdir)/src/tests $(BOOST_CPPFLAGS) $(CLIQUE_CPPFLAGS) $(CFITSIO_CPPFLAGS) $(OPENMP_CXXFLAGS)
AM_LDFLAGS =
AM_LIBS = $(top_builddir)/src/libharp/libharp.la $(BOOST_LDFLAGS) $(CFITSIO) $(CLIQUE) $(OPENMP_CXXFLAGS) -lm

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

