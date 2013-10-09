
EXTRA_DIST = testdata

AM_CPPFLAGS = -I$(top_srcdir)/src/libharp -I$(top_srcdir)/src/libharp/formats -I$(top_srcdir)/src/tests $(BOOST_CPPFLAGS) $(CLIQUE_CPPFLAGS) $(CFITSIO_CPPFLAGS) $(OPENMP_CXXFLAGS)
AM_LDFLAGS =
AM_LIBS = $(top_builddir)/src/libharp/libharp.la $(BOOST_SERIALIZATION_LDFLAGS) $(BOOST_PROGRAM_OPTIONS_LDFLAGS) $(BOOST_SERIALIZATION_LIBS) $(BOOST_PROGRAM_OPTIONS_LIB) $(CFITSIO) $(CLIQUE) $(OPENMP_CXXFLAGS) -lm

noinst_HEADERS = harp_test.hpp harp_test_serialize.hpp

bin_PROGRAMS = harp_test

TESTS = harp_test

harp_test_SOURCES = \
	harp_test.cpp \
	harp_test_serialize.cpp \
	harp_test_spec_specter.cpp \
	harp_test_spec_sim.cpp \
	harp_test_psf_gauss_sim.cpp \
	harp_test_psf_gauss.cpp \
	harp_test_image_sim.cpp \
	harp_test_image_fits.cpp

harp_test_LDADD = $(AM_LIBS)

# Manual dependency on the data directory

harp_test.$(OBJEXT) : testdata/test.txt

testdata/test.txt :
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

