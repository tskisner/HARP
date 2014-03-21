
EXTRA_DIST = testdata

AM_CPPFLAGS = -I$(top_srcdir)/src/libharp/math -I$(top_srcdir)/src/libharp/data -I$(top_srcdir)/src/libharp/extract -I$(top_srcdir)/src/libharp/plugins -I$(top_srcdir)/src/libharp -I$(top_srcdir)/src/tests $(BOOST_CPPFLAGS) $(CFITSIO_CPPFLAGS)
AM_LDFLAGS =
AM_LIBS = $(top_builddir)/src/libharp/libharp.la $(BOOST_SERIALIZATION_LDFLAGS) $(BOOST_PROGRAM_OPTIONS_LDFLAGS) $(BOOST_REGEX_LDFLAGS) $(BOOST_REGEX_LIBS) $(BOOST_SERIALIZATION_LIBS) $(BOOST_PROGRAM_OPTIONS_LIBS) $(CFITSIO) $(LAPACK_LIBS) $(BLAS_LIBS) $(OPENMP_CXXFLAGS) -lm
AM_CXXFLAGS = $(OPENMP_CXXFLAGS)

if ENABLE_AM_DEBUG
else
  AM_CXXFLAGS += -DNDEBUG
endif

if USE_AM_ALLSTATIC
  AM_LDFLAGS += -all-static
endif

noinst_HEADERS = harp_test.hpp harp_test_serialize.hpp

bin_PROGRAMS = harp_test

TESTS = harp_test

harp_test_SOURCES = \
	harp_test.cpp \
	harp_test_serialize.cpp \
	harp_test_linalg.cpp \
	harp_test_spec_simspecter.cpp \
	harp_test_psf_gauss_sim.cpp \
	harp_test_psf_gauss.cpp \
	harp_test_image_simfits.cpp \
	harp_test_specslice.cpp \
	harp_test_extract.cpp

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

