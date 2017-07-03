
# Initialize compile and linking flags

AM_CPPFLAGS =
AM_CXXFLAGS =
AM_LDFLAGS =
AM_LIBS =

# Append compile and linking flags based on configure variables
include $(top_srcdir)/src/build_options.am


EXTRA_DIST = testdata

noinst_HEADERS = harp_test.hpp harp_test_serialize.hpp

bin_PROGRAMS = harp_test

TESTS = harp_test

harp_test_SOURCES = \
	harp_test.cpp \
	harp_test_serialize.cpp \
	harp_test_linalg.cpp \
	harp_test_spec_simfits.cpp \
	harp_test_psf_gauss_sim.cpp \
	harp_test_psf_gauss.cpp \
	harp_test_image_simfits.cpp \
	harp_test_desi.cpp \
	harp_test_specslice.cpp \
	harp_test_extract.cpp

harp_test_LDADD = $(top_builddir)/src/libharp/libharp.la

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

