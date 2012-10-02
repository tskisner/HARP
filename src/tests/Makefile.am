CXX = @MPICXX@

EXTRA_DIST = testdata test_instruments.am harp_testcommands.cpp

AM_CPPFLAGS = -I$(top_srcdir)/src/libharp -I$(top_srcdir)/src/libharp/matrix -I$(top_srcdir)/src/tests $(BOOST_CPPFLAGS) $(ELEMENTAL_CPPFLAGS) $(CFITSIO_CPPFLAGS)
AM_LDFLAGS =
AM_LIBS = $(top_builddir)/src/libharp/libharp.la $(BOOST_LDFLAGS) $(ELEMENTAL) $(CFITSIO) -lm

noinst_HEADERS = harp_test.hpp

bin_PROGRAMS = harp_test
TESTS = harp_test
harp_test_SOURCES = harp_test.cpp harp_test_eigen.cpp harp_test_tinyKLT.cpp
harp_test_LDADD = $(AM_LIBS)

# Manual dependency on the data directory

harp_test.$(OBJEXT) : testdata/test_image.fits

testdata/test_image.fits :
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

include test_instruments.am

MAINTAINERCLEANFILES = test_instruments.am harp_testcommands.cpp harp_test_instruments.hpp

