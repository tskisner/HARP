
AM_CPPFLAGS = -I$(top_srcdir)/src/libharp -I$(top_srcdir)/src/apps $(MOAT_CPPFLAGS) $(CFITSIO_CPPFLAGS) $(BOOST_CPPFLAGS)
AM_LDFLAGS =
AM_LIBS = $(top_builddir)/src/libharp/libharp.la $(CFITSIO) $(MOAT) $(BOOST_PROGRAM_OPTIONS_LIB) $(BOOST_LDFLAGS) -lm

bin_PROGRAMS = harp_specstract

harp_specstract_SOURCES = harp_specstract.cpp
harp_specstract_LDADD = $(AM_LIBS)

