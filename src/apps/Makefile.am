CXX = @MPICXX@

AM_CPPFLAGS = -I$(top_srcdir)/src/libharp -I$(top_srcdir)/src/apps $(CFITSIO_CPPFLAGS) $(BOOST_CPPFLAGS) $(CLIQUE_CPPFLAGS) $(OPENMP_CXXFLAGS)

AM_LDFLAGS =

AM_LIBS = $(top_builddir)/src/libharp/libharp.la $(CFITSIO) $(BOOST_LDFLAGS) $(BOOST_PROGRAM_OPTIONS_LIB) $(CLIQUE) $(OPENMP_CXXFLAGS) -lm

bin_PROGRAMS = harp_extract

harp_extract_SOURCES = harp_extract.cpp
harp_extract_LDADD = $(AM_LIBS)

