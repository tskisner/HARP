CXX = @MPICXX@

AM_CPPFLAGS = -I$(top_srcdir)/src/libharp -I$(top_srcdir)/src/apps $(CFITSIO_CPPFLAGS) $(BOOST_CPPFLAGS) $(CLIQUE_CPPFLAGS) $(OPENMP_CXXFLAGS)
AM_LDFLAGS =
AM_LIBS = $(top_builddir)/src/libharp/libharp.la $(CFITSIO) $(CLIQUE) $(BOOST_LDFLAGS) $(BOOST_PROGRAM_OPTIONS_LIB) $(OPENMP_CXXFLAGS) -lm

if HAVE_AM_OPENMP
	AM_LIBS += $(OPENMP_CXXFLAGS)
endif

#bin_PROGRAMS = harp_specstract

#harp_specstract_SOURCES = harp_specstract.cpp
#harp_specstract_LDADD = $(AM_LIBS)

