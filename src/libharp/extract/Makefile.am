
SUBDIRS = 

AM_CPPFLAGS = -I$(top_srcdir)/src/libharp/math -I$(top_srcdir)/src/libharp/data -I$(top_srcdir)/src/libharp/plugins -I$(top_srcdir)/src/libharp/extract -I$(top_srcdir)/src/libharp $(CFITSIO_CPPFLAGS) $(BOOST_CPPFLAGS) $(OPENMP_CXXFLAGS)
AM_LDFLAGS =
AM_LIBS = $(OPENMP_CXXFLAGS)

AM_CXXFLAGS = $(OPENMP_CXXFLAGS)

if ENABLE_AM_DEBUG
else
  AM_CXXFLAGS += -DNDEBUG
endif

noinst_LTLIBRARIES = libharpextract.la

libharpextract_la_includedir = $(includedir)/harp

libharpextract_la_include_HEADERS = \
	harp/extraction.hpp \
	harp/spec_slice.hpp

libharpextract_la_SOURCES = \
	harp_extraction.cpp \
	harp_spec_slice.cpp

libharpextract_la_LIBADD = 

MAINTAINERCLEANFILES = 

