
SUBDIRS = 

AM_CPPFLAGS = -I$(top_srcdir)/src/libharp/math -I$(top_srcdir)/src/libharp/data -I$(top_srcdir)/src/libharp/plugin $(CFITSIO_CPPFLAGS) $(BOOST_CPPFLAGS) $(OPENMP_CXXFLAGS)
AM_LDFLAGS =
AM_LIBS = $(OPENMP_CXXFLAGS)

AM_CXXFLAGS = $(OPENMP_CXXFLAGS)

if ENABLE_AM_DEBUG
else
  AM_CXXFLAGS += -DNDEBUG
endif

noinst_HEADERS = 

noinst_LTLIBRARIES = libharpplugin.la

libharpplugin_la_includedir = $(includedir)/harp
libharpplugin_la_include_HEADERS = harp/plugin.hpp
libharpplugin_la_SOURCES = 
libharpplugin_la_LIBADD = 

EXTRA_DIST = harp/plugin.hpp.in plugin.am

MAINTAINERCLEANFILES = harp/plugin.hpp plugin.am

include plugin.am

