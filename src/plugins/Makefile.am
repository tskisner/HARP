
SUBDIRS = 

AM_CPPFLAGS = -I$(top_srcdir)/src/plugins -I$(top_srcdir)/src/libharp/math -I$(top_srcdir)/src/libharp/data -I$(top_srcdir)/src/libharp $(CFITSIO_CPPFLAGS) $(BOOST_CPPFLAGS) $(OPENMP_CXXFLAGS)
AM_LDFLAGS =
AM_LIBS = $(OPENMP_CXXFLAGS)

AM_CXXFLAGS = $(OPENMP_CXXFLAGS)

if ENABLE_AM_DEBUG
else
  AM_CXXFLAGS += -DNDEBUG
endif

noinst_HEADERS = 

lib_LTLIBRARIES = libharp-plugins.la

libharp_plugins_la_includedir = $(includedir)/harp
libharp_plugins_la_include_HEADERS = harp/static_plugins.hpp
libharp_plugins_la_SOURCES = 
libharp_plugins_la_LIBADD = 

EXTRA_DIST = harp/static_plugins.hpp.in plugin.am

MAINTAINERCLEANFILES = harp/static_plugins.hpp plugin.am plugin_register.cpp

CLEANFILES =

include plugin.am

