
SUBDIRS = 

AM_CPPFLAGS = -I$(top_srcdir)/src/libharp/plugins -I$(top_srcdir)/src/libharp/math -I$(top_srcdir)/src/libharp/data -I$(top_srcdir)/src/libharp/extract -I$(top_srcdir)/src/libharp $(CFITSIO_CPPFLAGS) $(BOOST_CPPFLAGS) $(OPENMP_CXXFLAGS)
AM_LDFLAGS =
AM_LIBS = $(OPENMP_CXXFLAGS)

AM_CXXFLAGS = $(OPENMP_CXXFLAGS)

if ENABLE_AM_DEBUG
else
  AM_CXXFLAGS += -DNDEBUG
endif

noinst_HEADERS = 

noinst_LTLIBRARIES = libharpplugins.la

libharpplugins_la_includedir = $(includedir)/harp
libharpplugins_la_include_HEADERS = harp/static_plugins.hpp
libharpplugins_la_SOURCES = 
libharpplugins_la_LIBADD = 

EXTRA_DIST = harp/static_plugins.hpp.in plugin.am

MAINTAINERCLEANFILES = harp/static_plugins.hpp plugin.am plugin_register.cpp

CLEANFILES =

include plugin.am

