
# Initialize compile and linking flags

AM_CPPFLAGS =
AM_CXXFLAGS =
AM_LDFLAGS =
AM_LIBS = 

# Append compile and linking flags based on configure variables
include $(top_srcdir)/src/build_options.am


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

