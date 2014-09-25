
# Initialize compile and linking flags

AM_CPPFLAGS =
AM_CXXFLAGS =
AM_LDFLAGS =
AM_LIBS = \
$(top_builddir)/src/libharp/libharp.la

# Append compile and linking flags based on configure variables
include $(top_srcdir)/src/build_options.am


bin_PROGRAMS = harp_test_python

TESTS = harp_test_python

harp_test_python_SOURCES = harp_test_python.cpp
harp_test_python_LDADD = $(AM_LIBS)

MAINTAINERCLEANFILES =

