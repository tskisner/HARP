
# Initialize compile and linking flags

AM_CPPFLAGS =
AM_CXXFLAGS =
AM_LDFLAGS =
AM_LIBS = 

# Append compile and linking flags based on configure variables
include $(top_srcdir)/src/build_options.am


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

