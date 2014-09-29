
SUBDIRS = math data extract plugins

# Initialize compile and linking flags

AM_CPPFLAGS =
AM_CXXFLAGS =
AM_LDFLAGS =
AM_LIBS = 

# Append compile and linking flags based on configure variables
include $(top_srcdir)/src/build_options.am

# Get config variables
include $(top_srcdir)/src/config.am


noinst_HEADERS = harp_internal.hpp
include_HEADERS = harp.hpp

noinst_SCRIPTS = gen_harpconfig.sh
bin_SCRIPTS = harpconfig


harpconfig : harpconfig.in
	@./gen_harpconfig.sh \
	"$(CONF_CXX)" \
	"$(CONF_CXXFLAGS)" \
	"$(CONF_CPPFLAGS)" \
	"$(CONF_LDFLAGS)" \
	"$(CONF_LIBS)" \
	"$(CONF_MPI_CXX)" \
	"$(CONF_MPI_CXXFLAGS)" \
	"$(CONF_MPI_CPPFLAGS)" \
	"$(CONF_MPI_LDFLAGS)" \
	"$(CONF_MPI_LIBS)" \
	"$(pyexecdir)" \
	"$(HARP_PLUGIN_COMPILE)" \
	"$(HARP_PLUGIN_LINK)" \
	"$(HARP_PLUGIN_EXT)"


lib_LTLIBRARIES = libharp.la

libharp_la_includedir = $(includedir)/harp

libharp_la_include_HEADERS =

libharp_la_SOURCES = harp.cpp

libharp_la_LDFLAGS = -version-info @HARP_LT_VERSION@

libharp_la_LIBADD = \
	$(top_builddir)/src/libharp/data/libharpdata.la \
	$(top_builddir)/src/libharp/math/libharpmath.la \
	$(top_builddir)/src/libharp/plugins/libharpplugins.la \
	$(top_builddir)/src/libharp/extract/libharpextract.la


CLEANFILES = harpconfig

EXTRA_DIST = harpconfig.in gen_harpconfig.sh

MAINTAINERCLEANFILES = 

