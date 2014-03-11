
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

if HAVE_AM_MPI
  AM_CPPFLAGS += -I$(top_srcdir)/src/libharp-mpi -I$(top_srcdir)/src/libharp-mpi/math $(ELEMENTAL_CPPFLAGS)
  
  lib_LTLIBRARIES += libharp-mpi-plugins.la

  libharp_mpi_plugins_la_includedir = $(includedir)/harp
  libharp_mpi_plugins_la_include_HEADERS = 
  libharp_mpi_plugins_la_SOURCES = 
  libharp_mpi_plugins_la_LIBADD = 
  libharp_mpi_plugins_la_CXX = @MPICXX@
  libharp_mpi_plugins_la_CXXFLAGS = $(AM_CXXFLAGS) -DUSE_MPI

endif

libharp_plugins_la_includedir = $(includedir)/harp
libharp_plugins_la_include_HEADERS = harp/static_plugins.hpp
libharp_plugins_la_SOURCES = 
libharp_plugins_la_LIBADD = 

EXTRA_DIST = harp/static_plugins.hpp.in plugin.am

MAINTAINERCLEANFILES = harp/static_plugins.hpp plugin.am

CLEANFILES =

include plugin.am

