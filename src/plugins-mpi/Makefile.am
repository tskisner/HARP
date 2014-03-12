CXX = @MPICXX@

SUBDIRS = 

AM_CPPFLAGS = -I$(top_srcdir)/src/plugins -I$(top_srcdir)/src/libharp/math -I$(top_srcdir)/src/libharp/data -I$(top_srcdir)/src/libharp -I$(top_srcdir)/src/libharp-mpi -I$(top_srcdir)/src/libharp-mpi/math $(ELEMENTAL_CPPFLAGS) $(CFITSIO_CPPFLAGS) $(BOOST_CPPFLAGS) $(OPENMP_CXXFLAGS)
AM_LDFLAGS =
AM_LIBS = $(OPENMP_CXXFLAGS)

AM_CXXFLAGS = $(OPENMP_CXXFLAGS)

if ENABLE_AM_DEBUG
else
  AM_CXXFLAGS += -DNDEBUG
endif
  
lib_LTLIBRARIES = libharp-mpi-plugins.la

libharp_mpi_plugins_la_SOURCES = 
libharp_mpi_plugins_la_LIBADD = 
libharp_mpi_plugins_la_CXXFLAGS = $(AM_CXXFLAGS) -DUSE_MPI

EXTRA_DIST = plugin.am

MAINTAINERCLEANFILES = plugin.am

CLEANFILES =

include plugin.am

