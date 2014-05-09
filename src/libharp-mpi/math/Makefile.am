CXX = @MPICXX@

AM_CPPFLAGS = -I$(top_srcdir)/src/libharp/math $(CFITSIO_CPPFLAGS) $(BOOST_CPPFLAGS) -I$(top_srcdir)/src/libharp-mpi/math $(ELEMENTAL_CPPFLAGS)
AM_LDFLAGS = 
AM_LIBS = $(OPENMP_CXXFLAGS)

AM_CXXFLAGS = $(OPENMP_CXXFLAGS)

if ENABLE_AM_DEBUG
else
  AM_CXXFLAGS += -DNDEBUG
endif

EXTRA_DIST =

noinst_HEADERS = harp_mpi_math_internal.hpp

noinst_LTLIBRARIES = libharpmath-mpi.la

libharpmath_mpi_la_includedir = $(includedir)/harp
libharpmath_mpi_la_include_HEADERS = harp/mpi_math.hpp harp/mpi_linalg.hpp harp/mpi_comm.hpp
libharpmath_mpi_la_SOURCES = harp_mpi_linalg.cpp harp_mpi_math.cpp
libharpmath_mpi_la_LIBADD = $(AM_LIBS)
    
