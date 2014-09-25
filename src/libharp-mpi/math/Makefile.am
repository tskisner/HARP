
# Initialize compile and linking flags

AM_CPPFLAGS =
AM_CXXFLAGS =
AM_LDFLAGS = 
AM_LIBS = 

# Append compile and linking flags based on configure variables
include $(top_srcdir)/src/build_options_mpi.am


EXTRA_DIST =

noinst_HEADERS = harp_mpi_math_internal.hpp

noinst_LTLIBRARIES = libharpmath-mpi.la

libharpmath_mpi_la_includedir = $(includedir)/harp
libharpmath_mpi_la_include_HEADERS = harp/mpi_math.hpp harp/mpi_linalg.hpp harp/mpi_comm.hpp
libharpmath_mpi_la_SOURCES = harp_mpi_linalg.cpp harp_mpi_math.cpp
libharpmath_mpi_la_LIBADD = 
    
