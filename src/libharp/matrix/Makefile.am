
AM_CPPFLAGS = -I$(top_srcdir)/src/libharp/matrix $(BOOST_CPPFLAGS)
AM_LDFLAGS = -L$(top_builddir)/src/libharp/matrix
AM_LIBS = $(BOOST_LDFLAGS) -lm


lib_LTLIBRARIES = libharpmatrix.la

libharpmatrix_la_includedir = $(includedir)/harp
libharpmatrix_la_noinstdir = $(includedir)/harp

libharpmatrix_la_noinst_HEADERS = harp/CholeskyDecomposition.hpp harp/EigenvalueDecomposition.hpp harp/LUDecomposition.hpp harp/QRDecomposition.hpp harp/SingularValueDecomposition.hpp

libharpmatrix_la_include_HEADERS = harp/matrix.hpp harp/lanczos.hpp

libharpmatrix_la_SOURCES = CholeskyDecomposition.cpp EigenvalueDecomposition.cpp LUDecomposition.cpp QRDecomposition.cpp

