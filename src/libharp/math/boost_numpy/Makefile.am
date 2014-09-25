
# Initialize compile and linking flags

AM_CPPFLAGS =
AM_CXXFLAGS =
AM_LDFLAGS =
AM_LIBS = 

# Append compile and linking flags based on configure variables
include $(top_srcdir)/src/build_options.am


noinst_LTLIBRARIES = libboost_numpy.la

#harp_test_numpy_dtype_mod.la harp_test_numpy_indexing_mod.la harp_test_numpy_ndarray_mod.la harp_test_numpy_shapes_mod.la harp_test_numpy_templates_mod.la harp_test_numpy_ufunc_mod.la

libboost_numpy_la_includedir = $(includedir)/boost

libboost_numpy_la_include_HEADERS = \
	boost/numpy.hpp \
	boost/numpy/internal.hpp \
	boost/numpy/dtype.hpp \
	boost/numpy/invoke_matching.hpp \
	boost/numpy/matrix.hpp \
	boost/numpy/ndarray.hpp \
	boost/numpy/numpy_object_mgr_traits.hpp \
	boost/numpy/scalars.hpp \
	boost/numpy/ufunc.hpp
	
libboost_numpy_la_SOURCES = \
	dtype.cpp \
	matrix.cpp \
	ndarray.cpp \
	numpy.cpp \
	scalars.cpp \
	ufunc.cpp

libboost_numpy_la_LDFLAGS = 

libboost_numpy_la_LIBADD = 


#harp_test_numpy_dtype_mod_la_LDFLAGS = -module
#harp_test_numpy_dtype_mod_la_LIBADD = $(top_builddir)/src/libharp/math/boost_numpy/libboost_numpy.la $(AM_LIBS)
#harp_test_numpy_dtype_mod_la_SOURCES = dtype_mod.cpp

#harp_test_numpy_indexing_mod_la_LDFLAGS = -module
#harp_test_numpy_indexing_mod_la_LIBADD = $(top_builddir)/src/libharp/math/boost_numpy/libboost_numpy.la $(AM_LIBS)
#harp_test_numpy_indexing_mod_la_SOURCES = dtype_mod.cpp

#harp_test_numpy_ndarray_mod_la_LDFLAGS = -module
#harp_test_numpy_ndarray_mod_la_LIBADD = $(top_builddir)/src/libharp/math/boost_numpy/libboost_numpy.la $(AM_LIBS)
#harp_test_numpy_ndarray_mod_la_SOURCES = dtype_mod.cpp

#harp_test_numpy_shapes_mod_la_LDFLAGS = -module
#harp_test_numpy_shapes_mod_la_LIBADD = $(top_builddir)/src/libharp/math/boost_numpy/libboost_numpy.la $(AM_LIBS)
#harp_test_numpy_shapes_mod_la_SOURCES = dtype_mod.cpp

#harp_test_numpy_templates_mod_la_LDFLAGS = -module
#harp_test_numpy_templates_mod_la_LIBADD = $(top_builddir)/src/libharp/math/boost_numpy/libboost_numpy.la $(AM_LIBS)
#harp_test_numpy_templates_mod_la_SOURCES = dtype_mod.cpp

#harp_test_numpy_ufunc_mod_la_LDFLAGS = -module
#harp_test_numpy_ufunc_mod_la_LIBADD = $(top_builddir)/src/libharp/math/boost_numpy/libboost_numpy.la $(AM_LIBS)
#harp_test_numpy_ufunc_mod_la_SOURCES = dtype_mod.cpp



