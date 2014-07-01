AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = $(BOOST_CPPFLAGS) $(PYTHON_CPPFLAGS) $(NUMPY_CPPFLAGS)
AM_LDFLAGS = 
AM_LIBS = $(BOOST_PYTHON_LDFLAGS) $(BOOST_PYTHON_LIBS) $(PYTHON_LDFLAGS)
AM_CXXFLAGS = 

if ENABLE_AM_DEBUG
else
  AM_CXXFLAGS += -DNDEBUG
endif

EXTRA_DIST = \
	boost/numpy/internal.hpp \
	boost/numpy/dtype.hpp \
	boost/numpy/invoke_matching.hpp \
	boost/numpy/matrix.hpp \
	boost/numpy/ndarray.hpp \
	boost/numpy/numpy_object_mgr_traits.hpp \
	boost/numpy/scalars.hpp \
	boost/numpy/ufunc.hpp

lib_LTLIBRARIES = libboost_numpy.la

libboost_numpy_la_includedir = $(includedir)/boost

libboost_numpy_la_include_HEADERS = \
	boost/numpy.hpp

libboost_numpy_la_SOURCES = \
	dtype.cpp \
	matrix.cpp \
	ndarray.cpp \
	numpy.cpp \
	scalars.cpp \
	ufunc.cpp

libboost_numpy_la_LDFLAGS = 

libboost_numpy_la_LIBADD = 


install-data-local:
	@$(mkinstalldirs) $(DESTDIR)$(includedir)/boost/numpy; \
	cp -r $(srcdir)/boost/numpy/* $(DESTDIR)$(includedir)/boost/numpy/; \
	chmod -R a=rX,u+w $(DESTDIR)$(includedir)/boost/numpy

uninstall-local:
	@rm -rf $(DESTDIR)$(includedir)/boost/numpy; \
	contents=`ls -A $(DESTDIR)$(includedir)/boost/`; \
	if [ x"$$contents" = "x" ]; then \
		rm -rf $(DESTDIR)$(includedir)/boost; \
	fi

