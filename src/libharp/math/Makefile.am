
SUBDIRS =

if HAVE_AM_PYTHON
  SUBDIRS += boost_numpy
endif

# Initialize compile and linking flags

AM_CPPFLAGS =
AM_CXXFLAGS =
AM_LDFLAGS =
AM_LIBS = 

# Append compile and linking flags based on configure variables
include $(top_srcdir)/src/build_options.am


EXTRA_DIST = boost

noinst_HEADERS = harp_math_internal.hpp

noinst_LTLIBRARIES = libharpmath.la

libharpmath_la_includedir = $(includedir)/harp
libharpmath_la_include_HEADERS = harp/math.hpp harp/linalg.hpp
libharpmath_la_SOURCES = harp_linalg.cpp harp_math.cpp
libharpmath_la_LIBADD = 

if HAVE_AM_PYTHON
  libharpmath_la_LIBADD += $(top_builddir)/src/libharp/math/boost_numpy/libboost_numpy.la
endif

install-data-local:
	@$(mkinstalldirs) $(DESTDIR)$(includedir)/boost; \
	cp -r $(srcdir)/boost/* $(DESTDIR)$(includedir)/boost/; \
	chmod -R a=rX,u+w $(DESTDIR)$(includedir)/boost

uninstall-local:
	@rm -rf $(DESTDIR)$(includedir)/boost/numeric/bindings; \
	numcontents=`ls -A $(DESTDIR)$(includedir)/boost/numeric`; \
	if [ x"$$numcontents" = "x" ]; then \
		rm -rf $(DESTDIR)$(includedir)/boost; \
	fi; \
	rm -rf $(DESTDIR)$(includedir)/harp

