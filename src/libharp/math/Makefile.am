
SUBDIRS = boost_numpy

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
libharpmath_la_LIBADD = \
$(top_builddir)/src/libharp/math/boost_numpy/libboost_numpy.la

install-data-local:
	@$(mkinstalldirs) $(DESTDIR)$(includedir)/boost; \
	cp -r $(srcdir)/boost/* $(DESTDIR)$(includedir)/boost/; \
	chmod -R a=rX,u+w $(DESTDIR)$(includedir)/boost

uninstall-local:
	@rm -rf $(DESTDIR)$(includedir)/boost/numeric/bindings; \
	rm -f $(DESTDIR)$(includedir)/boost/archive/polymorphic_portable_archive.hpp; \
	rm -f $(DESTDIR)$(includedir)/boost/archive/portable_archive*; \
	rm -f $(DESTDIR)$(includedir)/boost/archive/portable_iarchive.hpp; \
	rm -f $(DESTDIR)$(includedir)/boost/archive/portable_oarchive.hpp; \
	numcontents=`ls -A $(DESTDIR)$(includedir)/boost/numeric`; \
	archcontents=`ls -A $(DESTDIR)$(includedir)/boost/archive`; \
	if [ x"$$numcontents$$archcontents" = "x" ]; then \
		rm -rf $(DESTDIR)$(includedir)/boost; \
	fi; \
	rm -rf $(DESTDIR)$(includedir)/harp

