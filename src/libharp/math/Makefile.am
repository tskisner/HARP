
AM_CPPFLAGS = -I$(top_srcdir)/src/libharp/math $(CFITSIO_CPPFLAGS) $(BOOST_CPPFLAGS) $(OPENMP_CXXFLAGS)
AM_LDFLAGS = 
AM_LIBS =

EXTRA_DIST = boost

noinst_HEADERS = harp_math_internal.hpp

noinst_LTLIBRARIES = libharpmath.la

libharpmath_la_includedir = $(includedir)/harp
libharpmath_la_include_HEADERS = harp/math.hpp harp/linalg.hpp
libharpmath_la_SOURCES = harp_linalg.cpp harp_math.cpp
libharpmath_la_LIBADD = $(AM_LIBS)

install-data-local:
	@$(mkinstalldirs) $(DESTDIR)$(includedir)/boost; \
	cp -r $(srcdir)/boost/* $(DESTDIR)$(includedir)/boost/; \
	chmod -R a=rX,u+w $(DESTDIR)$(includedir)/boost

uninstall-local:
	@rm -rf $(DESTDIR)$(includedir)/boost; \
	rm -rf $(DESTDIR)$(includedir)/harp



    
