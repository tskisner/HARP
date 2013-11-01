CXX = @MPICXX@

AM_CPPFLAGS = -I$(top_srcdir)/src/libharp/math -I$(top_srcdir)/src/libharp -I$(top_srcdir)/src/apps $(CFITSIO_CPPFLAGS) $(BOOST_CPPFLAGS) $(OPENMP_CXXFLAGS)

AM_LDFLAGS =

AM_LIBS = $(top_builddir)/src/libharp/libharp.la $(CFITSIO) $(BOOST_SERIALIZATION_LDFLAGS) $(BOOST_PROGRAM_OPTIONS_LDFLAGS) $(BOOST_SERIALIZATION_LIBS) $(BOOST_PROGRAM_OPTIONS_LIB) $(LAPACK_LIBS) $(BLAS_LIBS) $(FCLIBS) $(OPENMP_CXXFLAGS) -lm

bin_SCRIPTS =
BUILT_SOURCES =
CLEANFILES =

if HAVE_AM_PYTHON
	bin_SCRIPTS += harp_slice_psf.py harp_slice_spec.py
	BUILT_SOURCES += harp_slice_psf.py harp_slice_spec.py
	CLEANFILES += harp_slice_psf.py harp_slice_spec.py
endif

harp_slice_psf.py : harp_slice_psf.py.in ../python/gen_python.sh
	@../python/gen_python.sh $(PYTHON) harp_slice_psf.py harp_slice_psf.py.in

harp_slice_spec.py : harp_slice_spec.py.in ../python/gen_python.sh
	@../python/gen_python.sh $(PYTHON) harp_slice_spec.py harp_slice_spec.py.in

