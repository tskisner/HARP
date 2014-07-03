
AM_CPPFLAGS = -I$(top_srcdir)/src/libharp/math -I$(top_srcdir)/src/libharp/data -I$(top_srcdir)/src/libharp/extract -I$(top_srcdir)/src/libharp/plugins -I$(top_srcdir)/src/libharp -I$(top_srcdir)/src/apps $(CFITSIO_CPPFLAGS) $(BOOST_CPPFLAGS)

AM_LDFLAGS =

AM_LIBS = $(top_builddir)/src/libharp/libharp.la $(BOOST_PROGRAM_OPTIONS_LDFLAGS) $(BOOST_PROGRAM_OPTIONS_LIBS) $(BOOST_REGEX_LDFLAGS) $(BOOST_REGEX_LIBS)

if HAVE_AM_PYTHON
  AM_LIBS += $(top_builddir)/src/python/boost_numpy/libboost_numpy.la $(BOOST_PYTHON_LDFLAGS) $(BOOST_PYTHON_LIBS) $(BOOST_SERIALIZATION_LDFLAGS) $(BOOST_SERIALIZATION_LIBS) $(PYTHON_LDFLAGS) $(PYTHON_LIBS)
else
  AM_LIBS += $(BOOST_SERIALIZATION_LDFLAGS) $(BOOST_SERIALIZATION_LIBS)
endif

AM_LIBS += $(CFITSIO) $(LAPACK_LIBS) $(BLAS_LIBS) $(FCLIBS) $(OPENMP_CXXFLAGS) -lm

AM_CXXFLAGS = $(OPENMP_CXXFLAGS)

if ENABLE_AM_DEBUG
else
  AM_CXXFLAGS += -DNDEBUG
endif

if USE_AM_ALLSTATIC
  AM_LDFLAGS += -all-static
endif

bin_SCRIPTS =
BUILT_SOURCES =
CLEANFILES =

bin_PROGRAMS = harp_example_imgstack harp_extract harp_validate

harp_example_imgstack_SOURCES = harp_example_imgstack.cpp
harp_example_imgstack_LDADD = $(AM_LIBS)

harp_extract_SOURCES = harp_extract.cpp
harp_extract_LDADD = $(AM_LIBS)

harp_validate_SOURCES = harp_validate.cpp
harp_validate_LDADD = $(AM_LIBS)


if HAVE_AM_PYTHON
	bin_SCRIPTS += harp_slice_psf.py harp_slice_spec.py
	BUILT_SOURCES += harp_slice_psf.py harp_slice_spec.py
	CLEANFILES += harp_slice_psf.py harp_slice_spec.py
endif

harp_slice_psf.py : harp_slice_psf.py.in ../python/gen_python.sh
	@../python/gen_python.sh $(PYTHON) harp_slice_psf.py harp_slice_psf.py.in

harp_slice_spec.py : harp_slice_spec.py.in ../python/gen_python.sh
	@../python/gen_python.sh $(PYTHON) harp_slice_spec.py harp_slice_spec.py.in

