
# Initialize compile and linking flags

AM_CPPFLAGS = 
AM_CXXFLAGS =
AM_LDFLAGS =
AM_LIBS =

# Append compile and linking flags based on configure variables
include $(top_srcdir)/src/build_options.am


bin_SCRIPTS =

BUILT_SOURCES =

CLEANFILES =

bin_PROGRAMS = harp_example_imgstack harp_extract harp_validate

harp_example_imgstack_SOURCES = harp_example_imgstack.cpp
harp_example_imgstack_LDADD = $(top_builddir)/src/libharp/libharp.la

harp_extract_SOURCES = harp_extract.cpp
harp_extract_LDADD = $(top_builddir)/src/libharp/libharp.la

harp_validate_SOURCES = harp_validate.cpp
harp_validate_LDADD = $(top_builddir)/src/libharp/libharp.la


if HAVE_AM_PYTHON
	bin_SCRIPTS += harp_slice_psf.py harp_slice_spec.py
	BUILT_SOURCES += harp_slice_psf.py harp_slice_spec.py
	CLEANFILES += harp_slice_psf.py harp_slice_spec.py
endif

harp_slice_psf.py : harp_slice_psf.py.in ../python/gen_python.sh
	@../python/gen_python.sh $(PYTHON) harp_slice_psf.py harp_slice_psf.py.in

harp_slice_spec.py : harp_slice_spec.py.in ../python/gen_python.sh
	@../python/gen_python.sh $(PYTHON) harp_slice_spec.py harp_slice_spec.py.in

