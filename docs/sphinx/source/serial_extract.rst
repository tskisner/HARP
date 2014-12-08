
.. _serial_extract:

Serial Spectral Extraction
==================================

The serial spectral extraction tools in HARP are useful debugging new data plugins or testing algorithms on small toy problems.  If you are performing production spectral extraction, you should use the :ref:`parallel extraction tools <parallel_extract>`.


.. doxygenclass:: harp::spec_slice_region
	:members:


.. doxygenclass:: harp::spec_slice
	:members:


.. doxygenfunction:: harp::sub_spec


.. doxygenfunction:: harp::accum_spec


.. doxygenfunction:: harp::noise_weighted_spec


.. doxygenfunction:: harp::inverse_covariance


.. doxygenfunction:: harp::resolution


.. doxygenfunction:: harp::accum_diag_resolution


.. doxygenfunction:: harp::extract


.. doxygenfunction:: harp::extract_slices




