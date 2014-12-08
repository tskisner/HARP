
.. _serial_utils:

Serial Utilities
==================================

The serial HARP library provides some useful math functions and a C++ interface to BLAS/LAPACK operations.  It also provides helper functions for reading / writing FITS data structures.


.. _serial_utils_math:

Low-Level Math
------------------


C++ BLAS / LAPACK Interface
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

http://svn.boost.org/svn/boost/sandbox/numeric_bindings/libs/numeric/bindings/doc/html/index.html


Linear Algebra Helper Functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doxygenfunction:: harp::eigen_decompose

.. doxygenfunction:: harp::eigen_compose

.. doxygenfunction:: harp::column_norm

.. doxygenfunction:: harp::apply_norm

.. doxygenfunction:: harp::apply_inverse_norm

.. doxygenfunction:: harp::norm

.. doxygenfunction:: harp::sparse_mv_trans


.. _serial_utils_fits:

FITS Helper Tools
---------------------
