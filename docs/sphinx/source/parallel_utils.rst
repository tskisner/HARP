
.. _parallel_utils:

Parallel Utilities
==================================

The parallel HARP library uses the Elemental library for dense linear algebra, and also includes a custom class for handling distributed sparse matrices.


.. _parallel_utils_comm:

Communication Helper Functions
----------------------------------



.. _parallel_utils_math:

Linear Algebra
------------------


.. doxygenfunction:: harp::mpi_matrix_zero

.. doxygenfunction:: harp::local_matrix_zero

.. doxygenfunction:: harp::mpi_eigen_decompose

.. doxygenfunction:: harp::mpi_eigen_compose

.. doxygenfunction:: harp::mpi_column_norm

.. doxygenfunction:: harp::mpi_apply_norm

.. doxygenfunction:: harp::mpi_apply_inverse_norm

.. doxygenfunction:: harp::mpi_norm

.. doxygenfunction:: harp::mpi_sparse_mv_trans

.. doxygenfunction:: harp::mpi_gang_distribute

.. doxygenfunction:: harp::mpi_gang_accum


.. doxygenclass:: harp::mpi_matrix_sparse_block
	:members:


.. doxygenclass:: harp::mpi_matrix_sparse
	:members:


