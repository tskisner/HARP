
.. _parallel_data:

Parallel Data I/O
==================================

HARP has distributed versions of the spec, targets, image, and psf classes.  These simply read the object on one process and broadcast a full copy to all other processes.



.. doxygenclass:: harp::mpi_spec
	:members:



.. doxygenclass:: harp::mpi_image
	:members:



.. doxygenclass:: harp::mpi_psf
	:members:



.. doxygenclass:: harp::mpi_targets
	:members:



