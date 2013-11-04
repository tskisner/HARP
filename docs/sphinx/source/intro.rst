
.. _intro:

Introduction
===============

The High performance Astrophysical Reconstruction and Processing (HARP) software suite is designed to serve as both a platform for production...  HARP has a number of design goals, including:

* Abstraction of instrument-specific features (data formats and data selection) into modular pieces so that it is easy to extend support to multiple receivers.
* Fast simulations of timestream data and processing, to enable Monte Carlo analyses.
* Scalability to the largest production supercomputers.
* Provide a rich collection of map making operations which are important for modern experiments.

.. _introserial:

Serial Operations
---------------------

HARP provides several base classes for reading and writing relevant data objects.  In these terms, "reading" might include simulating data as needed.

As a toy example, HARP includes serial spectral extraction functions.  These are mainly for testing algorithms and are not performant enough for real data processing.  See the :ref:`extractmpi` section for parallel tools for spectral extraction.


.. _intrompi:

Parallel Operations
------------------------



