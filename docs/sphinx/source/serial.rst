
.. _serial:

Serial Operations
==================================

The serial HARP library provides basic classes for reading and writing objects used in astrophysical data processing.  It also provides some useful math functions and a C++ interface to BLAS/LAPACK operations.  If python and numpy are available at compile time, python bindings to most classes will be created.


.. _serial-math:

Low-Level Math
------------------


C++ BLAS / LAPACK Interface
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

http://svn.boost.org/svn/boost/sandbox/numeric_bindings/libs/numeric/bindings/doc/html/index.html


.. _serial-io:

Data I/O
--------------

For each type of object, HARP provides several derived classes which support specific data "formats".  Some of these classes actually generate data on demand, or simulate the data at construction and cache it internally.  There are two ways in which these classes should be instantiated, depending on your use case:

#.  If you are writing general tools that can work with any format for a particular type, use the factory method for that type.  Always wrap the returned raw pointer in a shared_ptr (see examples) to ensure proper destruction of the object.
#.  If you are instantiating a particular format of an object (for example so that you can create it and write data out), then just declare the variable.

The first of the two use cases above allows us to implement `inversion of control <http://en.wikipedia.org/wiki/Inversion_of_control>`_.  In other words, it allows us to write general analysis tools that work with any kind of data that provides the correct interface.


.. _serial-spec:

Spectra
-------------

A collection of spectra are handled in HARP by the "spec" class.  This represents one or more spectra with a common wavelength solution.  This class also specifies which spectra represent "sky" signal.

.. todo:
	Should we generalize the "sky" concept to an enumerated type spanning many object types?


.. _serial-image:

Images
-------------





