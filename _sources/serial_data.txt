
.. _serial_data:

Serial Data I/O
==================================

For each type of object, HARP provides several derived classes which support specific data "formats".  Some of these classes actually generate data on demand, or simulate the data at construction and cache it internally.  There are two ways in which these classes should be instantiated, depending on your use case:

#.  If you are writing general tools that can work with any format for a particular type, use the factory method for that type.  Always wrap the returned raw pointer in a shared_ptr (see examples) to ensure proper destruction of the object.
#.  If you are instantiating a particular format of an object (for example so that you can create it and write data out), then just declare the variable.

The first of the two use cases above allows us to implement `inversion of control <http://en.wikipedia.org/wiki/Inversion_of_control>`_.  In other words, it allows us to write general analysis tools that work with any kind of data that provides the correct interface.


.. _serial_data_spec:

Spectra
-------------

A collection of spectra are handled in HARP by the "spec" class.  This represents one or more spectra with a common wavelength solution.  This class also specifies which spectra represent "sky" signal.

.. todo:
	Should we generalize the "sky" concept to an enumerated type spanning many object types?


.. doxygenclass:: harp::spec
	:members:


.. _serial_data_image:

Images
------------


.. doxygenclass:: harp::image
	:members:




.. _serial_data_psf:

Point Spread Functions
----------------------------

The HARP "psf" class represents the point spread function of the instrument and describes the projection of flux from individual wavelength points of each spectrum onto an image.


.. doxygenclass:: harp::psf
	:members:


.. _serial_data_targets:

Targets
----------------------------

The HARP "targets" class represents a collection of objects.


.. doxygenclass:: harp::targets
	:members:



