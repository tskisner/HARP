
.. _install:

Installation
==================================

HARP uses an autotools-based build system.  If you are installing from a distribution tarball, then you do not need to actually have autotools installed.  If you are installing from a git checkout, then you need to have recent versions of automake, autoconf, and libtool installed.

.. _install_deps:

Dependencies
----------------

The first step to installing HARP is to make sure that you have all of the dependencies installed.  You should probably use a package manager to install these dependencies (see below).  The minimal required dependencies for the serial toolkit are CFITSIO, LAPACK, and a full installation of BOOST:

| http://heasarc.gsfc.nasa.gov/fitsio/
| http://www.netlib.org/lapack/
| http://www.boost.org

Note that many High Performance Computing (HPC) systems have vendor-tuned BLAS/LAPACK libraries, which should be used if possible.  If you would like to build the Python language interface, then you additionally must have installed python and numpy:

| http://www.python.org
| http://www.numpy.org

If you would like to use the parallel toolkit, then obviously you need an MPI implementation (e.g. MPICH or OpenMPI) that has been built with or is compatible with the serial compilers you are using.  You must also install Elemental:

| http://libelemental.org

The table below shows these dependencies and what optional features are enabled by each.  Dependencies detected at configure time.

+----------------------+--------------------------+------------------------+------------------------+------------------------+------------------------+
|                      | LAPACK / BOOST / CFITSIO | Python / Numpy         | MPI / boost::mpi       | Elemental              | mpi4py                 |
+======================+==========================+========================+========================+========================+========================+
| Serial Tools         | .. image:: yesmark.png   | .. image:: noxmark.png | .. image:: noxmark.png | .. image:: noxmark.png | .. image:: noxmark.png |
+----------------------+--------------------------+------------------------+------------------------+------------------------+------------------------+
| Python Bindings      | .. image:: yesmark.png   | .. image:: yesmark.png | .. image:: noxmark.png | .. image:: noxmark.png | .. image:: noxmark.png |
+----------------------+--------------------------+------------------------+------------------------+------------------------+------------------------+
| MPI Tools            | .. image:: yesmark.png   | .. image:: noxmark.png | .. image:: yesmark.png | .. image:: yesmark.png | .. image:: noxmark.png |
+----------------------+--------------------------+------------------------+------------------------+------------------------+------------------------+
| MPI Python Bindings  | .. image:: yesmark.png   | .. image:: yesmark.png | .. image:: yesmark.png | .. image:: yesmark.png | .. image:: yesmark.png |
+----------------------+--------------------------+------------------------+------------------------+------------------------+------------------------+


.. _install_ubuntu:

Installing Dependencies on Ubuntu
------------------------------------

On Ubuntu 14.04, you can get all dependencies except for Elemental by doing:

	| $> sudo apt-get install build-essential libopenblas-dev liblapack-dev libboost-all-dev python-numpy python-mpi4py libcfitsio3-dev

Then download and install Elemental.  For default options, you can just specify the install prefix and remember to do an out-of-source build.  In this example we downloaded and extracted the 0.84 tarball and made a build directory in the same place as that extracted directory:

	| $> mkdir build
	| $> cd build
	| $> cmake -D CMAKE_INSTALL_PREFIX=/home/<blah>/software ../Elemental-0.84-p1
	| $> make
	| $> make install


.. _install_osx:

Installing Dependencies on OS X
----------------------------------

The easiest way by far to get things like boost installed on OS X is to install macports:

	| https://www.macports.org

Then use macports to install all dependencies except Elemental:

	| $> sudo port install cfitsio py27-numpy boost +python27 +openmpi py27-mpi4py -mpich +openmpi

Then download and install Elemental.  For default options, you can just specify the install prefix and remember to do an out-of-source build.  In this example we downloaded and extracted the 0.84 tarball and made a build directory in the same place as that extracted directory:

	| $> mkdir build
	| $> cd build
	| $> cmake -D CMAKE_INSTALL_PREFIX=/Users/<blah>/software ../Elemental-0.84-p1
	| $> make
	| $> make install


.. _install_config:

Configuring HARP
--------------------

Before configuring HARP, you should determine:

#.  where you want to install HARP.
#.  which compilers you are using.
#.  whether any dependencies are installed in a non-standard place (and where that place is).







