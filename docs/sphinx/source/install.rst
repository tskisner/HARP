
.. _install:

Installation
==================================

HARP uses an autotools-based build system.  If you are installing from a distribution tarball, then you do not need to actually have autotools installed.  If you are installing from a git checkout, then you need to have recent versions of automake, autoconf, and libtool installed.


The table below shows the dependencies that TOAST requires and what optional features are enabled by each.  The minimal dependencies are BOOST (full installation), LAPACK, CFITSIO, and the MOAT library.  Additional dependencies detected at configure time will enable optional features.

+----------------------+------------------------+------------------------+------------------------+------------------------+------------------------+------------------------+------------------------+------------------------+
|                      | LAPACK / BOOST         | Python / Numpy         | wcslib                 | MPI / boost::mpi       | Elemental              | HEALPix                | GetData                | HDF5                   |
|                      | MOAT / CFITSIO         |                        |                        |                        |                        |                        |                        |                        |
+======================+========================+========================+========================+========================+========================+========================+========================+========================+
| Serial Map-making    | .. image:: yesmark.png | .. image:: noxmark.png | .. image:: noxmark.png | .. image:: noxmark.png | .. image:: noxmark.png | .. image:: noxmark.png | .. image:: noxmark.png | .. image:: noxmark.png | 
| Tools                |                        |                        |                        |                        |                        |                        |                        |                        |
+----------------------+------------------------+------------------------+------------------------+------------------------+------------------------+------------------------+------------------------+------------------------+
| Python Bindings      | .. image:: yesmark.png | .. image:: yesmark.png | .. image:: noxmark.png | .. image:: noxmark.png | .. image:: noxmark.png | .. image:: noxmark.png | .. image:: noxmark.png | .. image:: noxmark.png |
+----------------------+------------------------+------------------------+------------------------+------------------------+------------------------+------------------------+------------------------+------------------------+
| WCS Projections      | .. image:: yesmark.png | .. image:: noxmark.png | .. image:: yesmark.png | .. image:: noxmark.png | .. image:: noxmark.png | .. image:: noxmark.png | .. image:: noxmark.png | .. image:: noxmark.png |
+----------------------+------------------------+------------------------+------------------------+------------------------+------------------------+------------------------+------------------------+------------------------+
| MPI Map-making       | .. image:: yesmark.png | .. image:: noxmark.png | .. image:: noxmark.png | .. image:: yesmark.png | .. image:: noxmark.png | .. image:: noxmark.png | .. image:: noxmark.png | .. image:: noxmark.png |
| Tools                |                        |                        |                        |                        |                        |                        |                        |                        |
+----------------------+------------------------+------------------------+------------------------+------------------------+------------------------+------------------------+------------------------+------------------------+
| MPI Pixel            | .. image:: yesmark.png | .. image:: noxmark.png | .. image:: noxmark.png | .. image:: yesmark.png | .. image:: yesmark.png | .. image:: yesmark.png | .. image:: noxmark.png | .. image:: noxmark.png |
| Noise Tools          |                        |                        |                        |                        |                        |                        |                        |                        |
+----------------------+------------------------+------------------------+------------------------+------------------------+------------------------+------------------------+------------------------+------------------------+
| Internal Dirfile     | .. image:: yesmark.png | .. image:: noxmark.png | .. image:: noxmark.png | .. image:: noxmark.png | .. image:: noxmark.png | .. image:: noxmark.png | .. image:: yesmark.png | .. image:: noxmark.png |
| Formats              |                        |                        |                        |                        |                        |                        |                        |                        |
+----------------------+------------------------+------------------------+------------------------+------------------------+------------------------+------------------------+------------------------+------------------------+
| Planck Data          | .. image:: yesmark.png | .. image:: noxmark.png | .. image:: noxmark.png | .. image:: noxmark.png | .. image:: noxmark.png | .. image:: noxmark.png | .. image:: noxmark.png | .. image:: noxmark.png |
| Formats              |                        |                        |                        |                        |                        |                        |                        |                        |
+----------------------+------------------------+------------------------+------------------------+------------------------+------------------------+------------------------+------------------------+------------------------+
| PolarBear Data       | .. image:: yesmark.png | .. image:: noxmark.png | .. image:: noxmark.png | .. image:: noxmark.png | .. image:: noxmark.png | .. image:: noxmark.png | .. image:: noxmark.png | .. image:: yesmark.png |
| Formats              |                        |                        |                        |                        |                        |                        |                        |                        |
+----------------------+------------------------+------------------------+------------------------+------------------------+------------------------+------------------------+------------------------+------------------------+
| Blast Data           | .. image:: yesmark.png | .. image:: noxmark.png | .. image:: noxmark.png | .. image:: noxmark.png | .. image:: noxmark.png | .. image:: noxmark.png | .. image:: yesmark.png | .. image:: noxmark.png |
| Formats              |                        |                        |                        |                        |                        |                        |                        |                        |
+----------------------+------------------------+------------------------+------------------------+------------------------+------------------------+------------------------+------------------------+------------------------+
| EBEX Data            | .. image:: yesmark.png | .. image:: noxmark.png | .. image:: noxmark.png | .. image:: noxmark.png | .. image:: noxmark.png | .. image:: noxmark.png | .. image:: yesmark.png | .. image:: noxmark.png |
| Formats              |                        |                        |                        |                        |                        |                        |                        |                        |
+----------------------+------------------------+------------------------+------------------------+------------------------+------------------------+------------------------+------------------------+------------------------+

Additionally, in order to actually select data for the various experiment data formats, usually Python is required to enable data selection classes and interfaces written in Python.  The requirements above represent the run-time dependencies to load a previously dumped run file and load it into the library.  Additionally, the Fortran 2003 language bindings for TOAST require a compatible Fortran compiler.


.. _install_ubuntu:

Installation on Ubuntu
--------------------------




.. _install_osx:

Installation on OS X
--------------------------

Apple no longer ships a Fortran compiler, so these instructions assume you are using clang and clang++ with Fortran bindings disabled.  The easiest way to install TOAST dependencies is by using Macports ()





Installation
==================================

HARP uses an autotools-based build system.  If you are installing from a distribution tarball, then you do not need to actually have autotools installed.  If you are installing from a git checkout, then you need to have recent versions of automake, autoconf, and libtool installed.


Dependencies
----------------

The first step to installing HARP is to make sure that you have all of the dependencies installed.  The minimal required dependencies for the serial toolkit are CFITSIO, LAPACK, and a full installation of BOOST:

| http://heasarc.gsfc.nasa.gov/fitsio/
| http://www.netlib.org/lapack/
| http://www.boost.org

Note that many High Performance Computing (HPC) systems have vendor-tuned BLAS/LAPACK libraries, which should be used if possible.  If you would like to build the Python language interface, then you additionally must have installed python and numpy:

| http://www.python.org
| http://www.numpy.org

If you would like to use the parallel toolkit, then obviously you need an MPI implementation (e.g. MPICH or OpenMPI) that has been built with or is compatible with the serial compilers you are using.  You must also install Elemental and Clique:

| http://libelemental.org
| https://github.com/poulson/Clique


Instructions for Ubuntu GNU/Linux 12.04 LTS
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


Instructions for Apple OS X
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^




Configuring HARP
--------------------

Before configuring HARP, you should determine:

#.  where you want to install HARP.
#.  which compilers you are using.
#.  whether any dependencies are installed in a non-standard place (and where that place is).




