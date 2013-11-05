
.. _install:

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




