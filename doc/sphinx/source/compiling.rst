.. _compiling:

Compiling
=========
Building SOFT2 is a two-step process. First, you must download and build the
SOFT support library, called `softlib`. After that, you can go ahead and build
the actual SOFT2 executable.

Preliminaries
-------------
Before building SOFT2 (or `softlib`), you need to make sure that you have the
following software installed on your system:

- `CMake <https://cmake.org/>`_ version 3.9 or later
- A C++17 compliant compiler and standard library (i.e. `gcc` >= 7 or equivalent)
- `GNU Scientific Library <https://www.gnu.org/software/gsl/>`_ version 2.0 or later
- `HDF5 <https://www.hdfgroup.org/>`_

If you intend to compile SOFT2 with MPI support, you also need an MPI library,
such as `OpenMPI <https://www.open-mpi.org/>`_ (not to be confused with OpenMP,
which is required as part of your compiler).

Building softlib
----------------
`softlib <https://github.com/hoppe93/softlib>`_ is a SOFT2 support library that
contains various helper classes and routines. You should place the directory
containing `softlib` in the directory that holds the SOFT2 root directory
(*note*: not in the SOFT2 root directory though).

Summary
*******
The whole build process is described in more detail below, but can be summarized
using the following chain of commands::

   $ cd /path/to/softlib
   $ mkdir build
   $ cd build
   $ cmake ../
   $ make

1. Create `build` directory
***************************
Once downloaded, navigate to the `softlib` directory and create a new directory
called `build`. Next, navigate into the newly created build directory.

2. Run CMake
************
Now it is time to run `cmake`::

   $ cmake ../

This tells CMake to read the file `CMakeLists.txt` located in the top directory
and execute all commands contained therein. There are a number of options that
can be given to CMake that can be used to configure `softlib`. **The default
options are the desired options in pretty much all builds.**

Each option is passed using the `-D` flag to CMake. The available options for `softlib` are

+------------------+----------------------------------------------------+-------------------+-------------+
| **Option**       | **Description**                                    | **Values**        | **Default** |
+------------------+----------------------------------------------------+-------------------+-------------+
| BUILD_TESTS      | Build unit tests                                   | ``ON`` or ``OFF`` | ``OFF``     |
+------------------+----------------------------------------------------+-------------------+-------------+
| DEBUG            | Build a debug version of `softlib`                 | ``ON`` or ``OFF`` | ``OFF``     |
+------------------+----------------------------------------------------+-------------------+-------------+
| OFFICIAL_MATLAB  | Link against official MATLAB libraries [#f1]_      | ``ON`` or ``OFF`` | ``OFF``     |
+------------------+----------------------------------------------------+-------------------+-------------+
| PRECISION_DOUBLE | Use double precision floating-point numbers [#f2]_ | ``ON`` or ``OFF`` | ``ON``      |
+------------------+----------------------------------------------------+-------------------+-------------+

.. [#f1] For SOFT2, ``PRECISION_DOUBLE`` should be turned on. Single-precision will yield very poor results.
.. [#f2] If this option is ``OFF``, MAT files can still be generated using the HDF5 API. If it is ``ON`` however, HDF5 can *not* be generated, and this requires an activated MATLAB installation to be present.

For example, to build ``softlib`` with an official MATLAB distribution, run::

   $ cmake ../ -DOFFICIAL_MATLAB=ON

instead of the ``cmake ../`` command shown above.

.. warning::

   When you run ``cmake ../``, CMake will output some information about which
   compiler and libraries it finds on your system. It can be worth to keep an
   eye on this output and verify that it agrees with your expectations. For
   example, if there are multiple compilers installed on the system, CMake may
   choose the wrong version.

.. note::

   To explicitly tell CMake which C/C++ compilers to use, run the following
   two commands before running CMake::

      $ export CC="/path/to/c-compiler"
      $ export CXX="/path/to/c++-compiler"

   (Note that you may have to clean up the ``build`` directory before running
   ``cmake`` again).

3. Run make
***********
In the final step of the build we run ``make``::

   $ make

Since the build can take some time, it is recommended to use the ``-j`` option
of ``make``. Passing this option to ``make`` along with a number ``N`` will cause
``make`` to run ``N`` compilations in parallel. For example, if your computer has
4 threads available, you can run::

   $ make -j 4

for optimal compilation speed.

Building SOFT2
--------------
Once ``softlib`` has been built it is time to download and build SOFT2.

Summary
*******
The whole build process is described in more detail below, but can be summarized
using the following chain of commands::

   $ cd /path/to/SOFT2/build
   $ cmake ../
   $ make

1. Run CMake
************
The SOFT2 code already comes with a ``build/`` directory. Navigate to that
directory and run::

   $ cmake ../

As with ``softlib``, there are a number of options that can be given to
``cmake`` to configure the SOFT2 build. The commands are specified to ``cmake``
using the ``-D`` command-line option. The following options are available:

+------------------+----------------------------------------------------+-------------------+-------------+
| **Option**       | **Description**                                    | **Values**        | **Default** |
+------------------+----------------------------------------------------+-------------------+-------------+
| BUILD_TESTS      | Build unit tests                                   | ``ON`` or ``OFF`` | ``OFF``     |
+------------------+----------------------------------------------------+-------------------+-------------+
| COLOR_TERMINAL   | Enable colored output                              | ``ON`` or ``OFF`` | ``ON``      |
+------------------+----------------------------------------------------+-------------------+-------------+
| DEBUG            | Include debug symbols in the binary                | ``ON`` or ``OFF`` | ``OFF``     |
+------------------+----------------------------------------------------+-------------------+-------------+
| OPTIMIZE_NATIVE  | Apply native compiler optimizations                | ``ON`` or ``OFF`` | ``ON``      |
+------------------+----------------------------------------------------+-------------------+-------------+
| PROFILING        | Compile with profiler flags                        | ``ON`` or ``OFF`` | ``OFF``     |
+------------------+----------------------------------------------------+-------------------+-------------+
| WITH_MPI         | Compile with MPI support                           | ``ON`` or ``OFF`` | ``OFF``     |
+------------------+----------------------------------------------------+-------------------+-------------+

For example, to disable colored terminal output (useful if you're redirecting
stdout to a text file for example), run ``cmake`` as::

   $ cmake ../ -DCOLOR_TERMINAL=OFF

.. warning::

   When you run ``cmake ../``, CMake will output some information about which
   compiler and libraries it finds on your system. It can be worth to keep an
   eye on this output and verify that it agrees with your expectations. For
   example, if there are multiple compilers installed on the system, CMake may
   choose the wrong version.

.. note::

   To explicitly tell CMake which C/C++ compilers to use, run the following
   two commands before running CMake::

      $ export CC="/path/to/c-compiler"
      $ export CXX="/path/to/c++-compiler"

   (Note that you may have to clean up the ``build`` directory before running
   ``cmake`` again).

2. Run make
***********
In the final step of the build we run ``make``::

   $ make

Since the build can take some time, it is recommended to use the ``-j`` option
of ``make``. Passing this option to ``make`` along with a number ``N`` will cause
``make`` to run ``N`` compilations in parallel. For example, if your computer has
4 threads available, you can run::

   $ make -j 4

for optimal compilation speed.

Complete Ubuntu example
-----------------------
On a clean Ubuntu 18.10 install, the complete installation procedure would
look as follows::

$ apt update
$ apt install build-essential libgsl-dev cmake libhdf5-serial-dev git
$ mkdir SOFT
$ cd SOFT
$ git clone https://github.com/hoppe93/softlib.git
$ git clone https://github.com/hoppe93/SOFT2.git
$ cd softlib
$ mkdir build
$ cd build
$ cmake ../
$ make
$ cd ../SOFT2/build
$ cmake ../
$ make

