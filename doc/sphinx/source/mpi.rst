
.. _mpi:

Running with MPI
================
SOFT2 is parallelized in two levels. Shared memory parallelization using the
OpenMP framework is a core part of SOFT that is always enabled. By default,
SOFT will run computations on as many threads as reported by OpenMP, which
permits full utilization of single-node computer systems, such as laptops and
desktops.

To run on distributed computer systems, however, the *Message Passing Interface*
(MPI) must be used to achieve full parallelization. MPI runs one SOFT2 process
on each node of the computer system, and allows them to communicate. Each
process can (and will by default) in turn run computations on several OpenMP
threads. In contrast to OpenMP, MPI is not enabled in SOFT by default and
requires an MPI library to be installed.

In what follows we will describe the steps necessary to run SOFT2 with MPI and
the options available for simulations.

Compiling
---------
Before compiling SOFT2, make sure that an MPI library is installed on your
system. Note that MPI is a standard, meaning that there are several
implementations available. SOFT2 doesn't care which implementation you use, and
so you are free to choose from any of the available. Some popular
implementations are

- `OpenMPI <https://www.open-mpi.org/>`_
- `MPICH <https://www.mpich.org/>`_
- `Intel MPI <https://software.intel.com/en-us/mpi-library>`_

If you are using the Intel C++ compiler, you should prefer to also use Intel
MPI.

To compile SOFT2 with MPI, simply run ``cmake`` with the flag ``-DWITH_MPI=ON``::

   cmake .. -DWITH_MPI=ON

This is the only modification required when compiling SOFT2, and you should
otherwise adhere to the steps described in :ref:`compiling`.

.. note::

   Note that MPI support need only be configured for SOFT2 and *not* for
   ``softlib``. Thus, ``softlib`` should be built in exactly the same way,
   regardless whether MPI support is desired or not.

Setting up a simulation
-----------------------
Any SOFT2 simulation that can be run without MPI should also be possible to run
when MPI support is enabled, without any modification to the simulation script.
To make the most of the MPI parallelization, SOFT2 also provides a few
additional options when compiled with MPI support. These options are however
primarily for advanced users, and unnecessary in most cases.

Specify order of phase space
****************************
When running with MPI, SOFT2 parallelizes the computation by dividing one of the
three phase space parameters (radius and momentum parameters) evenly among the
MPI processes. Which parameter to divide is determined by the option
:option:`@ParticleGenerator mpi_distribute_mode`, and by default this is set to
the parameter with the largest number of grid points.

Green's function output
***********************
The only type of :ref:`module-radiationoutput` that has a special MPI option is
the :ref:`module-ro-green` module. This is because the Green's function module
is typically the most memory-intensive module, and a special option has
therefore been added to help better utilize computer clusters.

By default, :ref:`module-ro-green` stores a version of the Green's function in
each MPI process and on each thread. At the end of the run, all local Green's
functions are reduced to a single Green's function on the main thread of the
root MPI process, which is then written to disk. This mode can also be manually
specified by setting ``mpi_mode=contiguous`` in the :ref:`module-ro-green`
module.

To better utilize the available memory however, the Green's function can be
split into several independent pieces, with each MPI process working on its own
piece of the function. At the end of the run, each process then writes its piece
to a separate file. To use the Green's function, the separate pieces should in
principle be loaded and put in sequential order in memory (although one can
think of other techniques for using the Green's function that doesn't require
all pieces to reside in the same memory space). This mode is designed for
generating very large Green's functions, and is selected via the option
``mpi_mode=chunked`` in the :ref:`module-ro-green` module.

.. important::

   Make sure that the phase space parameter that is divided among the MPI
   processes is a part of the Green's function, i.e. is in the Green's functions
   format string. Which parameter that is divided among the MPI processes is
   set using the :option:`@ParticleGenerator mpi_distribute_mode` option.

Running
-------
To run SOFT2 in MPI mode, simply use the ``mpirun`` command::

   $ mpirun -n N /path/to/soft inputfile

where ``N`` is the number of MPI processes to launch. Note that MPI only
provides any benefits when run on a computer cluster, in which case a scheduler
(such as SLURM) must be used. An example configuration for a SLURM system is
presented below, and a similar configuration would most likely be necessary with
other scheduler systems as well.

Using SLURM
***********
SLURM is one of the more popular so-called *schedulers* which is used to submit
jobs to a computer cluster. The example below illustrates how a SLURM
configuration for a SOFT2 run could look like. Make sure to thoroughly read the
documentation for the cluster you are running on though, as desired settings may
vary from system to system.

.. code-block:: bash

   #!/bin/bash
   #
   # Number of nodes to run on.
   #SBATCH -N 32
   #
   # Number of processor cores. If each node has 32 cores, and you have
   # requested 32 nodes above, then you should specify 32*32 = 1024 here.
   #SBATCH -n 1024
   #
   # Desired memory in MB per processor core.
   #SBATCH --mem-per-cpu=4000
   #
   # Maximum execution time of the job.
   #SBATCH --time=00:30:00
   #
   # Queue to submit job to (these names are usually specific to the cluster).
   #SBATCH -p sched_express
   #
   # Name of file to direct stdout/stderr to.
   #SBATCH -o soft-%j.out

   module load gsl hdf5 mpi

   # Run SOFT2
   mpirun /path/to/soft inputfile


