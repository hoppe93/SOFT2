.. _module-distribution-numerical:

(numerical)
-----------
This module provides a way for the user to specify an arbitrary phase space
distribution function :math:`f(r, p, \xi)`.

Summary of options
^^^^^^^^^^^^^^^^^^

+-----------------------------------+------------------------------------------------------------------------------+
| **Option**                        | **Description**                                                              |
+-----------------------------------+------------------------------------------------------------------------------+
| :option:`numerical flippitchsign` | If true, flips the sign of the particle pitch.                               |
+-----------------------------------+------------------------------------------------------------------------------+
| :option:`numerical interptype`    | Interpolation method to use.                                                 |
+-----------------------------------+------------------------------------------------------------------------------+
| :option:`numerical logarithmize`  | Whether or not to interpolate in the logarithm of the distribution function. |
+-----------------------------------+------------------------------------------------------------------------------+
| :option:`numerical name`          | Name of file containing numerical distribution function.                     |
+-----------------------------------+------------------------------------------------------------------------------+

Example configuration
^^^^^^^^^^^^^^^^^^^^^
The following example configures a SOFT numerical distribution function::

   @DistributionFunction ourDistributionFunction (numerical) {
       name = "/path/to/distribution.mat";
   }

File layout
^^^^^^^^^^^
The format for numerical distribution functions in SOFT is such that a set of
momentum-space distribution functions are specified at a number of radii. The
file therefore contains two levels -- a top-level where the radial grid is
specified along with the number of radial points ``nr``, and ``nr`` groups
(HDF5) or structs (MAT), each containing the momentum-space distribution
function corresponding to a particular point on the radial grid.

The following variables must/may be present in the top-level of the distribution
function file:

+--------------+--------------+--------------------------------------------------------------------------------+
| **Variable** | **Required** | **Description**                                                                |
+--------------+--------------+--------------------------------------------------------------------------------+
| ``r``        | **Yes**      | (Dimensionful) minor radius (one-dimensional).                                 |
+--------------+--------------+--------------------------------------------------------------------------------+
| ``type``     | *No*         | Type of distribution function (must be ``distribution/soft``)                  |
+--------------+--------------+--------------------------------------------------------------------------------+

Additionally, on the top-level, ``nr`` (number of points in ``r``) number of
groups/structs named ``rX`` must be present, each containing the variables:

+--------------+--------------+--------------------------------------------------------------------------------+
| **Variable** | **Required** | **Description**                                                                |
+--------------+--------------+--------------------------------------------------------------------------------+
| ``f``        | **Yes**      | Distribution function :math:`f(r, p, \xi)` of size :math:`n_\xi \times n_p`.   |
+--------------+--------------+--------------------------------------------------------------------------------+
| ``p``        | **Yes**      | Momentum grid vector (one-dimensional).                                        |
+--------------+--------------+--------------------------------------------------------------------------------+
| ``xi``       | **Yes**      | Pitch (cosine of pitch angle) grid vector (one-dimensional).                   |
+--------------+--------------+--------------------------------------------------------------------------------+

.. note::

   In ``f``, the "inner" dimension (i.e. last index in C/C++/Python, first
   index in Matlab/Fortran) must be ``p``. In memory, the layout of the
   array should be::

      f(p0,xi0)  f(p1,xi0)  f(p2,xi0) ... f(pN,xi0)  f(p0,xi1)  f(p1,xi1) ...

All options
^^^^^^^^^^^

.. program:: numerical

.. option:: flippitchsign

   :Default value: ``no``
   :Allowed values: ``yes`` or ``no``

   If ``yes``, transforms :math:`f(p,\xi)\to f(p,-\xi)`. This is useful in case
   the distribution function is defined in an inconsistent way (i.e. if a
   positive parallel electric field accelerates particles in the anti-parallel
   direction).

.. option:: name

   :Default value: Nothing
   :Allowed values: Any valid file name

   Name of the file containing the distribution function.

.. option:: interptype

   :Default value: ``cubic``
   :Allowed values: ``cubic`` or ``linear``

   SOFT interpolates in the given distribution function to evaluate it at
   arbitrary points on the phase space grid. A linear interpolation scheme is
   always used to interpolate in the radial coordinate, but interpolation in
   the momentum coordinates (:math:`p` and :math:`\xi`) can either be done using
   bi-linear or bi-cubic splines.

.. option:: logarithmize

   :Default value: ``no``
   :Allowed values: ``yes`` or ``no``

   If ``yes``, interpolates in the logarithm of the distribution function
   instead of in the distribution function directly. This can aid in fitting
   sharply decaying ditsribution functions.

