.. _module-distribution-numerical:

(numerical)
-----------
This module provides a way for the user to specify an arbitrary phase space
distribution function :math:`f(r, p, \xi)`.

Summary of options
^^^^^^^^^^^^^^^^^^

+----------------------------------+------------------------------------------------------------------------------+
| **Option**                       | **Description**                                                              |
+----------------------------------+------------------------------------------------------------------------------+
| :option:`numerical name`         | Name of file containing numerical distribution function.                     |
+----------------------------------+------------------------------------------------------------------------------+
| :option:`numerical interptype`   | Interpolation method to use.                                                 |
+----------------------------------+------------------------------------------------------------------------------+
| :option:`numerical logarithmize` | Whether or not to interpolate in the logarithm of the distribution function. |
+----------------------------------+------------------------------------------------------------------------------+

Example configuration
^^^^^^^^^^^^^^^^^^^^^
The following example configures a SOFT numerical distribution function::

   @DistributionFunction ourDistributionFunction (numerical) {
       name = "/path/to/distribution.mat";
   }

File layout
^^^^^^^^^^^
The following variables must/may be present in the distribution function file.

+--------------+--------------+---------------------------------------------------------------------------+
| **Variable** | **Required** | **Description**                                                           |
+--------------+--------------+---------------------------------------------------------------------------+
| ``f``        | **Yes**      | Distribution function :math:`f(r, p, \xi)` of size :math:`n_p n_\xi n_r`. |
+--------------+--------------+---------------------------------------------------------------------------+
| ``fp0``      | *No*         | Verification vector; :math:`f(r_0, p, \xi_0)`.                            |
+--------------+--------------+---------------------------------------------------------------------------+
| ``fr0``      | *No*         | Verification vector; :math:`f(r, p_0, \xi_0)`.                            |
+--------------+--------------+---------------------------------------------------------------------------+
| ``fxi0``     | *No*         | Verification vector; :math:`f(r_0, p_0, \xi)`.                            |
+--------------+--------------+---------------------------------------------------------------------------+
| ``p``        | **Yes**      | Momentum grid vector with :math:`n_p` elements.                           |
+--------------+--------------+---------------------------------------------------------------------------+
| ``punits``   | **Yes**      | Unit of the given momentum.                                               |
+--------------+--------------+---------------------------------------------------------------------------+
| ``r``        | **Yes**      | Radial grid vector with :math:`n_r` elements.                             |
+--------------+--------------+---------------------------------------------------------------------------+
| ``xi``       | **Yes**      | Pitch (cosine of pitch angle) grid vector with :math:`n_\xi` elements.    |
+--------------+--------------+---------------------------------------------------------------------------+

All options
^^^^^^^^^^^

.. program:: numerical

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

