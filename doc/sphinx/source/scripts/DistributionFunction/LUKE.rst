.. _module-distribution-luke:

(luke)
------
LUKE solves the bounce-averaged Fokker-Planck equation in one spatial dimension
and two momentum dimensions.

Summary of options
^^^^^^^^^^^^^^^^^^

+----------------------------------+------------------------------------------------------------------------------+
| **Option**                       | **Description**                                                              |
+----------------------------------+------------------------------------------------------------------------------+
| :option:`luke interptype`        | Interpolation method to use.                                                 |
+----------------------------------+------------------------------------------------------------------------------+
| :option:`dream logarithmize`     | Whether or not to interpolate in the logarithm of the distribution function. |
+----------------------------------+------------------------------------------------------------------------------+
| :option:`luke name`              | Name of file containing distribution function.                               |
+----------------------------------+------------------------------------------------------------------------------+

Example configuration
^^^^^^^^^^^^^^^^^^^^^

.. code-block::

   @DistributionFunction luke {
       name = "LUKEDistributionFunction.mat";
   }

File layout
^^^^^^^^^^^
The file containing the LUKE distribution function must contain the variables
listed in the table below:

+----------------+------------------------------------------------------------------------------------------+
| **Variable**   | **Description**                                                                          |
+----------------+------------------------------------------------------------------------------------------+
| ``betath_ref`` | Reference normalized thermal speed :math:`\beta_\mathrm{th,ref}` (for de-normalization). |
+----------------+------------------------------------------------------------------------------------------+
| ``f``          | Distribution function :math:`f(\psi,\xi,p)`.                                             |
+----------------+------------------------------------------------------------------------------------------+
| ``mhu``        | Pitch grid in the variable :math:`\xi = \cos\theta_{\mathrm{p}}`.                        |
+----------------+------------------------------------------------------------------------------------------+
| ``pn``         | Momentum grid, with momentum normalized to :math:`\beta_{\mathrm{th,ref}}`.              |
+----------------+------------------------------------------------------------------------------------------+
| ``xrhoG``      | Radial grid in the variable :math:`\psi/\psi_a`, i.e. normalized poloidal flux.          |
+----------------+------------------------------------------------------------------------------------------+

All options
^^^^^^^^^^^

.. program:: luke

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

.. option:: name

   :Default value: Nothing
   :Allowed values: Any valid file name

   Name of the file containing the distribution function.

