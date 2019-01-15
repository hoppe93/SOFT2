.. _module-equation:

@Equation
*********
This module specifies which equations of motion to solve, and the details of how
to solve them.

Types
^^^^^
This module has two secondary types. The secondary type, given in parentheses
after the configuration block name, specifies which set of equations to solve.
Note that the module settings listed under :ref:`module-equation-all-options`
are available regardless of the secondary type.

+--------------------+-----------------------------------------------------------+
| **Option**         | **Description**                                           |
+--------------------+-----------------------------------------------------------+
| ``guiding-center`` | Specifies the guiding-center equations of motion          |
+--------------------+-----------------------------------------------------------+
| ``particle``       | Specifies the particle ("full-orbit") equations of motion |
+--------------------+-----------------------------------------------------------+

Summary of options
^^^^^^^^^^^^^^^^^^

+-------------------------------+---------------------------------------------------------------+
| **Option**                    | **Description**                                               |
+-------------------------------+---------------------------------------------------------------+
| :option:`@Equation method`    | Numerical method used to solve equations of motion.           |
+-------------------------------+---------------------------------------------------------------+
| :option:`@Equation tolerance` | Relative tolerance of ODE solver (if adaptive solver is used) |
+-------------------------------+---------------------------------------------------------------+

Example configuration
^^^^^^^^^^^^^^^^^^^^^

This example illustrates how to configure the guiding-center equations of motion::

   @Equation gc (guiding-center) {
       method    = rkdp45;
       tolerance = 1e-9;
   }

.. _module-equation-all-options:

All options
^^^^^^^^^^^

.. program:: @Equation

.. option:: method

   :Default value: ``rkdp45``
   :Allowed values: ``rkdp45``

   Specifies the numerical method to use for solving the ODE that consititutes
   the guiding-center equations of motion. Currently, only one such numerical
   method is available, namely an adaptive Runge-Kutta solver of order 4, with
   an error estimator of order 5. The method uses Dormand-Prince coefficients,
   and is hence often known as the RKDP45 method [#rkdp45]_.

.. option:: tolerance

   :Default value: ``1e-9``
   :Allowed values: :math:`\epsilon < \text{tolerance} < 1`

   Sets the relative tolerance to use if an adaptive ODE solver is used to solve
   the equations of motion. The tolerance must be smaller than one, and greater
   than the machine epsilon :math:`\epsilon`, which is
   :math:`\approx 2\cdot 10^{-16}` using double-precision floating point numbers
   and :math:`\approx 10^{-7}` using single-precision.
 
.. [#rkdp45] Runge-Kutta-Dormand-Prince (RKDP45), Chapter XX: Integration of differential equations, *Numerical Recipes*, 3rd edition.

