.. _module-equation-particle:

@EquationParticle
*****************
The module specifying the particle equations of motion. When this equation is
used by the :ref:`module-particlepusher`, SOFT follows particles (i.e. not
guiding-centers). This means that regular (cartesian) particle coordinates are
used for evaluating quantities such as detected radiation.

Summary of options
^^^^^^^^^^^^^^^^^^

+---------------------------------------+---------------------------------------------------------------+
| **Option**                            | **Description**                                               |
+---------------------------------------+---------------------------------------------------------------+
| :option:`@EquationParticle method`    | Numerical method used to solve equations of motion.           |
+---------------------------------------+---------------------------------------------------------------+
| :option:`@EquationParticle tolerance` | Relative tolerance of ODE solver (if adaptive solver is used) |
+---------------------------------------+---------------------------------------------------------------+

Example configuration
^^^^^^^^^^^^^^^^^^^^^

This example illustrates how to configure the particle equations of motion::

   @EquationParticle particle {
       method    = rkdp45;
       tolerance = 1e-7;
   }

All options
^^^^^^^^^^^

.. program:: @EquationParticle

.. option:: method

   :Default value: ``rkdp45``
   :Allowed values: ``rkdp45``

   Specifies the numerical method to use for solving the ODE that consititutes
   the particle equations of motion. Currently, only one such numerical
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

