.. _module-particlepusher:

@ParticlePusher
***************
The particle pusher generates orbits---both particle and guiding-center
orbits---which can be fed to the various tools of SOFT that use them to, for
example, compute radiation images. To output the orbits, use the
:ref:`module-orbits` tool.

Summary of options
^^^^^^^^^^^^^^^^^^

+----------------------------------------------------+--------------------------------------------------------------------------+
| **Option**                                         | **Description**                                                          |
+----------------------------------------------------+--------------------------------------------------------------------------+
| :option:`@ParticlePusher equation`                 | Specifies which equation to solve.                                       |
+----------------------------------------------------+--------------------------------------------------------------------------+
| :option:`@ParticlePusher force_numerical_jacobian` | Forces numerical computation of the guiding-center Jacobian.             |
+----------------------------------------------------+--------------------------------------------------------------------------+
| :option:`@ParticlePusher nt`                       | Number of time steps to take.                                            |
+----------------------------------------------------+--------------------------------------------------------------------------+
| :option:`@ParticlePusher nudgevalue`               | Override :math:`\mathrm{d}\rho` used to compute GC Jacobian numerically. |
+----------------------------------------------------+--------------------------------------------------------------------------+
| :option:`@ParticlePusher time`                     | Duration to follow orbits for.                                           |
+----------------------------------------------------+--------------------------------------------------------------------------+
| :option:`@ParticlePusher timeunit`                 | Unit of the time specified to :option:`@ParticlePusher time`.            |
+----------------------------------------------------+--------------------------------------------------------------------------+

Example configuration
^^^^^^^^^^^^^^^^^^^^^

The following is the default configuration of the particle pusher object. It
causes the pusher to generate guiding-center orbits that are traced during one
poloidal transit time. The resulting orbit consists of 1000 time steps::

   @ParticlePusher __default__ {
       equation                 = guiding-center;
       timeunit                 = poloidal;
       time                     = 1;
       nt                       = 1e3;
       force_numerical_jacobian = no;
       nudgevalue               = __default__;
   }

Options
^^^^^^^

.. program:: @ParticlePusher

.. option:: equation

   :Default value: ``guiding-center``
   :Allowed values: Name of any :ref:`module-equation` configuration block.

   Specifies which equation to use. The two equations ``guiding-center`` and
   ``particle`` are pre-defined by this module with default settings. The
   settings of these equations can be overridden by adding configuration blocks
   for them.

.. option:: force_numerical_jacobian

   :Default value: ``no``
   :Allowed values: A boolean value; ``yes`` or ``no``.

   Force the guiding-center Jacobian to be computed numerically. If ``no``, an
   analytical expression will instead be used for the Jacobian. The numerical
   approach is slower, prone to instabilities and in general discouraged.

.. option:: nt

   :Default value: ``1e3``, i.e. 1000 points
   :Allowed values: Any positive integer.

   The number of time steps per orbit. The orbit quantities (particle position
   and momentum) are given in a uniformly distributed set of time points between
   0 and the maximum time, set by the :option:`@ParticlePusher time` parameter.

.. option:: nudgevalue

   :Default value: ``__default__`` (see description below)
   :Allowed values: Any real number.

   To compute the guiding-center Jacobian numerically, the derivatives of the
   particle position with respect to the initial radial location must be taken.
   This parameter is the distance :math:`\Delta\rho` by which each orbit is
   nudged in order to evaluate the derivative using a finite difference method.

.. option:: time

   :Default value: ``1``
   :Allowed values: Any positive real number.

   Final time point in which to evaluate orbit.

   **Note:** The :ref:`module-radiation` module expects this parameter to be
   set to ``1``, and the :option:`@ParticlePusher timeunit` parameter to be set
   to ``poloidal``. Note also that the :ref:`module-radiation` automatically
   discards the final time point to prevent double-counting.

.. option:: timeunit

   :Default value: ``poloidal``, i.e. poloidal transit time
   :Allowed values: ``poloidal`` and ``seconds``.

   The unit of the :option:`@ParticlePusher time` parameter. If this parameter
   is set to ``poloidal``, the :option:`@ParticlePusher time` gives the number
   of poloidal transits for which each particle should be followed. If this
   parameter is set to ``seconds``, :option:`@ParticlePusher time` gives the
   number of seconds to follow each orbit.

   *The option 'poloidal' also works for particle/full orbits, even though the
   concept of poloidal transit time is much less well-defined for in such cases.
   To overcome this poor definition, SOFT first solves the corresponding
   guiding-center orbit in order to determine the poloidal transit time, before
   solving the actual particle orbit.*

   **Note:** The :ref:`module-radiation` module expects this parameter to be
   set to ``poloidal``, and the :option:`@ParticlePusher time` parameter to be
   set to ``1``.

