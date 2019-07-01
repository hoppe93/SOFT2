.. _module-distribution-pitch:

(pitch)
-------
The pitch distribution function is an exponential in the cosine of the pitch
angle:

.. math::

   f(\xi) = \exp\left( C\xi \right),

where :math:`\xi = \cos\theta_{\rm p}` and :math:`C` is an input parameter.
This is the pitch part of the analytical avalanche distribution function
(see :ref:`module-distribution-avalanche`).

Summary of options
^^^^^^^^^^^^^^^^^^
+-------------------------+-----------------------------------------------------------------+
| **Option**              | **Description**                                                 |
+-------------------------+-----------------------------------------------------------------+
| :option:`pitch C`       | Distribution function exponent.                                 |
+-------------------------+-----------------------------------------------------------------+
| :option:`pitch radprof` | Name of configuration block defining the radial profile to use. |
+-------------------------+-----------------------------------------------------------------+

Example configuration
^^^^^^^^^^^^^^^^^^^^^
The following example defines a pitch distribution function and uses it for the
simulation::

   distribution_function = ourDistribution;

   @DistributionFunction ourDistribution (pitch) {
       C = 100;
   }

If a pitch distribution with a varying radial distribution function is desired,
this can be specified as follows::

   distribution_function = ourDistribution2;

   @DistributionFunction ourDistribution2 (pitch) {
       radprof = ourRadialProfile;
   }

   @RadialProfile ourRadialProfile (type) {
       ...
   }

(see the page about :ref:`module-radialprofile` for details about how to
configure the radial profile).

All options
^^^^^^^^^^^
The pitch distribution function has no options.

.. program:: pitch

.. option:: C

   :Default value: None
   :Allowed values: Any positive real number

   Specifies the value of the exponent in the pitch distribution.

.. option:: radprof

   :Default value: Uniform radial profile
   :Allowed values: Name of any defined :ref:`module-radialprofile`

   Specifies the radial profile object to use to generate a radial profile.

