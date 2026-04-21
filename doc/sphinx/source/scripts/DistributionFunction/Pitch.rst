.. _module-distribution-pitch:

(pitch)
-------
The pitch distribution function is an exponential in the cosine of the pitch
angle:

.. math::

   f(p,\xi) = \exp\left[ C(p) \left(\left|\xi\right|-1\right) \right],

where :math:`\xi = \cos\theta_{\rm p}` and

.. math::

   C(p) = C_0 + C_1p

describes the slope of the distribution function. This is the pitch part of the
analytical avalanche distribution function (see
:ref:`module-distribution-avalanche`).


Summary of options
^^^^^^^^^^^^^^^^^^
+-------------------------+-----------------------------------------------------------------+
| **Option**              | **Description**                                                 |
+-------------------------+-----------------------------------------------------------------+
| :option:`pitch C0`      | Constant term in exponent.                                      |
+-------------------------+-----------------------------------------------------------------+
| :option:`pitch C1`      | Linear term pre-factor in exponent.                             |
+-------------------------+-----------------------------------------------------------------+
| :option:`pitch radprof` | Name of configuration block defining the radial profile to use. |
+-------------------------+-----------------------------------------------------------------+

.. note::

   The default value for ``C1`` is zero, so that the pitch spectrum becomes
   a simple exponential.

Example configuration
^^^^^^^^^^^^^^^^^^^^^
The following example defines a pitch distribution function and uses it for the
simulation::

   distribution_function = ourDistribution;

   @DistributionFunction ourDistribution (pitch) {
       C0 = 100;
   }

A linear variation with momentum can also be obtained by setting::

   distribution_function = ourDistribution;

   @DistributionFunction ourDistribution (pitch) {
       C0 = -573;
       C1 = 16.4;
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

.. option:: C0

   :Default value: None
   :Allowed values: Any positive real number

   Specifies the constant term in the exponent of the pitch distribution.


.. option:: C1

   :Default value: 0
   :Allowed values: Any positive real number

   Specifes the pre-factor for the linear term in the exponent of the pitch
   distribution.

.. option:: radprof

   :Default value: Uniform radial profile
   :Allowed values: Name of any defined :ref:`module-radialprofile`

   Specifies the radial profile object to use to generate a radial profile.

