.. _module-distribution-unit:

(unit)
------
The unit distribution function is one (``1``) everywhere in momentum space. If
no radial profile is specified, it is also one at every radius.

Summary of options
^^^^^^^^^^^^^^^^^^
+----------------------------------+------------------------------------------------------------------------------+
| **Option**                       | **Description**                                                              |
+----------------------------------+------------------------------------------------------------------------------+
| :option:`unit radprof`           | Name of configuration block defining the radial profile to use.              |
+----------------------------------+------------------------------------------------------------------------------+

Example configuration
^^^^^^^^^^^^^^^^^^^^^
The following example defines a unit distribution function and uses it for the
simulation::

   distribution_function = ourDistribution;

   @DistributionFunction ourDistribution (unit) {}

An alternative method to specify a unit distribution function is to set the
``distribution_function`` option to the pre-defined
``__unit_distribution_function__``::

   distribution_function = __unit_distribution_function__;

If a unit momentum distribution with a varying radial distribution function is
desired, this can be specified as follows::

   distribution_function = ourDistribution2;

   @DistributionFunction ourDistribution2 (unit) {
       radprof = ourRadialProfile;
   }

   @RadialProfile ourRadialProfile (type) {
       ...
   }

(see the page about :ref:`module-radialprofile` for details about how to
configure the radial profile).

All options
^^^^^^^^^^^
The unit distribution function has no options.

.. program:: unit

.. option:: radprof

   :Default value: Uniform radial profile
   :Allowed values: Name of any defined :ref:`module-radialprofile`

   Specifies the radial profile object to use to generate a radial profile.

