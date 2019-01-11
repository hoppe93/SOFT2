.. _module-distribution-unit:

(unit)
------
The unit distribution function is one (``1``) everywhere, and has no options.

Example configuration
^^^^^^^^^^^^^^^^^^^^^
The following example defines a unit distribution function and uses it for the
simulation::

   distribution_function = ourDistribution;

   @DistributionFunction ourDistribution (unit) {}

Note that no options can be specified in the distribution function.

An alternative method to specify a unit distribution function is to set the
``distribution_function`` option to the pre-defined
``__unit_distribution_function__``::

   distribution_function = __unit_distribution_function__;

All options
^^^^^^^^^^^
The unit distribution function has no options.

