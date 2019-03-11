.. _module-radialprofile-uniform:

(uniform)
*********
This module provides a radial profile that is one (1) everywhere.

Summary of options
^^^^^^^^^^^^^^^^^^
This module has no options.

Example configuration
^^^^^^^^^^^^^^^^^^^^^

This example illustrates how a uniform radial profile can be defined and used
together with the :ref:`module-distribution-avalanche` distribution function in
SOFT2 simulations::

   # Global parameter specifying which distribution function to use
   distribution_function = analytical_avalanche;

   @DistributionFunction analytical_avalanche (avalanche) {
       EHat     = 10;            # Electric field strength (normalized to
                                 # the Connor-Hastie critical electric field)
       lnLambda = 17;            # Coulomb logarithm
       Zeff      = 4;             # Effective plasma charge

       # Set the radial profile
       radprof   = ourradprof;
   }

   @RadialProfile ourradprof (uniform) {}

Note that if the option ``radprof`` of the distribution function is not
specified, a uniform radial profile is automatically generated. Thus, the
following configuration yields the same result as the above::

   # Global parameter specifying which distribution function to use
   distribution_function = analytical_avalanche;

   @DistributionFunction analytical_avalanche (avalanche) {
       EHat     = 10;            # Electric field strength (normalized to
                                 # the Connor-Hastie critical electric field)
       lnLambda = 17;            # Coulomb logarithm
       Zeff      = 4;             # Effective plasma charge
   }

