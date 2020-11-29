
.. _module-radialprofile-gaussian:

(gaussian)
**********
This module provides a Gaussian shaped radial profile. The radial profile is
defined according to

.. math::

   f_r(r) = a \exp\left[ -\left(\frac{\rho-b}{c}\right)^2 \right]

This radial profile has been observed to match well with hollow runaway beams.

.. note::

   This module was kindly contributed by
   `Giorgio Ghillardi <https://github.com/giorgioghillardi>`_.

Summary of options
^^^^^^^^^^^^^^^^^^
The following parameters can be set on a gaussian radial profile.

+-----------------------------------+-----------------------------------------------------------------+
| **Option**                        | **Description**                                                 |
+-----------------------------------+-----------------------------------------------------------------+
| :option:`radprof-gaussian amax`   | Normalized minor radius of the outer edge of the electron beam. |
+-----------------------------------+-----------------------------------------------------------------+
| :option:`radprof-gaussian amin`   | Normalized minor radius of the inner edge of the electron beam. |
+-----------------------------------+-----------------------------------------------------------------+
| :option:`radprof-gaussian rhomax` | Major radius of the outer edge of the electron beam.            |
+-----------------------------------+-----------------------------------------------------------------+
| :option:`radprof-gaussian rhomin` | Major radius of the inner edge of the electron beam.            |
+-----------------------------------+-----------------------------------------------------------------+
| :option:`radprof-gaussian rmax`   | Minor radius of the outer edge of the electron beam.            |
+-----------------------------------+-----------------------------------------------------------------+
| :option:`radprof-gaussian rmin`   | Minor radius of the inner edge of the electron beam.            |
+-----------------------------------+-----------------------------------------------------------------+
| :option:`radprof-gaussian gau_a`  | Gaussian scale factor.                                          |
+-----------------------------------+-----------------------------------------------------------------+
| :option:`radprof-gaussian gau_b`  | Gaussian shift parameter.                                       |
+-----------------------------------+-----------------------------------------------------------------+
| :option:`radprof-gaussian gau_c`  | Gaussian width parameter.                                       |
+-----------------------------------+-----------------------------------------------------------------+

Example configuration
^^^^^^^^^^^^^^^^^^^^^

This example illustrates how a gaussian radial profile can be defined and used
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

   @RadialProfile ourradprof (gaussian) {
       # Set normalized minor radius
       # These settings creates a beam with radius
       # 2/3 of the device minor radius.
       amin = 0;
       amax = 0.99;

       # Create gaussian profile f(r) = a * exp(-(r-b)^2/c^2)
       gau_a = 0.05;
       gau_b = 1.80;
       gau_c = 0.04;
   }

All options
^^^^^^^^^^^

.. program:: radprof-gaussian

.. option:: amax

.. option:: amin

.. option:: rhomax

.. option:: rhomin

.. option:: rmax

.. option:: rmin

   :Default value: ``amin = 0`` and ``amax = 1``.
   :Allowed values: Any radial position that is inside the plasma and on the outboard side.

   Specifies the inner and outer edges of the electron beam. The prefix (a*,
   r*, rho*) specifies whether the edge is given in normalized minor radius,
   regular minor radius or major radius coordinates.

.. option:: gau_a

   :Default value: Nothing
   :Allowed values: Any real number.

   Scale factor in front of radial profile.

.. option:: gau_b

   :Default value: Nothing
   :Allowed values: Any real number.

   Shift of the peak in radius of the Gaussian profile.

.. option:: gau_c

   :Default value: Nothing
   :Allowed values: Any real number.

   Width of the Gaussian profile.

