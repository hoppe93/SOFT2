.. _module-radialprofile-power:

(power)
*******
This module provides a radial profile the decays according to a power-law from
a radius :math:`r_{\rm min}` to a radius :math:`r_{\rm max}` and is zero outside
of that radial interval. More precisely, on the interval
:math:`r_{\rm min}\leq r \leq r_{\rm max}` the radial profile is

.. math::
   :label: eq-radprof-power

   f_r(r) = 1 - \left( \frac{r - r_{\rm min}}{r_{\rm max} - r_{\rm min}} \right)^b,

and it is identically zero everywhere else. The parameter :math:`b` must be
specified by the user.

Summary of options
^^^^^^^^^^^^^^^^^^
The following parameters can be set on a power radial profile.

+----------------------------------+-----------------------------------------------------------------+
| **Option**                       | **Description**                                                 |
+----------------------------------+-----------------------------------------------------------------+
| :option:`radprof-power exponent` | Exponent in radial profile.                                     |
+----------------------------------+-----------------------------------------------------------------+
| :option:`radprof-power amax`     | Normalized minor radius of the outer edge of the electron beam. |
+----------------------------------+-----------------------------------------------------------------+
| :option:`radprof-power amin`     | Normalized minor radius of the inner edge of the electron beam. |
+----------------------------------+-----------------------------------------------------------------+
| :option:`radprof-power rhomax`   | Major radius of the outer edge of the electron beam.            |
+----------------------------------+-----------------------------------------------------------------+
| :option:`radprof-power rhomin`   | Major radius of the inner edge of the electron beam.            |
+----------------------------------+-----------------------------------------------------------------+
| :option:`radprof-power rmax`     | Minor radius of the outer edge of the electron beam.            |
+----------------------------------+-----------------------------------------------------------------+
| :option:`radprof-power rmin`     | Minor radius of the inner edge of the electron beam.            |
+----------------------------------+-----------------------------------------------------------------+

Example configuration
^^^^^^^^^^^^^^^^^^^^^

This example illustrates how a power radial profile can be defined and used
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

   @RadialProfile ourradprof (power) {
       # Set normalized minor radius
       # These settings creates a beam with radius
       # 2/3 of the device minor radius.
       amin = 0;
       amax = 0.67;

       # Create quadratically decreasing profile
       exponent = 2;
   }

All options
^^^^^^^^^^^

.. program:: radprof-power

.. option:: exponent

   :Default value: None (i.e. must be defined)
   :Allowed values: Any real value.

   Specifies the exponent :math:`b` in :eq:`eq-radprof-power`.

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

