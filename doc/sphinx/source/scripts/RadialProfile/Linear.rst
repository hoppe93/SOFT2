
.. _module-radialprofile-linear:

(linear)
********
This module provides a radial profile the decays linearly from a radius
:math:`r_{\rm min}` to a radius :math:`r_{\rm max}` and is zero outside of that
radial interval. More precisely, on the interval
:math:`r_{\rm min}\leq r \leq r_{\rm max}` the radial profile is

.. math::

   f_r(r) = 1 - \frac{r - r_{\rm min}}{r_{\rm max} - r_{\rm min}},

and it is identically zero everywhere else.

Summary of options
^^^^^^^^^^^^^^^^^^
The following parameters can be set on a linear radial profile.

+---------------------------------+-----------------------------------------------------------------+
| **Option**                      | **Description**                                                 |
+---------------------------------+-----------------------------------------------------------------+
| :option:`radprof-linear amax`   | Normalized minor radius of the outer edge of the electron beam. |
+---------------------------------+-----------------------------------------------------------------+
| :option:`radprof-linear amin`   | Normalized minor radius of the inner edge of the electron beam. |
+---------------------------------+-----------------------------------------------------------------+
| :option:`radprof-linear rhomax` | Major radius of the outer edge of the electron beam.            |
+---------------------------------+-----------------------------------------------------------------+
| :option:`radprof-linear rhomin` | Major radius of the inner edge of the electron beam.            |
+---------------------------------+-----------------------------------------------------------------+
| :option:`radprof-linear rmax`   | Minor radius of the outer edge of the electron beam.            |
+---------------------------------+-----------------------------------------------------------------+
| :option:`radprof-linear rmin`   | Minor radius of the inner edge of the electron beam.            |
+---------------------------------+-----------------------------------------------------------------+

Example configuration
^^^^^^^^^^^^^^^^^^^^^

This example illustrates how a linear radial profile can be defined and used
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

   @RadialProfile ourradprof (linear) {
       # Set normalized minor radius
       # These settings creates a beam with radius
       # 2/3 of the device minor radius.
       amin = 0;
       amax = 0.67;
   }

All options
^^^^^^^^^^^

.. program:: radprof-linear

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

