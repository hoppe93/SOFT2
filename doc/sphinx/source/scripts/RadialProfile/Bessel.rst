
.. _module-radialprofile-bessel:

(bessel)
********
This module provides a radial profile that decays according to a zeroth-order
Bessel function of the first kind, from a radius
:math:`r_{\rm min}` to a radius :math:`r_{\rm max}` and is zero outside of that
radial interval. More precisely, on the interval
:math:`r_{\rm min}\leq r \leq r_{\rm max}` the radial profile is

.. math::

   f_r(r) = J_0\left( \frac{r-r_{\rm min}}{r_{\rm max}-r_{\rm min}} x_0 \right),

and it is identically zero everywhere else. Here, :math:`J_0(x)` denotes a
zeroth-order Bessel function of the first kind and :math:`x_0` is its first
zero.

This type of radial profile is expected if runaway electrons are radially 
transported diffusively in the plasma, with a uniform diffusion coefficient.

Summary of options
^^^^^^^^^^^^^^^^^^
The following parameters can be set on a bessel radial profile.

+---------------------------------+-----------------------------------------------------------------+
| **Option**                      | **Description**                                                 |
+---------------------------------+-----------------------------------------------------------------+
| :option:`radprof-bessel amax`   | Normalized minor radius of the outer edge of the electron beam. |
+---------------------------------+-----------------------------------------------------------------+
| :option:`radprof-bessel amin`   | Normalized minor radius of the inner edge of the electron beam. |
+---------------------------------+-----------------------------------------------------------------+
| :option:`radprof-bessel rhomax` | Major radius of the outer edge of the electron beam.            |
+---------------------------------+-----------------------------------------------------------------+
| :option:`radprof-bessel rhomin` | Major radius of the inner edge of the electron beam.            |
+---------------------------------+-----------------------------------------------------------------+
| :option:`radprof-bessel rmax`   | Minor radius of the outer edge of the electron beam.            |
+---------------------------------+-----------------------------------------------------------------+
| :option:`radprof-bessel rmin`   | Minor radius of the inner edge of the electron beam.            |
+---------------------------------+-----------------------------------------------------------------+

Example configuration
^^^^^^^^^^^^^^^^^^^^^

This example illustrates how a bessel radial profile can be defined and used
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

   @RadialProfile ourradprof (bessel) {
       # Set normalized minor radius
       # These settings creates a beam with radius
       # 2/3 of the device minor radius.
       amin = 0;
       amax = 0.67;
   }

All options
^^^^^^^^^^^

.. program:: radprof-bessel

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

