.. _module-distribution-avalanche:

(avalanche)
***********
In a spatially uniform fully ionized plasma with constant electric field, when
the runaway generation is dominated by the avalanche mechanism---i.e. by
multiplication through large-angle collisions---an analytical expression can be
found for the resulting runaway electron distribution function. This
distribution function, which was first derived in [#fulop2006]_, takes the form

.. math::

   f(p, \xi) = \frac{A(p)}{2\pi m_e c\gamma_0 p^2}
   \frac{\exp\left[ -\gamma/\gamma_0 - A(p)(1+\xi) \right]}{1 - \exp\left[-2A(p)\right]}

with

.. math::

   A(p) &= \frac{\hat{E} + 1}{Z_{\rm eff} + 1} \gamma,\\
   \gamma_0 &= \log\Lambda\sqrt{Z_{\rm eff}+5}

where :math:`p` is the electron momentum, :math:`\xi = \cos\theta_{\rm p}` the
electron pitch with respect to the magnetic field, :math:`\gamma = \sqrt{1 + p^2}`,
:math:`\hat{E}` the electric field strength normalized to the threshold electric
field, :math:`Z_{\rm eff}` the plasma effective charge, :math:`\log\Lambda` the
plasma Coloumb logarithm, :math:`m_e` the electron rest mass and :math:`c` the
speed of light in vacuum.

.. [#fulop2006] Fülöp et al., 2006 "Destabilization of magnetosonic-whistler waves by a relativistic runaway beam". *Physics of Plasmas* **13** (6), 062506 `doi:10.1063/1.2208327 <https://doi.org/10.1063/1.2208327>`_.

Summary of options
^^^^^^^^^^^^^^^^^^
The following parameters can be set on an avalanche distribution function.

+----------------------------------+----------------------------------------------------------------------+
| **Option**                       | **Description**                                                      |
+----------------------------------+----------------------------------------------------------------------+
| :option:`avalanche EHat`         | Electric field strength, normalized to the threshold electric field. |
+----------------------------------+----------------------------------------------------------------------+
| :option:`avalanche lnLambda`     | Coloumb logarithm of the plasma.                                     |
+----------------------------------+----------------------------------------------------------------------+
| :option:`avalanche Zeff`         | Plasma effective charge.                                             |
+----------------------------------+----------------------------------------------------------------------+
| :option:`avalanche radprof`      | Name of configuration block for radial profile.                      |
+----------------------------------+----------------------------------------------------------------------+

Example configuration
^^^^^^^^^^^^^^^^^^^^^

This example illustrates how an analytical avalanche distribution function can
be defined and used in SOFT2 simulations::

   # Global parameter specifying which distribution function to use
   distribution_function = analytical_avalanche;

   @DistributionFunction analytical_avalanche (avalanche) {
       EHat     = 10;            # Electric field strength (normalized to
                                 # the Connor-Hastie critical electric field)
       lnLambda = 17;            # Coulomb logarithm
       Zeff      = 4;             # Effective plasma charge
   }

All options
^^^^^^^^^^^

.. program:: avalanche

.. option:: EHat

   :Default value: None
   :Allowed values: Any real number

   Electric field strength, normalized to the Connor-Hastie critical electric field.

.. option:: lnLambda

   :Default value: None
   :Allowed values: Any real numer

   Value of the Coulomb logarithm.

.. option:: radprof

   :Default value: Uniform radial profile
   :Allowed values: Name of any defined :ref:`module-radialprofile`

   Specifies the radial profile object to use to generate a radial profile.

.. option:: Zeff

   :Default value: None
   :Allowed values: Any real number

   Effective plasma charge.

