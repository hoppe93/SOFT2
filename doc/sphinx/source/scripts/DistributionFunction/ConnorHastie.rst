.. _module-distribution-connor:

(connorhastie)
**************
When the runaway generation is dominated by the primary (or Dreicer) mechanism,
the distribution function is predicted to take the form given by
Connor & Hastie [#connor1975]_:

.. math::

   f(p, \xi) = \frac{1}{p\xi}\exp\left[ -\frac{\hat{E}+1}{2(Z_{\rm eff}+1)}\frac{1-\xi^2}{\xi}p \right],

where :math:`p` is the electron momentum, :math:`\xi = \cos\theta_{\rm p}` the
electron pitch with respect to the magnetic field, :math:`\hat{E}` the electric
field strength normalized to the threshold electric field, and :math:`Z_{\rm eff}`
the plasma effective charge. The threshold electric field :math:`E_{\rm c}` is given by

.. math::

   E_{\rm c} = \frac{n_ee^3\ln\Lambda}{4\pi\epsilon_0^2 m_ec^2},

where :math:`n_e` is the electron density, :math:`e` the elementary charge,
:math:`\ln\lambda` the Coulomb logarithm, :math:`\epsilon_0` the permittivity of
free space, :math:`m_e` the electron mass and :math:`c` the speed of light in vacuum.

Since the synchrotron radiation emitted by these electrons primarily comes from
those electrons with the highest energy, the distribution above can be combined
with a Gaussian-shaped envelope for :math:`p>p_{\rm max}`, for some
:math:`p_{\rm max}`, that cuts off the distribution function at high energies:

.. math::

   f(p, \xi) \to f(p,\xi) \exp\left[ - \left(\frac{p - p_{\rm max}}{\Delta p}\right)^2 \right].

Note that this cut-off is only applied for :math:`p\geq p_{\rm max}`, and so for
lower values of momentum the standard Connor-Hastie distribution is used.

.. [#connor1975] Connor and Hastie, 1975 "Relativistic limitations on runaway electrons". *Nuclear Fusion* **15** (3), 415 `doi:10.1088/0029-5515/15/3/007 <https://doi.org/10.1088/0029-5515/15/3/007>`_.

Summary of options
^^^^^^^^^^^^^^^^^^
The following parameters can be set on a Connor-Hastie distribution function.

+----------------------------------+----------------------------------------------------------------------+
| **Option**                       | **Description**                                                      |
+----------------------------------+----------------------------------------------------------------------+
| :option:`connorhastie deltaP`    | Width of the Gaussian envelope to apply to cut-off the distribution. |
+----------------------------------+----------------------------------------------------------------------+
| :option:`connorhastie pMax`      | Momentum value above which to apply the Gaussian cut-off.            |
+----------------------------------+----------------------------------------------------------------------+
| :option:`connorhastie EHat`      | Electric field strength, normalized to the threshold electric field. |
+----------------------------------+----------------------------------------------------------------------+
| :option:`connorhastie Zeff`      | Plasma effective charge.                                             |
+----------------------------------+----------------------------------------------------------------------+
| :option:`connorhastie radprof`   | Name of configuration block for radial profile.                      |
+----------------------------------+----------------------------------------------------------------------+

Example configuration
^^^^^^^^^^^^^^^^^^^^^

This example illustrates how a Connor-Hastie distribution function can
be defined and used in SOFT2 simulations::

   # Global parameter specifying which distribution function to use
   distribution_function = connor-hastie;

   @DistributionFunction connor-hastie (connorhastie) {
       EHat     = 10;            # Electric field strength (normalized to
                                 # the Connor-Hastie critical electric field)
       Zeff      = 4;            # Effective plasma charge
   }

All options
^^^^^^^^^^^

.. program:: connorhastie

.. option:: deltaP

   :Default value: None
   :Allowed values: Any positive real number

   Width of the Gaussian envelope to apply for :math:`p > p_{\rm max}` in order
   to cut off the distribution function at high energies.

.. option:: pMax

   :Default value: None
   :Allowed values: Any positive real number

   Momentum above which the Connor-Hastie distribution function will be modified
   by a Gaussian envelope function. The purpose of this is to cut-off the
   distribution at high energies, so that ``pMax`` effectively corresponds to
   the "maximum" momentum value in the distribution function.

.. option:: EHat

   :Default value: None
   :Allowed values: Any real number

   Electric field strength, normalized to the Connor-Hastie critical electric field
   (see above).

.. option:: radprof

   :Default value: Uniform radial profile
   :Allowed values: Name of any defined :ref:`module-radialprofile`

   Specifies the radial profile object to use to generate a radial profile.

.. option:: Zeff

   :Default value: None
   :Allowed values: Any real number

   Effective plasma charge.

