.. _module-radiationmodel:

@RadiationModel
***************
This module determines how the *angular dependence* of the emitted radiation
should be treated. We can conveniently express the position of an observer
relative to the emitting particle using a spherical coordinate system
:math:`(r, \theta, \phi)`. Since the *direction* in which the observer is
located relative to the emitter is fully determined by the two angles
:math:`\theta` and :math:`\phi`, and since the amount of radiation *emitted*
towards an observer is independent of how far away the observer is located,
we can speak of the *angular distribution* of emitted radiation.

To not assume anything about the angular distribution of the emitted radiation,
the :ref:`module-rm-angulardistribution` should be used. It accounts for the
full angular distribution of all types of radiation, and is therefore also
somewhat heavier to run.

Directed radiation, such as bremsstrahlung and synchrotron radiation, are
almost exactly emitted along the velocity vector of the emitting particle, with
an angular spread that is proportional to :math:`\gamma^{-1}`, where
:math:`\gamma` denotes the relativistic factor. As an approximation, we may
therefore assume that all radiation is emitted *exactly* along the velocity
vector of the particle. We refer to this approximation as the *cone
approximation*, and it is implemented in the :ref:`module-rm-cone` module.
It corresponds to assuming

.. math::

   \frac{\mathrm{d} P}{\mathrm{d}\Omega} = \frac{P_0}{2\pi}\delta\left( \hat{\boldsymbol{v}}\cdot\hat{\boldsymbol{n}} - 1 \right),

where :math:`\frac{\mathrm{d} P}{\mathrm{d}\Omega}` is the angular distribution
of radiation (see the synthetic radiation diagnostic equation on the page about
:ref:`module-radiation`), :math:`P_0` is the total amount of emitted radiation
(possibly in a limited wavelength range), :math:`\delta` is a Dirac delta
function, :math:`\hat{\boldsymbol{v}}` denotes the direction of motion of the
particle and :math:`\hat{\boldsymbol{n}}` the unit vector pointing out the
direction along which the particle is being observed.

The cone model is significantly faster to run than the
:ref:`module-rm-angulardistribution` model, and can produce accurate results if
the energies of the emitting particles are sufficiently high.

Finally, we may also assume that the emitted radiation is independent of the
two angles :math:`\theta` and :math:`\phi` and that the radiation is emitted
uniformly in all directions. A special radiation model has been implemented
for this case, called :ref:`module-rm-isotropic`, which is primarily used for
benchmarking purposes and producing pretty images.

.. toctree::
   :hidden:

   RadiationModel/AngularDistribution
   RadiationModel/Cone
   RadiationModel/Isotropic

Available models
----------------

+--------------------------------------+--------------------------------------------------------------------------+
| **Model**                            | **Description**                                                          |
+--------------------------------------+--------------------------------------------------------------------------+
| :ref:`module-rm-angulardistribution` | No approximations. Consider full angular distribution of radiation.      |
+--------------------------------------+--------------------------------------------------------------------------+
| :ref:`module-rm-cone`                | Cone approximation. All radiation emitted exactly along velocity vector. |
+--------------------------------------+--------------------------------------------------------------------------+
| :ref:`module-rm-isotropic`           | Isotropically emitted radiation.                                         |
+--------------------------------------+--------------------------------------------------------------------------+

