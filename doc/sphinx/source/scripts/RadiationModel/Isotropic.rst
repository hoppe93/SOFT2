.. _module-rm-isotropic:

(isotropic)
***********
The isotropic radiation model assumes that all particles emit radiation
uniformly in all directions. Mathematically, this radiation model takes the form

.. math::
   :label: eq-rm-isotropic

   \frac{\mathrm{d} P}{\mathrm{d}\Omega} = P_0,

where :math:`P_0` is a user-defined constant. Note that the angular distribution
of radiation is independent of any particle or background plasma properties.

Summary of options
------------------

+--------------------------------------------+----------------------------------------------------+
| **Option**                                 | **Description**                                    |
+--------------------------------------------+----------------------------------------------------+
| :option:`@RadiationModel(isotropic) value` | The value of :math:`P_0` in :eq:`eq-rm-isotropic`. |
+--------------------------------------------+----------------------------------------------------+

Example configuration
---------------------
The following example configures an isotropic model::

   @RadiationModel ourModel (isotropic) {
       value = 1;
   }

All options
-----------

.. program:: @RadiationModel(isotropic)

.. option:: value

   :Default value: 1.
   :Allowed values: Any real value.

   Power per unit solid angle emitted by the particle.

