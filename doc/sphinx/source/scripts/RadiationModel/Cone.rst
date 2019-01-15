.. _module-rm-cone:

(cone)
******
The cone radiation model is based on the *cone approximation*, which makes the
assumption that all radiation is emitted *exactly* along the velocity vector of
the particle. This assumption can be very useful for studying the radiation from
highly relativistic particles.

Summary of options
------------------

+--------------------------------------------+---------------------------------------------------------+
| **Option**                                 | **Description**                                         |
+--------------------------------------------+---------------------------------------------------------+
| :option:`@RadiationModel(cone) emission`   | Type of radiation emission to model.                    |
+--------------------------------------------+---------------------------------------------------------+
| :option:`@RadiationModel(cone) projection` | Which projection method to use for the cone model.      |
+--------------------------------------------+---------------------------------------------------------+
| :option:`@RadiationModel(cone) zeff`       | (Bremsstrahlung emission only) Effective plasma charge. |
+--------------------------------------------+---------------------------------------------------------+

Example configuration
---------------------
The following example configures a cone model for bremsstrahlung in a plasma
with effective charge :math:`Z_{\rm eff} = 4`::

   @RadiationModel ourModel (cone) {
       emission = bremsstrahlung;
       zeff     = 4;
   }

All options
-----------

.. program:: @RadiationModel(cone)

.. option:: emission

   :Default value: None
   :Allowed values: ``bremsstrahlung``, ``synchrotron`` or ``unit``

   Specifies the type of radiation emission to model. The three available types
   of radiation are ``bremsstrahlung``, ``synchrotron`` and ``unit``. The two
   former model bremsstrahlung and synchrotron respectively, while ``unit`` is a
   type of radiation that is one everywhere. Note that the radiation is assumed
   to be emitted *exactly* along the velocity vector of the particle so that
   ``unit`` is very different from the :ref:module-rm-isotropic` radiation
   model. The ``unit`` emission type rather has more similarities with the
   ``bremsstrahlung`` emission type.

.. option:: projection

   :Default value: ``reverse``
   :Allowed values: ``reverse`` (or ``original``)

   The cone model relies on a projection of either the guiding-center cone onto
   the detector plane, or the detector surface onto the guiding-center velocity
   orthogonal plane. The ``original`` projection does the former, while the
   ``reverse`` projection does the latter.

   Note that the ``original`` model at the moment contains a bug and should be
   avoided.

.. option:: zeff

   :Default value: ``1``
   :Allowed values: Any real value.

   Specifies the effective plasma charge to use when calculating the emitted
   bremsstrahlung.

