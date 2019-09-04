.. _module-rm-cone:

(cone)
******
The cone radiation model is based on the *cone approximation*, which makes the
assumption that all radiation is emitted *exactly* along the velocity vector of
the particle. This assumption can be very useful for studying the radiation from
highly relativistic particles.

Summary of options
------------------

+--------------------------------------------+--------------------------------------------------------------------------+
| **Option**                                 | **Description**                                                          |
+--------------------------------------------+--------------------------------------------------------------------------+
| :option:`@RadiationModel(cone) emission`   | Type of radiation emission to model.                                     |
+--------------------------------------------+--------------------------------------------------------------------------+
| :option:`@RadiationModel(cone) n`          | (Bremsstrahlung emission only) List of plasma species densities.         |
+--------------------------------------------+--------------------------------------------------------------------------+
| :option:`@RadiationModel(cone) projection` | Which projection method to use for the cone model.                       |
+--------------------------------------------+--------------------------------------------------------------------------+
| :option:`@RadiationModel(cone) Z`          | (Bremsstrahlung emission only) List of plasma species atomic numbers.    |
+--------------------------------------------+--------------------------------------------------------------------------+
| :option:`@RadiationModel(cone) Z0`         | (Bremsstrahlung emission only) List of plasma species net charges.       |
+--------------------------------------------+--------------------------------------------------------------------------+

Example configuration
---------------------
The following example configures a cone model for bremsstrahlung in a
deuterium-tritium-tungsten plasma::

   @RadiationModel ourModel (cone) {
       emission = bremsstrahlung;
       Z        = 1, 1, 74;
       Z0       = 1, 1, 72;
       n        = 5e19;
   }

All options
-----------

.. program:: @RadiationModel(cone)

.. option:: emission

   :Default value: None
   :Allowed values: ``bremsstrahlung``, ``bremsstrahlung_screened``, ``synchrotron`` or ``unit``

   Specifies the type of radiation emission to model. The four available types
   of radiation are ``bremsstrahlung``, ``bremsstrahlung_screened``,
   ``synchrotron`` and ``unit``. The first two model bremsstrahlung, the latter
   taking the effect of screening of the nuclei into account. The third model
   synchrotron radiation, while ``unit`` is a type of radiation that is one
   everywhere. Note that the radiation is assumed to be emitted *exactly* along
   the velocity vector of the particle so that ``unit`` is very different from
   the :ref:module-rm-isotropic` radiation model. The ``unit`` emission type
   rather has more similarities with the ``bremsstrahlung`` emission type.

.. option:: n

   :Default value: None
   :Allowed values: A list of real values.

   Specifies the density to use for each plasma ion species when calculating the
   emitted electron-ion bremsstrahlung.

.. option:: projection

   :Default value: ``reverse``
   :Allowed values: ``reverse`` (or ``original``)

   The cone model relies on a projection of either the guiding-center cone onto
   the detector plane, or the detector surface onto the guiding-center velocity
   orthogonal plane. The ``original`` projection does the former, while the
   ``reverse`` projection does the latter.

   Note that the ``original`` model at the time of writing contains a bug and
   should be avoided.

.. option:: Z

   :Default value: None
   :Allowed values: A list of real values.

   Specifies the effective plasma charge to use for each plasma ion species when
   calculating the emitted electron-ion bremsstrahlung.

.. option:: Z0

   :Default value: None
   :Allowed values: A list of real values.

   Specifies the plasma net charge to use for each plasma ion species when
   calculating the emitted electron-ion bremsstrahlung.

