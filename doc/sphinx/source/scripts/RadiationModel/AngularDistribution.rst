.. _module-rm-angulardistribution:

(angdist)
*********
The angular distribution radiation model takes the full angular distribution
of radiation into account. In contrast to the cone model---which assumes that
all radiation is emitted *exactly* along the particle velocity vector---this
model allows for arbitrary shapes of the angular distribution of radiation.

Summary of options
------------------
The following options are available in the ``angdist`` radiation model.

+---------------------------------------------+-------------------------------------------------------------------+
| **Option**                                  | **Description**                                                   |
+---------------------------------------------+-------------------------------------------------------------------+
| :option:`@RadiationModel(angdist) emission` | Type of radiation emission to model.                              |
+---------------------------------------------+-------------------------------------------------------------------+
| :option:`@RadiationModel(angdist) nsamples` | Number of points in each dimension of the detector surface.       |
+---------------------------------------------+-------------------------------------------------------------------+
| :option:`@RadiationModel(angdist) qrule2d`  | Quadrature rule to use for integrating over the detector surface. |
+---------------------------------------------+-------------------------------------------------------------------+

Example configuration
---------------------
The following example configures a synchrotron radiation model with the detector
surface discretized using a total of 16 points::

   @RadiationModel ourModel (angdist) {
       emission = synchrotron;
       nsamples = 4;
   }

All options
-----------

.. program:: @RadiationModel(angdist)

.. option:: emission

   :Default value: None
   :Allowed values: ``bremsstrahlung`` or ``synchrotron``

   Specifies the type of radiation emission to model. Currently, the two
   available options are ``bremsstrahlung`` and ``synchrotron``. The
   appropriate formulas to use are chosen automatically, depending on whether
   spectral dependence is considered, and if guiding-center drifts are included
   or not.

.. option:: nsamples

   :Default value: 1
   :Allowed values: Any positive integer.

   Number of points in each dimension to discretize the detector surface with.
   The total number of points on the detector surface is therefore
   ``nsamples^2``.

.. option:: qrule2d

   :Default value: ``simpson``
   :Allowed values: ``simpson``

   Specifies the quadrature rule to use for integrating over the detector
   surface. Currently, only Simpson's rule has been implemented, with the
   exception of the special case ``nsamples = 1``, in which case the function
   is merely evaluated in the detector center-point.

