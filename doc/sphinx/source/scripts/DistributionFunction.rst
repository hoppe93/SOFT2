.. _module-distribution:

@DistributionFunction
*********************
The distribution function module, which defines a distribution function to
weigh results with. A number of different types of distribution functions can be
generated or loaded in SOFT2, including the analytical avalanche distribution [#fulop2006]_,
direct output from the kinetic solver CODE, as well as regular SOFT phase-space
distribution functions.

.. [#fulop2006] Fülöp et al., 2006 "Destabilization of magnetosonic-whistler waves by a relativistic runaway beam". *Physics of Plasmas* **13** (6), 062506 `doi:10.1063/1.2208327 <https://doi.org/10.1063/1.2208327>`_.

Summary of options
^^^^^^^^^^^^^^^^^^

+----------------------------------+-----------------------------------------------------------------------+
| **Option**                       | **Distribution function type**                                        |
+----------------------------------+-----------------------------------------------------------------------+
| :option:`@Distribution type`     | All                                                                   |
+----------------------------------+-----------------------------------------------------------------------+
| :option:`avalanche EHat`         | :ref:`module-distribution-avalanche`                                  |
+----------------------------------+-----------------------------------------------------------------------+
| :option:`avalanche lnLambda`     | :ref:`module-distribution-avalanche`                                  |
+----------------------------------+-----------------------------------------------------------------------+
| :option:`avalanche Zeff`         | :ref:`module-distribution-avalanche`                                  |
+----------------------------------+-----------------------------------------------------------------------+
| :option:`avalanche radprof`      | :ref:`module-distribution-avalanche`, :ref:`module-distribution-code` |
+----------------------------------+-----------------------------------------------------------------------+
| :option:`code name`              | :ref:`module-distribution-code`, :ref:`module-distribution-numerical` |
+----------------------------------+-----------------------------------------------------------------------+
| :option:`code interptype`        | :ref:`module-distribution-code`, :ref:`module-distribution-numerical` |
+----------------------------------+-----------------------------------------------------------------------+
| :option:`code time`              | :ref:`module-distribution-code`                                       |
+----------------------------------+-----------------------------------------------------------------------+
| :option:`numerical logarithmize` | :ref:`module-distribution-numerical`                                  |
+----------------------------------+-----------------------------------------------------------------------+

Example configuration
^^^^^^^^^^^^^^^^^^^^^

This example illustrates how an analytical avalanche distribution function can
be defined and used in SOFT2 simulations::

   @DistributionFunction analytical_avalanche {
       type     = "avalanche";

       # Parameters specific to the analytical avalanche function
       EHat     = 10;            # Electric field strength (normalized to
                                 # the Connor-Hastie critical electric field)
       lnLambda = 17;            # Coulomb logarithm
       Zeff     = 4;             # Effective plasma charge
   }

Common settings
^^^^^^^^^^^^^^^

.. program:: @Distribution

.. option:: type

   :Default value: None
   :Allowed values: ``avalanche``, ``code`` or ``numerical``

   Specifies the type of the distribution function. Depending on which type is
   specified here, different options can be set for the distribution function,
   as detailed below for the different types.

.. _module-distribution-avalanche:

avalanche
^^^^^^^^^

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

.. _module-distribution-code:

code
^^^^

.. program:: code

.. option:: interptype

   :Default value: ``cspline``
   :Allowed values: ``akima``, ``akima_periodic``, ``cspline``, ``cspline_periodic``, ``linear``, ``polynomial``, ``steffen``

   Determines which interpolation method to use for interpolating in the
   momentum dimension.

.. option:: name

   :Default value: None
   :Allowed values: String

   Specifies the name of the file containing the CODE distribution function
   to load.

.. option:: radprof

   :Default value: Uniform radial profile
   :Allowed values: Name of any defined :ref:`module-radialprofile`

   Specifies the radial profile object to use to generate a radial profile.

.. option:: time

   :Default value: ``-1`` (last timestep)
   :Allowed values: Any integer with absolute value less than the number of time points in the distribution function

   Selects the index of the time step to take the distribution function from.
   Negative indices count from the back of the array, so that ``-1`` corresponds
   to the last timestep, ``-2`` to the next-to-last etc.


.. _module-distribution-numerical:

numerical
^^^^^^^^^

.. program:: numerical

.. option:: interptype

   :Default value: ``cubic``
   :Allowed values: ``cubic``, ``linear``

   Determines which interpolation method to use for interpolating in
   momentum-space. Support is available for bilinear and bicubic splines.

.. option:: logarithmize

   :Default value: ``no``
   :Allowed values: ``yes``, ``no``

   If ``yes``, SOFT will logarithmize the distribution function before
   interpolation. Thus, interpolation is done in :math:`\log f` rather than
   the distribution function itself.

.. option:: name

   :Default value: None
   :Allowed values: String

   Specifies the name of the file containing the CODE distribution function
   to load.


.. _module-distribution-unit:

unit
^^^^

The unit distribution function is one (``1``) everywhere, and has no options.

