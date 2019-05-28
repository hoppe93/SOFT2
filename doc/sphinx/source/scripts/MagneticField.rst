.. _module-magneticfield:

@MagneticField
**************
The magnetic field module defines a magnetic field and domain (walls) to use
in the simulation. For setting up a quick simulation, the analytical magnetic
field can be used to generate a magnetic field with circular flux surfaces.
For more advanced applications, SOFT2 is able to take numeric 2D magnetic
fields as input.

The type of the magnetic field to load is specified as a secondary type of the
block (i.e. in parenthesis after the block name). At the moment, two different
types of magnetic fields are available: ``analytical`` and ``numeric``.

More details about the magnetic field in SOFT2 can be found on the page
:ref:`magnetic-fields`.

Relating safety factor and current profile
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The safety factor type ``current`` has been derived from the current density

.. math::

   \boldsymbol{j} = \sigma_I\hat{\boldsymbol{\varphi}}j(r)
   = \sigma_I\hat{\boldsymbol{\varphi}} j_0\left[ 1 - \left(\frac{r}{r_0}\right)^{q_{a2}} \right],

where :math:`j_0` is the central current density and :math:`r_0` is the plasma
radius (see :ref:`magnetic-fields` for details about the notation). The
corresponding safety factor is

.. math::
   
   q(r) = \frac{q_0}{1 - \frac{2}{n+2}\left( \frac{r}{r_0} \right)^{q_{a2}}}.

From this, the plasma current :math:`I` can be related to the central safety
factor :math:`q_0\equiv q(0)`:

.. math::
    
    I = \frac{2\pi B_0}{\mu_0 q_0 R_{\rm m}} r_0^2 \left( 1 - \frac{2}{n+2} \right),

or conversely

.. math::
   
   q_0 = \frac{2\pi B_0}{\mu_0 R_{\rm m} I} r_0^2 \left( 1 - \frac{2}{n+2} \right),

where :math:`\mu_0\approx 4\pi\times 10^{-7}\,\text{H/m}` is the vacuum
permeability. With :math:`B_0`, :math:`\mu_0`, :math:`R_{\rm m}` and :math:`r_0`
in SI units, and the total plasma current in mega-ampere, the central safety
factor can be computed from

.. math::
   
   q_0 = \frac{5n}{n+2} \frac{B_0 r_0^2}{R_{\rm m} I}.

Summary of options
^^^^^^^^^^^^^^^^^^

+-------------------------------+----------------------------------------+---------------------------------------------------------------------------------------------+
| **Option**                    | **Type**                               | **Description**                                                                             |
+-------------------------------+----------------------------------------+---------------------------------------------------------------------------------------------+
| :option:`analytical B0`       | :ref:`module-magneticfield-analytical` | Magnetic field strength (on-axis) of analytical magnetic field                              |
+-------------------------------+----------------------------------------+---------------------------------------------------------------------------------------------+
| :option:`analytical Rm`       | :ref:`module-magneticfield-analytical` | Major radius of analytical magnetic field                                                   |
+-------------------------------+----------------------------------------+---------------------------------------------------------------------------------------------+
| :option:`analytical zm`       | :ref:`module-magneticfield-analytical` | Vertical position of magnetic axis                                                          |
+-------------------------------+----------------------------------------+---------------------------------------------------------------------------------------------+
| :option:`analytical rminor`   | :ref:`module-magneticfield-analytical` | Minor radius of analytical magnetic field                                                   |
+-------------------------------+----------------------------------------+---------------------------------------------------------------------------------------------+
| :option:`analytical sigmaB`   | :ref:`module-magneticfield-analytical` | Toroidal field direction of analytical magnetic field                                       |
+-------------------------------+----------------------------------------+---------------------------------------------------------------------------------------------+
| :option:`analytical sigmaI`   | :ref:`module-magneticfield-analytical` | Current direction of analytical magnetic field                                              |
+-------------------------------+----------------------------------------+---------------------------------------------------------------------------------------------+
| :option:`analytical qtype`    | :ref:`module-magneticfield-analytical` | Safety factor type: ``current``, ``constant``, ``linear``, ``qudratic`` or ``exponential``  |
+-------------------------------+----------------------------------------+---------------------------------------------------------------------------------------------+
| :option:`analytical qa1`      | :ref:`module-magneticfield-analytical` | First safety factor parameter of analytical magnetic field                                  |
+-------------------------------+----------------------------------------+---------------------------------------------------------------------------------------------+
| :option:`analytical qa2`      | :ref:`module-magneticfield-analytical` | Second safety factor parameter of analytical magnetic field                                 |
+-------------------------------+----------------------------------------+---------------------------------------------------------------------------------------------+
| :option:`numeric filename`    | :ref:`module-magneticfield-numeric`    | Name of file containing numeric magnetic field                                              |
+-------------------------------+----------------------------------------+---------------------------------------------------------------------------------------------+
| :option:`numeric filetype`    | :ref:`module-magneticfield-numeric`    | Override filetype of file containing numeric magnetic field                                 |
+-------------------------------+----------------------------------------+---------------------------------------------------------------------------------------------+

Example configuration
^^^^^^^^^^^^^^^^^^^^^
The following defines an analytical magnetic field, similar to Alcator C-Mod::

   magnetic_field = analytical-field;

   @MagneticField analytical-field (analytical) {
       B0     = 5.0;
       Rm     = 0.68;
       rminor = 0.22;
       qtype  = quadratic;
       qa1    = 3.0;
       qa2    = 1.0;
   }

An example of a numeric magnetic field::

   magnetic_field = numeric-field;

   @MagneticField numeric-field (numeric) {
       filename = "/path/to/magnetic-field.mat";
       # SOFT will automatically identify the
       # above file as a 'MAT' file, due to its
       # '.mat' filename extension.
   }

Common settings
^^^^^^^^^^^^^^^

.. program:: @MagneticField

.. _module-magneticfield-analytical:

analytical
^^^^^^^^^^

.. program:: analytical

.. option:: B0

   :Default value: None
   :Allowed values: Any positive real number

   The magnetic field strength on-axis (in Tesla).

.. option:: Rm

   :Default value: None
   :Allowed values: Any positive real number (greater than :option:`analytical rminor`)

   The tokamak major radius (in meters).

.. option:: zm

   :Default value: ``0``
   :Allowed values: Any real number

   The vertical position of the magnetic axis.

.. option:: rminor

   :Default value: None
   :Allowed values: Any positive real number (less than :option:`analytical Rm`)

   The tokamak minor radius (in meters).

.. option:: sigmaB

.. option:: sigmaI

   :Default value: ``cw``
   :Allowed values: ``cw`` / ``-``, or ``ccw`` / ``+``

   Sign of the toroidal field component (``sigmaB``) and plasma current
   (``sigmaI``). The value ``cw`` corresponds to the toroidal component being
   oriented in the clock-wise direction, when seen from above the tokamak, while
   ``ccw`` corresponds to the toroidal component being oriented in the counter
   clock-wise direction, when seen from above.

   Instead of specifying the direction, the sign may be given directly as either
   ``-`` (clock-wise) or ``+`` (counter clock-wise).

.. option:: qa1

.. option:: qa2

   :Default value: 1.0
   :Allowed values: Any real number

   Constants used to define the safety factor. See :option:`analytical qtype`
   for details about how exactly they affect the different safety factor.

.. option:: qtype

   :Default value: None
   :Allowed values: ``constant``, ``exponential``, ``linear``, ``quadratic``

   Specifies the radial dependence of the safety factor. The functional forms
   for the various options are given in terms of the normalized minor radius
   :math:`a` (normalized to the value of ``rminor``) in the table below. The
   constants :math:`q_{a1}` and :math:`q_{a2}` are specified using the
   :option:`analytical qa1` and :option:`analytical qa2` options.

   +-------------+----------------------------------------------------------------+
   | **qtype**   | **Functional form**                                            |
   +-------------+----------------------------------------------------------------+
   | current     | :math:`q(a) = q_{a1}/\left( 1-\frac{2}{n+2}a^{q_{a2}} \right)` |
   +-------------+----------------------------------------------------------------+
   | constant    | :math:`q(a) = q_{a1}`                                          |
   +-------------+----------------------------------------------------------------+
   | linear      | :math:`q(a) = q_{a1} a + q_{a2}`                               |
   +-------------+----------------------------------------------------------------+
   | quadratic   | :math:`q(a) = q_{a1} a^2 + q_{a2}`                             |
   +-------------+----------------------------------------------------------------+
   | exponential | :math:`q(a) = e^{q_{a1}} + q_{a2}`                             |
   +-------------+----------------------------------------------------------------+

.. _module-magneticfield-numeric:

numeric
^^^^^^^

.. program:: numeric

.. option:: filename

   :Default value: None
   :Allowed values: String

   Specifies the name of the file that contains the magnetic field to load.

.. option:: filetype

   :Default value: Auto-determine filetype based on filename
   :Allowed values: ``h5``, ``hdf5``, ``mat``, ``out``, ``sdt``

   Overrides the filetype of the given file. Usually SOFT tries to determine
   which filetype a given file is of based on its filename extension. By
   explicitly setting this option, this check is overriden allows you to use
   non-standard filename extensions for the input file.

