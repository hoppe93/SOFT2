.. _module-magneticfield:

@MagneticField
**************
The magnetic field module defines a magnetic field and domain (walls) to use
in the simulation. For setting up a quick simulation, the analytical magnetic
field can be used to generate a magnetic field with circular flux surfaces.
For more advanced applications, SOFT2 is able to take numeric 2D magnetic
fields as input.

More details about the magnetic field in SOFT2 can be found on the page
:ref:`magnetic-fields`.

Summary of options
^^^^^^^^^^^^^^^^^^

+-------------------------------+----------------------------------------+--------------------------------------------------------------------------------+
| **Option**                    | **Type**                               | **Description**                                                                |
+-------------------------------+----------------------------------------+--------------------------------------------------------------------------------+
| :option:`@MagneticField type` | :ref:`module-magneticfield-analytical` | Type of magnetic field: ``analytical`` or ``numeric``                          |
+-------------------------------+----------------------------------------+--------------------------------------------------------------------------------+
| :option:`analytical B0`       | :ref:`module-magneticfield-analytical` | Magnetic field strength (on-axis) of analytical magnetic field                 |
+-------------------------------+----------------------------------------+--------------------------------------------------------------------------------+
| :option:`analytical Rm`       | :ref:`module-magneticfield-analytical` | Major radius of analytical magnetic field                                      |
+-------------------------------+----------------------------------------+--------------------------------------------------------------------------------+
| :option:`analytical rminor`   | :ref:`module-magneticfield-analytical` | Minor radius of analytical magnetic field                                      |
+-------------------------------+----------------------------------------+--------------------------------------------------------------------------------+
| :option:`analytical sigmaB`   | :ref:`module-magneticfield-analytical` | Toroidal field direction of analytical magnetic field                          |
+-------------------------------+----------------------------------------+--------------------------------------------------------------------------------+
| :option:`analytical sigmaI`   | :ref:`module-magneticfield-analytical` | Current direction of analytical magnetic field                                 |
+-------------------------------+----------------------------------------+--------------------------------------------------------------------------------+
| :option:`analytical qtype`    | :ref:`module-magneticfield-analytical` | Safety factor type: ``constant``, ``linear``, ``qudratic`` or ``exponential``  |
+-------------------------------+----------------------------------------+--------------------------------------------------------------------------------+
| :option:`analytical qa1`      | :ref:`module-magneticfield-analytical` | First safety factor parameter of analytical magnetic field                     |
+-------------------------------+----------------------------------------+--------------------------------------------------------------------------------+
| :option:`analytical qa2`      | :ref:`module-magneticfield-analytical` | Second safety factor parameter of analytical magnetic field                    |
+-------------------------------+----------------------------------------+--------------------------------------------------------------------------------+
| :option:`numeric filename`    | :ref:`module-magneticfield-numeric`    | Name of file containing numeric magnetic field                                 |
+-------------------------------+----------------------------------------+--------------------------------------------------------------------------------+
| :option:`numeric filetype`    | :ref:`module-magneticfield-numeric`    | Override filetype of file containing numeric magnetic field                    |
+-------------------------------+----------------------------------------+--------------------------------------------------------------------------------+

Example configuration
^^^^^^^^^^^^^^^^^^^^^
The following defines an analytical magnetic field, similar to Alcator C-Mod::

   magnetic_field = analytical-field;

   @MagneticField analytical-field {
       type   = analytical;
       B0     = 5.0;
       Rm     = 0.68;
       rminor = 0.22;
       qtype  = quadratic;
       qa1    = 3.0;
       qa2    = 1.0;
   }

An example of a numeric magnetic field::

   magnetic_field = numeric-field;

   @MagneticField numeric-field {
       type     = numeric;
       filename = "/path/to/magnetic-field.mat";
       # SOFT will automatically identify the
       # above file as a 'MAT' file, due to its
       # '.mat' filename extension.
   }

Common settings
^^^^^^^^^^^^^^^

.. program:: @MagneticField

.. option:: type

   :Default value: None
   :Allowed values: ``analytical`` or ``numeric``

   Specifies the type of magnetic field to define. There are two types of
   magnetic fields in SOFT2: analytical and numerical. The former are defined
   through an equation presented on the page :ref:`magnetic-fields`. The numeric
   magnetic field is loaded from file which must have certain variables in it.
   The details about the layout of magnetic field files in SOFT can be found
   on the page about :ref:`magnetic-fields`.

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

   +-------------+------------------------------------+
   | **qtype**   | **Functional form**                |
   +-------------+------------------------------------+
   | constant    | :math:`q(a) = q_{a1}`              |
   +-------------+------------------------------------+
   | linear      | :math:`q(a) = q_{a1} a + q_{a2}`   |
   +-------------+------------------------------------+
   | quadratic   | :math:`q(a) = q_{a1} a^2 + q_{a2}` |
   +-------------+------------------------------------+
   | exponential | :math:`q(a) = e^{q_{a1}} + q_{a2}` |
   +-------------+------------------------------------+

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

