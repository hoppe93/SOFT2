.. _module-detector:

@Detector
**********
The ``@Detector`` module is a submodule of the :ref:`module-radiation` tool.
This module configures the radiation detector to use, and sets properties such
as its position, viewing direction, field-of-view, spectral range and more.

When specifying detector position, the following illustration of the coordinate
system assumed in SOFT can be useful:

.. image:: ../figs/detector_coordinates.svg
   :width: 70%
   :align: center

Summary of options
^^^^^^^^^^^^^^^^^^

+------------------------+------------------------------------------------------------+
| **Option**             | **Description**                                            |
+------------------------+------------------------------------------------------------+
| :option:`aperture`     | Aperture size / square side length                         |
+------------------------+------------------------------------------------------------+
| :option:`direction`    | Viewing direction / detector plane normal vector           |
+------------------------+------------------------------------------------------------+
| :option:`position`     | Position relative to tokamak point-of-symmetry             |
+------------------------+------------------------------------------------------------+
| :option:`vision_angle` | Detector vision angle (field-of-view (half) opening angle) |
+------------------------+------------------------------------------------------------+
| :option:`spectrum`     | Spectral range configuration of detector                   |
+------------------------+------------------------------------------------------------+

Example configuration
^^^^^^^^^^^^^^^^^^^^^
This example configuration of a detector corresponds to wide-angle visible camera
located in the midplane of a medium-sized tokamak::

   @Detector example_detector {
       aperture     = 0.006;              # in meters
       direction    = 0, 1, 0;            # x,y,z
       position     = 0, 1.7, 0;          # x,y,z (relative to point of symmetry)
       vision_angle = 1.25 fov;           # Field-of-view (half) opening angle
       spectrum     = 440e-9, 790e-9, 40; # Lower wavelength (m), Upper (m), Number of points
   }

Options
^^^^^^^

.. option:: aperture

   | **Default value:** None
   | **Allowed values:** Any positive real number

   Size of detector aperture. All detectors are modeled as squares in SOFT, with
   the aperture specified here corresponding to the square side length.

.. option:: direction

   | **Default value:** None
   | **Example line:** ``direction = 1.3, -0.25, 0;``
   | **Allowed values:** Any real 3-vector except null

   Detector viewing direction, i.e. normal vector of the detector plane. This
   vector is normalized internally by SOFT to become a unit vector, and does not
   have to specified as a unit vector.

.. option:: position

   | **Default value:** None
   | **Example line:** ``position = 0, -1.069, 0;``
   | **Allowed values:** ``

   Detector position relative to the tokamak point-of-symmetry. Units used for
   the vector components are meters.

.. option:: vision_angle

   | **Default value:** None
   | **Example line:** ``vision_angle = 0.75 image;``
   | **Allowed values:** Real number, optionally followed by either ``fov`` or ``image``

   Specifies the half opening angle of the field-of-view. If only a number is
   given, then it is implied that the ``fov`` opening angle is specified.

   A camera image will be a square inscribed in the circular field-of-view. By
   appending ``image`` after the real number, you indicate that the given vision
   angle is the minimum angle between a vector extending from the detector to the
   side of the image, and the detector normal. In a 2D top view, if a cone with
   the specified ``image`` vision angle was plotted, this would correspond exactly
   to the edges of the observed image.

.. option:: spectrum

   | **Default value:** ``no``
   | **Example line:** ``spectrum = 785e-9, 795e-9, 10;``
   | **Allowed values:** (i) Vector of two positive real numbers and one positive integer, or (ii) ``no``

   Sets the limits and resolution of the detector spectral range. If the number
   of spectral points is 0, or if this parameter is assigned the boolean ``no``
   value, the detector will have an infinite spectral range.

   The numbers given to the spectral range specify (i) the lower spectral bound,
   (ii) the upper spectral bound, and (iii) the number of points on the interval.
   The units of the bounds depend on the emission model used. The synchrotron
   emission models will assume that the spectrum limits are *wavelengths* given
   in units of meters. The bremsstrahlung models assume that the spectrum limits
   are photon energies, given in normalized units (normalized to :math:`m_e c^2`,
   the electron mass times speed-of-light squared).

