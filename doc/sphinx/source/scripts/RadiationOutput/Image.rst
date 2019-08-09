.. _module-ro-image:

(image)
*******
The *image* radiation output generates radiation pixel images. The value of each
pixel is determined by the amount of radiation within a (usually very small)
pyramid-like field-of-view extending from the detector. A new feature of SOFT2
is that images can be rectangular, and not just square as in the old version of
SOFT.

Summary of options
^^^^^^^^^^^^^^^^^^

+------------------------------------------------+------------------------------------------------------------------------------------------------+
| **Option**                                     | **Description**                                                                                |
+------------------------------------------------+------------------------------------------------------------------------------------------------+
| :option:`@RadiationOutput(image) common`       | List of common quantities to include in the output file.                                       |
+------------------------------------------------+------------------------------------------------------------------------------------------------+
| :option:`@RadiationOutput(image) output`       | Sets the name of the output file.                                                              |
+------------------------------------------------+------------------------------------------------------------------------------------------------+
| :option:`@RadiationOutput(image) pixels`       | Specifies the number of pixels in each direction of the image.                                 |
+------------------------------------------------+------------------------------------------------------------------------------------------------+
| :option:`@RadiationOutput(image) stokesparams` | If ``yes``, measures polarization and stores on image per Stokes parameter (synchrotron only). |
+------------------------------------------------+------------------------------------------------------------------------------------------------+

Detector plane basis
^^^^^^^^^^^^^^^^^^^^
The detector plane is spanned by two vectors which we refer to as
:math:`\hat{e}_1` and :math:`\hat{e}_2`. Together with the detector viewing
direction (or normal vector) :math:`\hat{n}` they form an orthonormal basis in
space. The vectors :math:`\hat{e}_1` and :math:`\hat{e}_2` are used for several
purposes, but one of the more important purposes is for spanning the camera
images that this @RadiationOutput produces.

The vector :math:`\hat{e}_1` is defined so that it always lies in the horizontal
plane. Its mathematical definition is

.. math::

   \hat{e}_1 = \begin{cases} \hat{y}\cos\vartheta + \hat{z}\sin\vartheta, \quad&\text{ if } \hat{n}\cdot\hat{y} = 0,\\
   \left[ \hat{x}\left(\hat{n}\cdot\hat{y}\right) - \hat{y}\left( \hat{n}\cdot\hat{x} \right)\right]\cos\vartheta + \hat{z}\sin\vartheta,
   \quad&\text{ otherwise}.
   \end{cases}

where :math:`\vartheta` denotes the angle between :math:`\hat{e}_1` and the
horizontal plane (the "tilt angle"). From this, :math:`\hat{e}_2` is defined
so that :math:`(\hat{e}_1, \hat{e}_2, \hat{n})` form a right-handed orthonormal
basis:

.. math::

   \hat{e}_2 = \hat{n}\times\hat{e}_1.

The values of both of these vectors are output on ``stdout`` by SOFT during
execution.

Example configuration
^^^^^^^^^^^^^^^^^^^^^
The following generates a square radiation image that is 600x600 pixels in
size::

   @RadiationOutput ourImage (image) {
       output = "ourImage.h5";
       pixels = 600;
   }

In SOFT2, images can also be rectangular. The following example generates a
rectangular image with 600 pixels along the :math:`\hat{e}_1` axis (always in
the horizontal plane), and 450 pixels along the
:math:`\hat{e}_2` axis::

   @RadiationOutput ourRectangularImage (image) {
       output = "ourRectangularImage.h5";
       pixels = 600 450;
   }

One can also measure the polarization of incoming radiation. This feature is
only available for spectral-dependent synchrotron radiation. To enable
measurement of Stokes parameters, simply enabled the ``stokesparams`` option::

   @RadiationOutput ourStokesImage (image) {
       output       = "ourStokesImage.h5";
       pixels       = 600 450;
       stokesparams = yes;
   }

If ``stokesparams`` is ``yes``, then the output file will contain four variables
called ``StokesI``, ``StokesQ``, ``StokesU`` and ``StokesV``, instead of the
usual ``image`` variable. Each of the ``Stokes*`` variables is an image with
each pixel corresponding to the value of the Stokes parameters in its
field-of-view.

Output file structure
^^^^^^^^^^^^^^^^^^^^^
The output file contains the following variables:

+-----------------------+----------------------+---------------------------------------------------------+
| **Variable**          | **Included when**    | **Description**                                         |
+-----------------------+----------------------+---------------------------------------------------------+
| ``image``             | ``stokesparams=no``  | Radiation image matrix.                                 |
+-----------------------+----------------------+---------------------------------------------------------+
| ``StokesI``           | ``stokesparams=yes`` | Image of Stokes :math:`I` parameter.                    |
+-----------------------+----------------------+---------------------------------------------------------+
| ``StokesQ``           | ``stokesparams=yes`` | Image of Stokes :math:`Q` parameter.                    |
+-----------------------+----------------------+---------------------------------------------------------+
| ``StokesU``           | ``stokesparams=yes`` | Image of Stokes :math:`U` parameter.                    |
+-----------------------+----------------------+---------------------------------------------------------+
| ``StokesV``           | ``stokesparams=yes`` | Image of Stokes :math:`V` parameter.                    |
+-----------------------+----------------------+---------------------------------------------------------+

Common quantities
-----------------
By default, the following "common quantities" are also included in the output
file:

+-----------------------+---------------------------------------------------------+
| **Name**              | **Description**                                         |
+-----------------------+---------------------------------------------------------+
| ``detectorDirection`` | Unit vector representing viewing direction of detector. |
+-----------------------+---------------------------------------------------------+
| ``detectorPosition``  | Vector representing position of detector.               |
+-----------------------+---------------------------------------------------------+
| ``detectorVisang``    | (Full) FOV vision angle of the detector.                |
+-----------------------+---------------------------------------------------------+
| ``wall``              | Domain contour used for the simulation.                 |
+-----------------------+---------------------------------------------------------+

*For details about which other common quantities can be included in the output,
please consult the page about the* :ref:`module-radiationoutput` *class of
modules.*

.. note::

   The actual image is contained *either* in the ``image`` variable if the input
   parameter ``stokesparams=no``, or in the ``StokesI``, ``StokesQ``, ``StokesU``
   and ``StokesQ`` variables if ``stokesparams=yes``. These variables are
   matrices of the same size as the number of pixels specified in the input file.
   Note that the first dimension corresponds to :math:`\hat{e}_1`, which always
   lies in the horizontal plane, meaning the one often desires to transpose the
   image before showing it.

.. tip::

   If you are uncertain about the direction of the image, you can try to move
   the camera vertically upwards or downwards. When doing so, you should expect
   the radiation spot to move in the opposite direction in the image.

All options
^^^^^^^^^^^

.. program:: @RadiationOutput(image)

.. option:: common

   :Default value: ``none``
   :Allowed values: See the list on :ref:`module-radiationoutput`.

   Specifies which "common quantities" to include in the output file. A full
   list of possible options is given on :ref:`module-radiationoutput`.

.. option:: output

   :Default value: Nothing
   :Allowed values: Any valid file name.

   Specifies the name of the output file to generate. The file name extension
   determines the type of the output file.

.. option:: pixels

   :Default value: Nothing
   :Allowed values: One or two positive integers.

   Specifies the number of pixels in the image along the :math:`\hat{e}_1` and
   :math:`\hat{e}_2` directions respectively. Either one or two numbers can be
   given. If only one number is given, the number of pixels will be the same
   along both axes (and equal to the given number), yielding a square image. If
   two numbers are specified, the first number gives the number of pixels along
   the :math:`\hat{e}_1` axis and the second number gives the number of pixels
   along the :math:`\hat{e}_2` axis.

.. option:: stokesparams

   :Default value: ``no``
   :Allowed values: ``yes`` or ``no``

   If ``yes``, measures the polarization of the radiation and produces one
   image for each of the four Stokes parameters :math:`I`, :math:`Q`, :math:`U`
   and :math:`V`. This feature is only available for synchrotron radiation and
   requires the detector measure in a limited spectral range.

