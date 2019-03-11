.. _module-ro-topview:

(topview)
*********
The *topview* @RadiationOutput is used to generate top-view maps of where
radiation is emitted from. The options and output variables for this module are
very similar to that of the :ref:`module-ro-image` module and they mainly differ
in how the image data they generate are interpreted.

Summary of options
^^^^^^^^^^^^^^^^^^

+--------------------------------------------+------------------------------------------------------------------+
| **Option**                                 | **Description**                                                  |
+--------------------------------------------+------------------------------------------------------------------+
| :option:`@RadiationOutput(topview) output` | Sets the name of the output file.                                |
+--------------------------------------------+------------------------------------------------------------------+
| :option:`@RadiationOutput(topview) pixels` | Specifies the number of pixels in each direction of the topview. |
+--------------------------------------------+------------------------------------------------------------------+

Top view geometry
^^^^^^^^^^^^^^^^^
The top view geometry is independent of the detector properties. The global
(tokamak) x-axis is always along the horizontal in the image, while the global
(tokamak) y-axis is always along the vertical. The center of the image
corresponds to the axis of symmetry of the tokamak. The physical size of the
image is twice the maximum major radius of the domain used (i.e. either tokamak
wall or magnetic field separatrix).

An illustration of the geometry of the topview image is shown below. The dashed
square indicates the edges of the image and the dotted circle indicates the
location of the inner wall, i.e. the minimum radius of the domain used.

.. image:: ../../_static/figs/topview_ill.svg
   :width: 80%
   :align: center

Example configuration
^^^^^^^^^^^^^^^^^^^^^
The following generates a topview that is 600x600 pixels in size::

   @RadiationOutput ourTopview (topview) {
       output = "ourTopview.h5";
       pixels = 600;
   }

.. note::

   In contrast to :ref:`module-ro-image`, topviews are always square images.

Output file structure
^^^^^^^^^^^^^^^^^^^^^
The output file contains the following variables:

+-----------------------+---------------------------------------------------------+
| **Variable**          | **Description**                                         |
+-----------------------+---------------------------------------------------------+
| ``detectorDirection`` | Unit vector representing viewing direction of detector. |
+-----------------------+---------------------------------------------------------+
| ``detectorPosition``  | Vector representing position of detector.               |
+-----------------------+---------------------------------------------------------+
| ``detectorVisang``    | (Full) FOV vision angle of the detector.                |
+-----------------------+---------------------------------------------------------+
| ``image``             | Radiation image matrix.                                 |
+-----------------------+---------------------------------------------------------+
| ``wall``              | Domain contour used for the simulation.                 |
+-----------------------+---------------------------------------------------------+

Metadata
--------
Four variables containing a form of metadata are always present in the ouput
file. These are the ``detectorDirection``, ``detectorPosition``,
``detectorVisang`` and ``wall`` variables which give the detector viewing
direction, position, vision angle and the domain contour used in the simulation.
Note that the vision angle is given for the field-of-view, and is twice the
value given as input for backwards-compatibility reasons.

All options
^^^^^^^^^^^

.. program:: @RadiationOutput(topview)

.. option:: output

   :Default value: Nothing
   :Allowed values: Any valid file name.

   Specifies the name of the output file to generate. The file name extension
   determines the type of the output file.

.. option:: pixels

   :Default value: Nothing
   :Allowed values: Any positive integer.

   Specifies the number of pixels in each direction of the image. Thus, the
   total number of pixels in the image will be the square of this number.
   *Note that in contrast to :ref:`module-ro-image`, topviews are always square
   and as such only one number can be assigned to this option.*

