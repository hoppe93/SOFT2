.. _module-radiationoutput:

@RadiationOutput
****************
The ``@RadiationOutput`` modules allow one to consider various radiation
quantities, such as images, spectra, Green's functions and more. At the moment,
there are five different radiation output sub-modules available:
:ref:`module-ro-green`, :ref:`module-ro-image`, :ref:`module-ro-sovvolume`,
:ref:`module-ro-space3d`, :ref:`module-ro-spectrum` and :ref:`module-ro-topview`.

.. toctree::
   :hidden:

   RadiationOutput/Green
   RadiationOutput/Image
   RadiationOutput/SoVVolume
   RadiationOutput/Space3D
   RadiationOutput/Spectrum
   RadiationOutput/Topview

Available sub-modules
---------------------

+---------------------------------+---------------------------------------------+
| **Module name**                 | **Output description**                      |
+---------------------------------+---------------------------------------------+
| :ref:`module-ro-green`          | Green's/weight functions                    |
+---------------------------------+---------------------------------------------+
| :ref:`module-ro-image`          | Camera images                               |
+---------------------------------+---------------------------------------------+
| :ref:`module-ro-sovvolume`      | Calculate the surface-of-visibility volume. |
+---------------------------------+---------------------------------------------+
| :ref:`module-ro-space3d`        | 3D maps of radiation                        |
+---------------------------------+---------------------------------------------+
| :ref:`module-ro-spectrum`       | Radiation spectra                           |
+---------------------------------+---------------------------------------------+
| :ref:`module-ro-topview`        | Tokamak topviews of radiation               |
+---------------------------------+---------------------------------------------+

Include common data
--------------------
All ``@RadiationOutput`` modules support the ``common`` option, which allows you
to specify which common/meta data to include in the output file. The value
assigned to the option should be a list of names of the data objects to include
in the file. The available options are shown in the table below.

+------------------------------+-------------------------------------------------------------+
| **Name**                     | **Description**                                             |
+------------------------------+-------------------------------------------------------------+
| ``detectorAperture``         | Detector size/aperture.                                     |
+------------------------------+-------------------------------------------------------------+
| ``detectorDirection``        | Detector viewing direction.                                 |
+------------------------------+-------------------------------------------------------------+
| ``detectorPosition``         | Detector position (relative to point-of-symmetry).          |
+------------------------------+-------------------------------------------------------------+
| ``detectorVisang``           | Detector vision angle (FOV).                                |
+------------------------------+-------------------------------------------------------------+
| ``distribution``             | Distribution function, as evaluated by SOFT.                |
+------------------------------+-------------------------------------------------------------+
| ``domain``                   | Orbit solution domain. Either tokamak wall or separatrix.   |
+------------------------------+-------------------------------------------------------------+
| ``param1``                   | First momentum grid                                         |
+------------------------------+-------------------------------------------------------------+
| ``param1name``               | Name of first momentum grid parameter.                      |
+------------------------------+-------------------------------------------------------------+
| ``param2``                   | Second momentum grid                                        |
+------------------------------+-------------------------------------------------------------+
| ``param2name``               | Name of second momentum grid parameter.                     |
+------------------------------+-------------------------------------------------------------+
| ``r``                        | Radial grid (major radius).                                 |
+------------------------------+-------------------------------------------------------------+
| ``wall``                     | Alias for ``domain``.                                       |
+------------------------------+-------------------------------------------------------------+

In addition to these, it is also possible to specify any of the options in the
table below. Those options do not however represent single objects, such as
those in the table above, but instead enables or disables bulks of objects to
output.

+-------------+----------------------------------------------+
| **Name**    | **Description**                              |
+-------------+----------------------------------------------+
| ``all``     | Include all available objects in the output. |
+-------------+----------------------------------------------+
| ``default`` | Include the default objects in the output.   |
+-------------+----------------------------------------------+
| ``none``    | Do not include any common output.            |
+-------------+----------------------------------------------+

Additionally, it is possible to prefix any of the options in the first table
with either a ``+`` or ``-`` to indicate whether the object should be included
(``+``) or not (``-``). Options specified later overrides former ones, meaning
that the line::

   common = all -domain -wall

would include all the available objects in the output, except for the ``domain``
and ``wall`` objects. Similarly,

::

   common = default none +domain

would first enable all default objects, then undo that option and remove all
objects, and finally add the ``domain``. Thus, the only common object in the
output would be the ``domain`` object.

To see which common quantities are included with the ``default`` option, please
consult the page of the relevant ``@RadiationOutput`` module.

Example configuration
---------------------
Please, see the pages for each sub-module for examples of how to configure each
of them.

Note that several ``@RadiationOutput`` modules can be used with each 
:ref:`module-radiation` module. Just specify the name of each output module to
use, separated by spaces or commas, in the :ref:`module-radiation` configuration
block::

   @Radiation rad {
       ...
       output = outModule1 outModule2 outModule3;
       # ...or...
       output = outModule1,outModule2,outModule3;
       # ...or even...
       output = outModule1, outModule2, outModule3;
   }

To specify which common quantities to include, the option ``common`` should be
included in the appropriate ``@RadiationOutput`` configuration block::

   @RadiationOutput outModule1 (XXX) {
       ...
       common = default +domain;
       ...
   }

