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

