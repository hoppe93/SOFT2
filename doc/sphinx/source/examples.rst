Examples
========
This page collects a number of example SOFT simulations that you can use to
learn how to run SOFT. The example pages show you the expected result and gives
the necessary parameters to reproduce the result. You can then try to write a
SOFT configuration script that properly reproduces the image, spectrum or
figure. Each example comes with a complete solution as well.

All examples, along with plotting scripts, can be found in the ``examples/``
subdirectory of the `SOFT2 <https://github.com/hoppe93/SOFT2>`_ repository.

Orbits
------

.. toctree::
   :hidden:

   examples/Orbits

.. container:: floatblock

   .. container:: leftside

      .. figure:: _static/figs/examples/orbits-particle.png
         :height: 200px
         :target: examples/Orbits.html

   .. container:: rightside

      **Orbits**

      This example illustrates how guiding-center and particle orbits can be
      computed using the :ref:`module-orbits` module of SOFT.

      **Simulation time:** Short

      Check out: :ref:`example-orbits`.

Radiation images
----------------

.. toctree::
   :hidden:

   examples/Zhou2014
   examples/AngularDistribution

.. container:: floatblock

   .. container:: leftside

      .. figure:: _static/figs/examples/Zhou2014.png
         :height: 200px
         :target: examples/Zhou2014.html

   .. container:: rightside

      **Synchrotron images from Zhou et al. (2014)**

      An example showing how to reproduce Figs. 5, 6 and 7 of
      `Zhou et al., PoP 21 (2014) <https://doi.org/10.1063/1.4881469>`_.

      **Simulation time:** Short

      Check out: :ref:`example-zhou2014`.


.. container:: floatblock

   .. container:: leftside

      .. figure:: _static/figs/examples/AngularDistribution.png
         :height: 200px
         :target: examples/AngularDistribution.html

   .. container:: rightside
      
      **Full angular distribution of radiation**

      Run SOFT, taking the full angular distribution of radiation into account.

      **Simulation time:** Medium

      Check out :ref:`example-angulardistribution`.


.. container:: floatblock

   .. container:: leftside

      .. figure:: _static/figs/examples/Bremsstrahlung.png
         :height: 200px
         :target: https://www.google.se/

   .. container:: rightside

      **Bremsstrahlung images**

      Generate bremsstrahlung images

      **Simulation time:** Short

Special output
--------------

Spectrometer output
^^^^^^^^^^^^^^^^^^^

Green's function output
^^^^^^^^^^^^^^^^^^^^^^^

Advanced input
--------------

Numeric magnetic field
^^^^^^^^^^^^^^^^^^^^^^

Distribution function
^^^^^^^^^^^^^^^^^^^^^

