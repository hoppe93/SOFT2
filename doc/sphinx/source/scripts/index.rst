Configuration Scripts
=====================
The SOFT2 configuration scripting language differs in many ways from that of
the original version of SOFT. In this chapter you will find information about
the basic structure of SOFT2 ``pi`` file (aka configuration script) as well
as a list of available modules and options that can be configured.

Scripting language
------------------

Global options
--------------
.. toctree::
   :hidden:

   global

The following options are *global*, i.e. are set outside of any environment.
Click on the options below for further details.

+------------------------------------------+---------------------------------------------------------+
| **Option**                               | **Description**                                         |
+------------------------------------------+---------------------------------------------------------+
| :option:`global distribution_function`   | Specifies configuration of distribution function to use |
+------------------------------------------+---------------------------------------------------------+
| :option:`global include_drifts`          | Whether or not to include guiding-center drifts         |
+------------------------------------------+---------------------------------------------------------+
| :option:`global magnetic_field`          | Specifies configuration of magnetic field to use        |
+------------------------------------------+---------------------------------------------------------+
| :option:`global num_threads`             | Number of OpenMP threads to run SOFT on                 |
+------------------------------------------+---------------------------------------------------------+
| :option:`global particle_generator`      | Specifies particle generator (phase-space) to use       |
+------------------------------------------+---------------------------------------------------------+
| :option:`global particle_pusher`         | Specifies particle pusher (orbit follower) to use       |
+------------------------------------------+---------------------------------------------------------+

Configurable modules
--------------------

.. toctree::
   :hidden:

   Detector
   DistributionFunction
   Equation
   MagneticField
   Orbits
   ParticleGenerator
   ParticlePusher
   RadiationModel
   RadiationOutput
   RadialProfile
   Radiation

+---------------------------------+-----------------------------------------------------+
| **Module**                      | **Description**                                     |
+---------------------------------+-----------------------------------------------------+
| :ref:`module-distribution`      | Distribution functions                              |
+---------------------------------+-----------------------------------------------------+
| :ref:`module-equation`          | Equation to solve: "particle" or "guiding-center"   |
+---------------------------------+-----------------------------------------------------+
| :ref:`module-magneticfield`     | Magnetic field module                               |
+---------------------------------+-----------------------------------------------------+
| :ref:`module-particlegenerator` | Configuration of phase-space                        |
+---------------------------------+-----------------------------------------------------+
| :ref:`module-particlepusher`    | Particle pusher and time-stepping                   |
+---------------------------------+-----------------------------------------------------+
| :ref:`module-radialprofile`     | Radial profile                                      |
+---------------------------------+-----------------------------------------------------+

+---------------------------------+-----------------------------------------------------+
| **Tool**                        | **Description**                                     |
+---------------------------------+-----------------------------------------------------+
| :ref:`module-orbits`            | Simulate and output particle/guiding-center orbits  |
+---------------------------------+-----------------------------------------------------+
| :ref:`module-radiation`         | Synthetic radiation diagnostic tool                 |
+---------------------------------+-----------------------------------------------------+

+---------------------------------+-----------------------------------------------------+
| **Radiation sub-module**        | **Description**                                     |
+---------------------------------+-----------------------------------------------------+
| :ref:`module-detector`          | Synthetic detector configuration                    |
+---------------------------------+-----------------------------------------------------+
| :ref:`module-radiationmodel`    | Radiation model (cone/angular distribution)         |
+---------------------------------+-----------------------------------------------------+
| :ref:`module-radiationoutput`   | Radiation output module                             |
+---------------------------------+-----------------------------------------------------+

