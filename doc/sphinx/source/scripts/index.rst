Configuration Scripts
=====================
The SOFT2 configuration scripting language differs in many ways from that of
the original version of SOFT. In this chapter you will find information about
the basic structure of SOFT2 ``pi`` file (aka configuration script) as well
as a list of available modules and options that can be configured.

Scripting language
------------------
SOFT uses a custom language for configuring runs. The language is used for
assigning values to configuration parameters, and as such, the most common
construct found in a SOFT configuration file is::

   parameterName = value;

This assigns the value ``value`` to the parameter named ``parameterName``. Note
the semi-colon at the end of the line, which is mandatory.

Parameters exist on many levels---they can apply globally or to specific modules
used in the simulation---and are grouped using curly brackets, ``{`` and ``}``.
Everything between a pair of curly brackets is usually referred to as a *block*.
Assignments that are outside of a block apply globally, i.e. they do not belong
to any specific module, but are rather used by *all* modules.

Module parameters are assigned using the syntax shown below::

   @ModuleName blockName {
       parameterName1 = value1;
       ...
   }

First we specify which type of module we're configuring by stating its name,
prefixed with an ``@`` sign: ``@ModuleName``. A module is a component of the
simulation that provides specific functionality. Examples include
``@ParticleGenerator`` which generates particles for our simulations,
``@ParticlePusher`` which simulates guiding-center and particle orbits, and
``@Radiation`` which simulates a radiation detector.

The ``blockName`` is specific to this block of settings. It is possible to
define several configuration blocks for the same module ``@ModuleName``. Which
one(s) to use is specified by a parameter in the parent module. For example,
the ``@Equation`` module (which defines which equations of motion to solve) is
used by the ``@ParticlePusher`` module, and hence we could define the equation
to solve using the following code::

   @ParticlePusher PPusher {
       ...
       equation = gcEquation;
       ...
   }
   @EquationGuidingCenter gcEquation {
       tolerance = 1e-9;
   }
   @EquationParticle pEquation {
       tolerance = 1e-9;
   }

Here we define two separate equations but we only use one. This is due to that
we assign the configuration block ``gcEquation`` to the ``equation`` parameter
in the ``PPusher`` configuration block. The ``@ParticlePusher`` is an example of
a module that provides core functionality in SOFT and is therefore selected
using a global parameter.

Comments
********
SOFT configuration files support single-line comments. They are started with the
``#`` symbol and causes SOFT to ignore everything that appears until the next
newline character.

Include other files
*******************
In SOFT2, it is possible to source other configuration scripts. This is done by
placing the name of the file to include between ``<`` and ``>``. The
included configuration file will be loaded and all its settings applied before
proceeding in the parent configuration file. This means that whichever parameter
assignment comes last applies. As an example, consider the following two
configuration files:

**file1**::

   parameter1 = value1;
   parameter2 = value2;

**file2**::

   parameter1 = value3;
   
   <file1>

   parameter2 = value4;

After passing ``file2`` to SOFT, i.e. running ``soft file2``, then ``parameter1``
will be ``value1``, whereas ``parameter2`` (which was assigned after including
``file1``) will be ``value4``.

Multiple values and strings
***************************
SOFT2 allows multiple values to be assigned to a single parameter. An example
of when this is useful is for vector quantities, such as detector position,
which needs one value for each component of its position vector. SOFT achieves
this by splitting the assigned values at each whitespace and comma. For example,
the following lines assigns three different values to the parameter::

   parameterName = value1 value2 value3;
   parameterName = value1,value2,value3;

(Note: the two lines are shown to show off the possibilities in SOFT; the two
lines are equivalent). In some cases it may be desired to assign a value
containing a space to parameter though, such as when assigning the name of a
file. To overcome this, SOFT allows values to be surrounded by double quotes::

   parameterName = "value1 value2 value3";

The double quotes causes ``parameterName`` to now instead be assigned *a single*
value which contains spaces. Double quotes can be used everywhere to emphasize
what the value being assigned is (even surrounding numbers).

FAQ
***
Some frequently asked questions:

**Are newlines mandatory?**
No, newline characters are not mandatory in SOFT configuration files. The entire
configuration could be written on a single line. They are however very
convenient and highly recommended!

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
   EquationGuidingCenter
   EquationParticle
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
| :ref:`module-equation-gc`       | Guiding-center equations of motion                  |
+---------------------------------+-----------------------------------------------------+
| :ref:`module-equation-particle` | Particle equations of motion                        |
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

