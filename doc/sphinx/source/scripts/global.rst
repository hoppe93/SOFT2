.. _options-global:

List of global options
----------------------

.. program:: global

.. option:: distribution_function

   :Default value: A unit (=1) distribution function
   :Allowed values: Name of distribution function configuration object

   Specifies the distribution function to use. Distribution functions are configured
   separately in their own ``@DistributionFunction`` environments. Only one
   distribution function can be set per simulation.

   More details can be found on the page about :ref:`module-distribution`.

   I wonder if we can ref general stuff :ref:`some_option`.

.. option:: include_drifts

   :Default value: ``no``
   :Allowed values: Boolean

   When enabled, keeps terms up to and including order :math:`\epsilon` in
   all expressions. This means that guiding-center orbits drift, and that
   more accurate radiation models are used in the ``@Radiation`` tool.

.. option:: magnetic_field

   :Default value: ``__default__``
   :Allowed values: Name of magnetic field configuration object

   Specifies the magnetic field to use. Magnetic fields are configured separately
   in their own ``@MagneticField`` environments. Only one magnetic field can be
   set per simulation.

   If this parameter is left unspecified, then SOFT will try to select
   a magnetic field automatically. If exactly one magnetic field has been defined
   in the configuration file, then SOFT will select that magnetic field. Otherwise,
   if zero or more than one magnetic field has been defined, SOFT will throw an
   error and exit.

   More details can be found on the page about :ref:`module-magneticfield`.

.. option:: num_threads

   :Default value: ``OMP_NUM_THREADS``
   :Allowed values: Any positive integer

   Sets the number of OpenMP threads to run SOFT on. This value overrides the
   value specified in the environment variable ``OMP_NUM_THREADS``. If unspecified,
   this parameter is set ``OMP_NUM_THREADS``, which is typically equal to the
   number of threads available on your CPU, and will thus allow SOFT to utilize
   your CPU to the fullest.

.. option:: particle_generator

   :Default value: ``__default__``
   :Allowed values: Name of a particle generator configuration object

   Specifes the particle generator (phase-space) to use. Particle generators
   are configured separately in their own ``@ParticleGenerator`` environments.
   Only one particle generator can be set per simulation.

   If this parameter is left unspecified, then SOFT will try to select a
   particle generator automatically. If exactly one particle generator has been
   defined in the configuration file, then SOFT will select that particle
   generator. Otherwise, if zero or more than one particle generator has been
   defined, SOFT will throw an error and exit.

   More details can be found on the page about :ref:`module-particlegenerator`.

.. option:: particle_pusher

   :Default value: ``__default__``
   :Allowed values: Name of a particle pusher configuration object

   Specifies the particle pusher (orbit solver) to use. Particle pushers
   are configured separately in their own ``@ParticlePusher`` environments.
   Only one particle pusher can be set per simulation.

   If this parameter is left unspecified, then SOFT will try to select a
   particle pusher automatically. If exactly one particle pusher has been
   defined in the configuration file, then SOFT will select that particle
   pusher. Otherwise, if zero or more than one particle pusher has been
   defined, SOFT will throw an error and exit.

   More details can be found on the page about :ref:`module-particlepusher`.

