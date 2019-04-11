.. _module-particlegenerator:

@ParticleGenerator
******************
The particle generator is in charge of sampling particles from the phase space.
Thus, the particle generator must know which type of particle to sample (i.e.
electron, proton or some other charged particle?) as well as the bounds of the
phase space from which to pick the particle properties (energy, initial pitch
angle and initial position).

Summary of options
^^^^^^^^^^^^^^^^^^

+--------------------------------------------------+------------------------------------------------------------------+
| **Option**                                       | **Description**                                                  |
+--------------------------------------------------+------------------------------------------------------------------+
| :option:`@ParticleGenerator charge`              | Charge of particle (in units of elementary charge).              |
+--------------------------------------------------+------------------------------------------------------------------+
| :option:`@ParticleGenerator driftshifttol`       | Relative tolerance when computing the drift orbit shift.         |
+--------------------------------------------------+------------------------------------------------------------------+
| :option:`@ParticleGenerator mass`                | Rest mass of particle (in units of the electron rest mass).      |
+--------------------------------------------------+------------------------------------------------------------------+
| :option:`@ParticleGenerator mpi_distribute_mode` | Specifies which parameter to split among MPI processes.          |
+--------------------------------------------------+------------------------------------------------------------------+
| :option:`@ParticleGenerator position`            | Whether specified position refers to particle or guiding-center. |
+--------------------------------------------------+------------------------------------------------------------------+
| :option:`@ParticleGenerator progress`            | Instructs SOFT to output info on the simulation progress.        |
+--------------------------------------------------+------------------------------------------------------------------+

In addition to these options, the following phase space parameters can also be
specified. Exactly one spatial and two momentum parameters must be specified,
and the two momentum parameters chosen must span all of momentum space together.

+--------------------------------------+----------+-----------------------------------------------------------------------------------------------------+
| **Parameter**                        | **Type** | **Description**                                                                                     |
+--------------------------------------+----------+-----------------------------------------------------------------------------------------------------+
| :option:`@ParticleGenerator a`       | Spatial  | Normalized minor radius (0 on magnetic axis; 1 on last closed flux surface).                        |
+--------------------------------------+----------+-----------------------------------------------------------------------------------------------------+
| :option:`@ParticleGenerator r`       | Spatial  | Minor radius (meters).                                                                              |
+--------------------------------------+----------+-----------------------------------------------------------------------------------------------------+
| :option:`@ParticleGenerator rho`     | Spatial  | Major radius (meters).                                                                              |
+--------------------------------------+----------+-----------------------------------------------------------------------------------------------------+
| :option:`@ParticleGenerator gamma`   | Momentum | Particle energy, normalized to the particle rest mass (= Lorentz factor).                           |
+--------------------------------------+----------+-----------------------------------------------------------------------------------------------------+
| :option:`@ParticleGenerator p`       | Momentum | Particle momentum, normalized to the particle rest mass.                                            |
+--------------------------------------+----------+-----------------------------------------------------------------------------------------------------+
| :option:`@ParticleGenerator ppar`    | Momentum | Particle momentum parallel to magnetic field, normalized to the particle rest mass.                 |
+--------------------------------------+----------+-----------------------------------------------------------------------------------------------------+
| :option:`@ParticleGenerator pperp`   | Momentum | Particle momentum perpendicular to magnetic field, normalized to the particle rest mass.            |
+--------------------------------------+----------+-----------------------------------------------------------------------------------------------------+
| :option:`@ParticleGenerator thetap`  | Momentum | Particle pitch angle :math:`\theta_{\rm p}` (in radians).                                           |
+--------------------------------------+----------+-----------------------------------------------------------------------------------------------------+
| :option:`@ParticleGenerator ithetap` | Momentum | Pi-complement of particle pitch angle, :math:`\theta_{\rm p}' = \pi - \theta_{\rm p}` (in radians). |
+--------------------------------------+----------+-----------------------------------------------------------------------------------------------------+
| :option:`@ParticleGenerator xi`      | Momentum | Cosine of particle pitch angle, :math:`\xi = \cos\theta_{\rm p}`.                                   |
+--------------------------------------+----------+-----------------------------------------------------------------------------------------------------+

Example configuration
^^^^^^^^^^^^^^^^^^^^^

**Electron** --- The following example sets up a phase-space for an electron
with 100 points on the grid in each dimension. The mass and charge default to
those of an electron, and so do not have to be specified. We also instruct SOFT
to output a total of 100 progress messages during the run. Since we do not set
the meaning of the position explicitly, SOFT assumes that we specify the
position of the guiding-center::

   @ParticleGenerator PGen_electron {
       r        = 0, 0.30, 100;    # Minor radius (m)
       gamma    = 1.1, 3.0, 100;   # Energy       (mc^2)
       thetap   = 0.02, 0.8, 100;  # Pitch angle  (rad)
       progress = 100;             # Output 100 progress messages
   }

**Deuterium** --- The following example sets up a phase-space for a deuterium
ion with 100 points on the grid in each dimension. The proton-electron mass
ratio is approximately :math:`m_{\rm p} / m_{\rm e}\approx 1836`, and hence the
deuterium-electron mass ratio is approximately
:math:`m_{\rm D} / m_{\rm e}\approx 3672`. We explicitly state that we specify
the *particle* position using the ``position`` option::

   @ParticleGenerator PGen_deuterium {
       mass   = 3672;            # Electron masses
       charge = 2;               # Elementary charges
       a      = 0, 1, 100;       # Normalized minor radius
       p      = 1e-3, 1e-1, 100; # Momentum (mc)
       thetap = 0.1, 0.3, 100;   # Pitch angle (rad)
   }

**With MPI** --- The following example is near-identical to the electron example
above, but explicitly instructs SOFT to split the ``gamma`` (energy) parameter
among the MPI processes::

   @ParticleGenerator PGen_electron {
       r                   = 0, 0.30, 100;    # Minor radius (m)
       gamma               = 1.1, 3.0, 100;   # Energy       (mc^2)
       thetap              = 0.02, 0.8, 100;  # Pitch angle  (rad)
       progress            = 100;             # Output 100 progress messages
       mpi_distribute_mode = gamma;           # Divide the energy parameter among MPI processes
   }

Options
^^^^^^^

.. program:: @ParticleGenerator

.. option:: charge

   :Default value: ``-1``
   :Allowed values: Any non-zero real number.

   Charge of particle to simulate. The charge is given in units of the
   elementary charge so that a value of ``1`` corresponds to the *proton*
   charge and ``-1`` to the *electron* charge.

.. option:: driftshifttol

   :Default value: ``1e-4``
   :Allowed values: :math:`\epsilon < \text{tolerance} < 1`

   Tolerance for determining the guiding-center drift orbit shift (which is used
   to determine where to launch particles and where to sample the distribution
   function). In general, this parameter does not need to be changed.

.. option:: mass

   :Default value: ``1``
   :Allowed values: Any positive real number.

   Rest mass of particle to simulate. The mass is given in units of the electron
   rest mass. In these units, the proton mass is
   :math:`m_{\rm p}\approx 1836.15267389` [#wikimassratio]_.

.. option:: mpi_distribute_mode

   :Default value: ``auto``.
   :Allowed values: ``1``, ``2``, ``auto``, ``radius`` and all of the phase-space parameters listed under :ref:`partgen-phase-space-params`.

   When running in MPI mode (MPI = Message Passing Interface; distributed
   memory parallelization), this parameter can be used to specify how the phase
   space should be divided among the MPI processes. In contrast to regular
   OpenMP parallelization, which builds a queue of points in phase space, MPI
   requires one phase space parameter to be divided evenly among the processes.

   Which parameter to divide is specified by giving the name of the parameter,
   as listed under :ref:`partgen-phase-space-params`, or by giving one of
   ``1``, ``2`` and ``radius``. The former two cause SOFT to divide the first
   and second momentum parameter respectively (i.e. the alphabetically first
   and second momentum parameter), while the latter divides the radial
   parameter, whichever it may be.

   If ``auto`` is specified, SOFT2 chooses the phase space parameter with the
   largest number of grid points. This is the default setting.

.. option:: position

   :Default value: ``guiding-center``
   :Allowed values: ``gc``, ``guiding-center`` and ``particle``.

   Specifies whether the *guiding-center* or *particle* is given as input. If,
   for example, the particle position is specified, but guiding-center orbits
   are simulated, then the guiding-centers are offset from the given position
   by one Larmor radius, and vice versa for the opposite case.

.. option:: progress

   :Default value: ``no``
   :Allowed values: ``yes``, ``no`` or a positive integer.

   If ``yes`` or a positive integer ``n``, outputs a message reporting the
   progress of the simulation a total of ``n`` times during the run. The
   reports are split evenly accross the phase space, meaning that if the
   phase space consists of ``N`` total grid points, then SOFT reports progress
   roughly when the number of processed grid points is a multiple of ``N / n``.

.. [#wikimassratio] https://en.wikipedia.org/wiki/Proton-to-electron_mass_ratio

.. _partgen-phase-space-params:

Phase space parameters
^^^^^^^^^^^^^^^^^^^^^^

.. option:: a

   **Normalized minor radius** --- The initial minor radial location of the
   particle/guiding-center, normalized to the minor radial value of the last
   closed flux surface. Thus, :math:`a = 0` corresponds to the magnetic axis
   and :math:`a = 1` to the maximum radius of the last closed flux surface.

.. option:: r

   **Minor radius** --- The initial minor radial location of the
   particle/guiding-center, given in meters.

.. option:: rho

   **Major radius** --- The initial major radial location of the
   particle/guiding-center, given in meters.

.. option:: gamma

   **Energy** --- The energy of the particle/guiding-center, normalized to the
   particle rest mass :math:`mc^2`, where :math:`c` denotes the speed of light
   in vacuum. This quantity is also known as the *Lorentz factor* or
   *relativistic factor*, and can also be written
   :math:`\gamma = \left( 1 - v^2/c^2 \right)^{-1/2}`, where :math:`v` is the
   speed of the particle.

.. option:: p

   **Momentum** --- The momentum of the particle/guiding-center, normalized to
   the particle rest mass :math:`mc`, where :math:`c` denotes the speed of light
   in vacuum. This quantity is related to the particle energy/relativistic
   factor through :math:`\gamma^2 = 1 + p^2`.

.. option:: ppar

   **Parallel momentum** --- Momentum of the particle parallel to the magnetic
   field, normalized to the particle rest mass :math:`mc`, where :math:`c`
   denotes the speed of light in vacuum.

.. option:: pperp

   **Perpendicular momentum** --- Momentum of the particle perpendicular to the
   magnetic field, normalized to the particle rest mass :math:`mc`, where
   :math:`c` denotes the speed of light in vacuum.

.. option:: thetap

   **Pitch angle** --- Angle between the particle velocity vector and the
   magnetic field vector. Given in radians. The pitch angle ranges between
   :math:`0` and :math:`\pi`. A value greater than :math:`\pi/2` means that the
   particle is moving antiparallel to the magnetic field.

.. option:: ithetap

   **Complementary pitch angle** --- Same as :option:`@ParticleGenerator thetap`,
   except that it is defined as "pi minus :option:`@ParticleGenerator thetap`",
   i.e. :math:`\theta_{\rm p}' = \pi - \theta_{\rm p}`. This is useful when
   simulating particles with negative parallel momentum (moving in the
   antiaparallel direction of the magnetic field), since instead of specifying
   :math:`\theta_{\rm p} = 3.14159265359` we can then set
   :math:`\theta_{\rm p}' = 0`.

.. option:: xi

   **Cosine of pitch angle** --- Cosine of the pitch angle
   :math:`\theta_{\rm p}`, i.e. :math:`\xi = \cos\theta_{\rm p}`.

Jacobians
^^^^^^^^^

The following is a list of all the Jacobian determinants for transformations
from a Cartesian coordinate system :math:`(p_x, p_y, p_z)` to other coordinate
systems :math:`(p_1, p_2, \zeta)`, where :math:`\zeta` is the gyro angle.

**gamma / ppar** --- :math:`(\gamma, p_{\parallel})`

.. math::

   \mathrm{d}p_x\mathrm{d}p_y\mathrm{d}p_z = \gamma\mathrm{d}\gamma\mathrm{d}p_{\parallel}\mathrm{d}\zeta

**gamma / pperp** --- *Does not contain sufficient information to determine if
the guiding-center is travelling in the parallel or anti-parallel direction of
the magnetic field.*

**gamma / thetap** --- :math:`(\gamma, \theta_{\rm p})`

.. math::

   \mathrm{d}p_x\mathrm{d}p_y\mathrm{d}p_z = \gamma\sin\theta_{\rm p}\sqrt{\gamma^2-1}\,\mathrm{d}\gamma\mathrm{d}\theta_{\rm p}\mathrm{d}\zeta

**gamma / xi** --- :math:`(\gamma, \xi)`

.. math::

   \mathrm{d}p_x\mathrm{d}p_y\mathrm{d}p_z = \gamma\sqrt{\gamma^2-1}\,\mathrm{d}\gamma\mathrm{d}\xi\mathrm{d}\zeta

**p / ppar** --- :math:`(p, p_{\parallel})`

.. math::

   \mathrm{d}p_x\mathrm{d}p_y\mathrm{d}p_z = p\,\mathrm{d}p\mathrm{d}p_{\parallel}\mathrm{d}\zeta

**p / pperp** --- *Does not contain sufficient information to determine if
the guiding-center is travelling in the parallel or anti-parallel direction of
the magnetic field.*

**p / thetap** --- :math:`(p, \theta_{\rm p})`

.. math::

   \mathrm{d}p_x\mathrm{d}p_y\mathrm{d}p_z = p^2\sin\theta_{\rm p}\,\mathrm{d}p\mathrm{d}\theta_{\rm p}\mathrm{d}\zeta

**p / \xi** --- :math:`(p, \xi)`

.. math::

   \mathrm{d}p_x\mathrm{d}p_y\mathrm{d}p_z = p^2\,\mathrm{d}p\mathrm{d}\theta_{\rm p}\mathrm{d}\zeta

**ppar / pperp** --- :math:`(p_{\parallel}, p_{\perp})`

.. math::

   \mathrm{d}p_x\mathrm{d}p_y\mathrm{d}p_z = p_\perp\,\mathrm{d}p_{\parallel}\mathrm{d}p_{\perp}\mathrm{d}\zeta

**ppar / thetap** --- :math:`(p_{\parallel}, \theta_{\rm p})`

.. math::

   \mathrm{d}p_x\mathrm{d}p_y\mathrm{d}p_z = \frac{p_\parallel^2\sin\theta_{\rm p}}{\cos^3\theta_{\rm p}}\,\mathrm{d}p_{\parallel}\mathrm{d}\theta_{\rm p}\mathrm{d}\zeta

**ppar / xi** --- :math:`(p_{\parallel}, \xi)`

.. math::

   \mathrm{d}p_x\mathrm{d}p_y\mathrm{d}p_z = \frac{p_\parallel^2}{\xi^3}\,\mathrm{d}p_{\parallel}\mathrm{d}\xi\mathrm{d}\zeta

**pperp / thetap** --- :math:`(p_{\parallel}, \theta_{\rm p})`

.. math::

   \mathrm{d}p_x\mathrm{d}p_y\mathrm{d}p_z = \frac{p_\perp^2}{\sin^2\theta_{\rm p}}\,\mathrm{d}p_{\perp}\mathrm{d}\theta_{\rm p}\mathrm{d}\zeta

**pperp / xi** --- :math:`(p_\perp, \xi)`

.. math::

   \mathrm{d}p_x\mathrm{d}p_y\mathrm{d}p_z = \frac{p_\perp^2}{\left( 1 - \xi^2 \right)^{3/2}}\,\mathrm{d}p_{\perp}\mathrm{d}\xi\mathrm{d}\zeta

