.. _example-orbits:

Orbits
------
Example of how to generate guiding-center and particle orbits using SOFT.

.. image:: ../_static/figs/examples/orbits-comb.png
   :align: center

*Example of orbits computed using SOFT. The left figure shows guiding-center
orbits, while the right figure shows particle orbits.*

Important points
****************
SOFT always calculates orbits, but to explicitly output the orbits calculated,
the :ref:`module-orbits` module must be used.

One-file example
****************
The following is a one-file example of a configuration file that will generate
the guiding-center orbits in the figure at the top of this page::

   # Simulate orbits in a circular tokamak
   # magnetic field.
   ##############################

   magnetic_field     = mf;
   tools              = orbits;
   include_drifts     = yes;

   # Configuration of EAST-like magnetic equilibrium
   @MagneticField mf (analytical) {
       B0     = 5;     # On-axis field strength (T)
       Rm     = 0.68;  # Major radius (m)
       rminor = 0.22;  # Minor radius (m)

       # Safety-factor
       qtype  = linear;
       qa1    = 2;
       qa2    = 1;
       sigmaB = ccw;
   }

   # Phase space grid
   @ParticleGenerator PGen {
       a      = 0.1, 0.9, 9;  # Minor radius
       gamma  = 10, 10, 1;    # Energy (mc^2) (approx. 30 MeV)
       xi     = 0.3, 0.3, 1;  # Cosine of pitch angle

       position = guiding-center;
   }

   # Orbit generator
   @ParticlePusher PPusher {
       nt       = 100;        # Number of timesteps per orbit (resolution parameter)
   }

   @Orbits orbits {
       output = "data/guiding-center.mat";
   }

Configuration scripts
*********************
The following configuration scripts are available in the `SOFT2 GitHub
repository <https://github.com/hoppe93/SOFT2/>`_ and will generate the two
figures at the top of this page.

.. |baseline| replace:: **baseline**
.. _baseline: https://github.com/hoppe93/SOFT2/blob/master/examples/Orbits/baseline
.. |gc| replace:: gc
.. _gc: https://github.com/hoppe93/SOFT2/blob/master/examples/Orbits/gc
.. |particle| replace:: particle
.. _particle: https://github.com/hoppe93/SOFT2/blob/master/examples/Orbits/particle

+-------------+----------------------------------------------------+
| **File**    | **Description**                                    |
+-------------+----------------------------------------------------+
| |baseline|_ | Baseline configuration applied to all simulations. |
+-------------+----------------------------------------------------+
| |gc|_       | Simulate guiding-center orbits.                    |
+-------------+----------------------------------------------------+
| |particle|_ | Simulate particle orbits.                          |
+-------------+----------------------------------------------------+

