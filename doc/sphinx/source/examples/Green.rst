.. _example-green:

Green's function output
-----------------------
SOFT is capable of generating so-called *Green's functions*.

.. image:: ../_static/figs/examples/Green.png
   :width: 100%
   :align: center

Important points
****************
For Green's function output, the :ref:`module-ro-green` module must be used. It
requires at least two parameters to be specified:
:option:`@RadiationOutput(green) format` (specifying the dependences of the
Green's function), and :option:`@RadiationOutput(green) output` (name of output
file).

The format option is a string consisting of any sequence of the following six
characters:

+------------+---------------------------------------------+
| **Format** | **Description**                             |
+------------+---------------------------------------------+
| ``1``      | (Alphabetically) first momentum parameter.  |
+------------+---------------------------------------------+
| ``2``      | (Alphabetically) second momentum parameter. |
+------------+---------------------------------------------+
| ``i``      | Vertical pixel dimension.                   |
+------------+---------------------------------------------+
| ``j``      | Horizontal pixel dimension.                 |
+------------+---------------------------------------------+
| ``r``      | Radial parameter.                           |
+------------+---------------------------------------------+
| ``w``      | Radiation spectrum.                         |
+------------+---------------------------------------------+

Dominant particles
******************

Example configuration
*********************
The following generates a Green's function :math:`G(p_\parallel, p_\perp)`::

   # Generate a momentum-space
   # Green's function.
   ##############################

   magnetic_field     = mf;
   tools              = rad;

   # Configuration of EAST-like magnetic equilibrium
   @MagneticField mf (analytical) {
       B0     = 5;     # On-axis field strength (T)
       Rm     = 0.68;  # Major radius (m)
       rminor = 0.22;  # Minor radius (m)

       # Safety-factor
       qtype  = linear;
       qa1    = 2;
       qa2    = 1;
   }

   # Phase space grid
   @ParticleGenerator PGen {
       a      = 0.0, 0.95, 20; # Normalized minor radius
       ppar   = -2, -165, 40;     # Parallel momentum (mc)
       pperp  = 1, 15, 40;    # Perpendicular momentum (mc)

       progress = 10;
   }

   # Orbit generator
   @ParticlePusher PPusher {
       nt = 2000;        # Number of timesteps per orbit (resolution parameter)
   }

   # Radiation tool
   @Radiation rad {
       detector = det;

       ntoroidal   = 7000;    # No. of toroidal sections in tokamak (resolution parameter)
       model       = cone;    # Radiation model to use
       output      = green;   # List of configuration of output
   }

   # Detector properties
   # Set up a tangentially viewing HXR camera.
   @Detector det {
       aperture     = 0.006;
       position     = 0.68, -0.68, 0;
       direction    = 0, 1, 0;
       vision_angle = 0.78 fov;
       spectrum     = no;
   }

   # Radiation model
   @RadiationModel cone (cone) {
       emission = synchrotron;
   }

   @RadiationOutput green (green) {
       format = 12;
       output = "data/green.mat";
   }

