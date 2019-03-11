.. _module-radialprofile:

@RadialProfile
**************
To enable combinations of momentum-space only kinetic solvers with varying
radial distributions of electrons, SOFT2 provides the ``RadialProfile`` module
which applies a momentum-independent radial profile to a distribution function.
This module must be combined with any of the :ref:`module-distribution` modules
that support coupling.

.. toctree::
   :hidden:

   RadialProfile/Linear
   RadialProfile/Power
   RadialProfile/Uniform

Types
^^^^^
SOFT provides several different types of radial profiles, each with its own set
of parameters. Which type of radial profile to use is specified using the
secondary type of the configuration block, i.e. by giving the type in
parentheses after the block name. The available distribution function types are
listed in the table below.

+-------------------------------------+--------------------------------------------------------------------------------+
| **Type**                            | **Function**                                                                   |
+-------------------------------------+--------------------------------------------------------------------------------+
| :ref:`module-radialprofile-linear`  | :math:`1 - \frac{r - r_{\rm min}}{r_{\rm max} - r_{\rm min}}`                  |
+-------------------------------------+--------------------------------------------------------------------------------+
| :ref:`module-radialprofile-power`   | :math:`1 - \left( \frac{r - r_{\rm min}}{r_{\rm max} - r_{\rm min}} \right)^b` |
+-------------------------------------+--------------------------------------------------------------------------------+
| :ref:`module-radialprofile-uniform` | :math:`1`                                                                      |
+-------------------------------------+--------------------------------------------------------------------------------+

Here, :math:`r` denotes the minor radial location of the particle,
:math:`r_{\rm min}` and :math:`r_{\rm max}` are the minor radial coordinates
of the inner and outer edges of the electron beam, and :math:`b` is a model
parameter. Note that the default values :math:`r_{\rm min}` and
:math:`r_{\rm max}` are 0 and the separatrix radius respectively. The radii
may not be negative, and if :math:`r_{\rm min} > 0`, the electron beam will be
hollow.

Example configuration
^^^^^^^^^^^^^^^^^^^^^
Please, see the pages for the various radial profile types to view examples of
how to configure each.

