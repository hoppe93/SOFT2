.. _example-green:

Green's function output
-----------------------
SOFT is capable of generating so-called *Green's functions*.

Important points
****************
For Green's function output, the :ref:`module-ro-green` module must be used. It
requires at least two parameters to be specified:
:option:`@RadiationOutput(green) format` (specifying the dependences of the
Green's function), and :option:`@RadiationOutput(green) output` (name of output
file).

Dominant particles
******************

Example configuration
*********************
The following generates a Green's function :math:`G(p_\parallel, p_\perp)`::


