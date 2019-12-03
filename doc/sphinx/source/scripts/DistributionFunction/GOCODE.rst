.. _module-distribution-gocode:

(gocode)
--------
GO+CODE refers to the coupled fluid-kinetic framework consisting of the fluid
code GO [#smith2006]_ [#papp2013]_, which is used to solve for the toroidal
electric field evolution, and the kinetic solver CODE [#landreman2014]_
[#stahl2016]_ [#retoolsCODE]_, which solves the spatially homogeneous
Fokker--Planck equation for electrons. The framework calculates the electron
distribution function with one spatial dimension and two momentum dimensions.

.. [#smith2006] Smith et al., 2006 "Runaway electrons and the evolution of the plasma current in tokamak disruptions". *Physics of Plasmas* **13**, 102502 `doi:10.1063/1.2358110 <https://doi.org/10.1063/1.2358110>`_.
.. [#papp2013] Papp et al., 2013 "The effect of ITER-like wall on runaway electron generation in JET". *Nuclear Fusion* **53** (12), 123017 `doi:10.1088/0029-5515/53/12/123017 <https://doi.org/10.1088/0029-5515/53/12/123017>`_.
.. [#landreman2014] Landreman et al., 2014 "Numerical calculation of the runaway electron distribution function and associated synchrotron emission". *Computer Physics Communications* **185** (3), 847-855 `doi:10.1016/j.cpc.2013.12.004 <https://doi.org/10.1016/j.cpc.2013.12.004>`_.
.. [#stahl2016] Stahl et al., 2016 "Kinetic modelling of runaway electrons in dynamic scenarios", *Nuclear Fusion* **56** (11), 112009 `doi:10.1088/0029-5515/56/11/112009 <https://doi.org/10.1088/0029-5515/56/11/112009>`_.
.. [#retoolsCODE] http://ft.nephy.chalmers.se/retools/code/

Summary of options
^^^^^^^^^^^^^^^^^^

+----------------------------------+------------------------------------------------------------------------------+
| **Option**                       | **Description**                                                              |
+----------------------------------+------------------------------------------------------------------------------+
| :option:`gocode interptype`      | Interpolation method to use when interpolating in the distribution function. |
+----------------------------------+------------------------------------------------------------------------------+
| :option:`gocode name`            | Name of file containing distribution function.                               |
+----------------------------------+------------------------------------------------------------------------------+
| :option:`gocode time`            | Index of distribution function time slice to use for simulation.             |
+----------------------------------+------------------------------------------------------------------------------+

Example configuration
^^^^^^^^^^^^^^^^^^^^^

The following example shows how to load a GO+CODE distribution function and use
the final time step of the simulation::

   distribution_function = ourDistribution;

   @DistributionFunction ourDistribution (gocode) {
       name    = "/path/to/gocode/distribution.mat";
       time    = -1;         # -1 = last time step
   }

File layout
^^^^^^^^^^^
In contrast to pure CODE and LUKE distribution functions, a GO+CODE cannot be
directly saved to a MAT file. Instead, the output must be converted into a
special format, consisting of two levels of structs.

In the root path of the file, the following variables must be set:

+--------------+---------------------------------------------------------------------------------------------------------------------------------------------------+
| **Variable** | **Description**                                                                                                                                   |
+--------------+---------------------------------------------------------------------------------------------------------------------------------------------------+
| ``r``        | Radial grid on which the distribution is defined (corresponding to the dimensionful ``r`` variable in GO).                                        |
+--------------+---------------------------------------------------------------------------------------------------------------------------------------------------+
| ``t``        | List of times at which the distribution is available.                                                                                             |
+--------------+---------------------------------------------------------------------------------------------------------------------------------------------------+
| ``nt``       | Number of elements in ``t``.                                                                                                                      |
+--------------+---------------------------------------------------------------------------------------------------------------------------------------------------+
| ``tX``       | Struct, containing the distribution functions at time index ``X`` (i.e. ``X`` is an integer, ranging from ``0`` to ``nt-1``.                      |
+--------------+---------------------------------------------------------------------------------------------------------------------------------------------------+

Each struct ``tX`` (``t0``, ``t1`` etc.) in turn contains several variables
``rY``, with ``Y`` an integer ranging from ``0`` to ``nr-1``, where ``nr`` is
the number of elements in ``r``. Finally, each ``rY`` contains a momentum
distribution which has the format described on the page for the CODE
distribution function: :ref:`module-distribution-code`. We repeat it here for
convenience:

+--------------+---------------------------------------------------------------------------------------------------------------------------------------------------+
| **Variable** | **Description**                                                                                                                                   |
+--------------+---------------------------------------------------------------------------------------------------------------------------------------------------+
| ``delta``    | Normalization of momentum. Defined such that :math:`p = y\delta`, where :math:`p` is momentum normalized to the electron rest mass :math:`m_e c`. |
+--------------+---------------------------------------------------------------------------------------------------------------------------------------------------+
| ``f``        | Distribution function matrix. Must be of size :math:`n_t\times n_y n_\xi`.                                                                        |
+--------------+---------------------------------------------------------------------------------------------------------------------------------------------------+
| ``Nxi``      | Number of Legendre modes used in the CODE calculation to resolve the angular dependence.                                                          |
+--------------+---------------------------------------------------------------------------------------------------------------------------------------------------+
| ``nref``     | Reference density (used for un-normalising the distribution function).                                                                            |
+--------------+---------------------------------------------------------------------------------------------------------------------------------------------------+
| ``Tref``     | Reference temperature (used for un-normalising the distribution function).                                                                        |
+--------------+---------------------------------------------------------------------------------------------------------------------------------------------------+
| ``y``        | Momentum grid, given in thermal units.                                                                                                            |
+--------------+---------------------------------------------------------------------------------------------------------------------------------------------------+

All options
^^^^^^^^^^^

.. program:: gocode

.. option:: interptype

   :Default value: ``cspline``
   :Allowed values: ``akima``, ``akima_periodic``, ``cspline``, ``cspline_periodic``, ``linear``, ``polynomial``, ``steffen``

   Determines which interpolation method to use for interpolating in the
   momentum dimension.

.. option:: name

   :Default value: None
   :Allowed values: String

   Specifies the name of the file containing the GO+CODE distribution function
   to load.

.. option:: time

   :Default value: ``-1`` (last timestep)
   :Allowed values: Any integer with absolute value less than the number of time points in the distribution function

   Selects the index of the time step to take the distribution function from.
   Negative indices count from the back of the array, so that ``-1`` corresponds
   to the last timestep, ``-2`` to the next-to-last etc.

