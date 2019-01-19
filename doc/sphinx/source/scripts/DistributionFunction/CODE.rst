.. _module-distribution-code:

(code)
------
The numerical tool CODE (for *COllisional Distribution of Electrons*) [#landreman2014]_
[#stahl2016]_ [#retoolsCODE]_, solves the spatially homogeneous Fokker--Planck
equation for electrons. The code includes many effects important to runaway
electron dynamics, including electric field acceleration, time-dependent plasma
parameters, avalanche multiplication, bremsstrahlung and synchrotron radiation
losses, as well as partial plasma ionization. SOFT provides support for loading
distribution functions calculated with CODE directly through the present module.

.. [#landreman2014] Landreman et al., 2014 "Numerical calculation of the runaway electron distribution function and associated synchrotron emission". *Computer Physics Communications* **185** (3), 847-855 `doi:10.1016/j.cpc.2013.12.004 <https://doi.org/10.1016/j.cpc.2013.12.004>`_.
.. [#stahl2016] Stahl et al., 2016 "Kinetic modelling of runaway electrons in dynamic scenarios", *Nuclear Fusion* **56** (11), 112009 `doi:10.1088/0029-5515/56/11/112009 <https://doi.org/10.1088/0029-5515/56/11/112009>`_.
.. [#retoolsCODE] http://ft.nephy.chalmers.se/retools/code/

Summary of options
^^^^^^^^^^^^^^^^^^

+----------------------------------+------------------------------------------------------------------------------+
| **Option**                       | **Description**                                                              |
+----------------------------------+------------------------------------------------------------------------------+
| :option:`code interptype`        | Interpolation method to use when interpolating in the distribution function. |
+----------------------------------+------------------------------------------------------------------------------+
| :option:`code name`              | Name of file containing distribution function.                               |
+----------------------------------+------------------------------------------------------------------------------+
| :option:`code radprof`           | Name of configuration block defining the radial profile to use.              |
+----------------------------------+------------------------------------------------------------------------------+
| :option:`code time`              | Index of distribution function time slice to use for simulation.             |
+----------------------------------+------------------------------------------------------------------------------+

Example configuration
^^^^^^^^^^^^^^^^^^^^^

The following example shows how to load a CODE distribution function and use
the final time step of the CODE simulation. The configuration also sets a radial
profile to use with the CODE momentum-space distribution function. If
:option:`code radprof` is not explicitly set, a uniform radial profile is used::

   distribution_function = ourDistribution;

   @DistributionFunction ourDistribution (code) {
       name    = "/path/to/code/distribution.mat";
       radprof = ourRadProf;
       time    = -1;         # -1 = last time step
   }

   @RadialProfile ourRadProf (linear) {}

File layout
^^^^^^^^^^^
The file containing the CODE distribution function must contain the variables
listed in the table below. They are all contained in the CODE output struct.

+--------------+---------------------------------------------------------------------------------------------------------------------------------------------------+
| **Variable** | **Description**                                                                                                                                   |
+--------------+---------------------------------------------------------------------------------------------------------------------------------------------------+
| ``delta``    | Normalization of momentum. Defined such that :math:`p = y\delta`, where :math:`p` is momentum normalized to the electron rest mass :math:`m_e c`. |
+--------------+---------------------------------------------------------------------------------------------------------------------------------------------------+
| ``f``        | Distribution function matrix. Must be of size :math:`n_t\times n_y n_\xi`.                                                                        |
+--------------+---------------------------------------------------------------------------------------------------------------------------------------------------+
| ``Nxi``      | Number of Legendre modes used in the CODE calculation to resolve the angular dependence.                                                          |
+--------------+---------------------------------------------------------------------------------------------------------------------------------------------------+
| ``y``        | Momentum grid, given in thermal units.                                                                                                            |
+--------------+---------------------------------------------------------------------------------------------------------------------------------------------------+

All options
^^^^^^^^^^^

.. program:: code

.. option:: interptype

   :Default value: ``cspline``
   :Allowed values: ``akima``, ``akima_periodic``, ``cspline``, ``cspline_periodic``, ``linear``, ``polynomial``, ``steffen``

   Determines which interpolation method to use for interpolating in the
   momentum dimension.

.. option:: name

   :Default value: None
   :Allowed values: String

   Specifies the name of the file containing the CODE distribution function
   to load.

.. option:: radprof

   :Default value: Uniform radial profile
   :Allowed values: Name of any defined :ref:`module-radialprofile`

   Specifies the radial profile object to use to generate a radial profile.

.. option:: time

   :Default value: ``-1`` (last timestep)
   :Allowed values: Any integer with absolute value less than the number of time points in the distribution function

   Selects the index of the time step to take the distribution function from.
   Negative indices count from the back of the array, so that ``-1`` corresponds
   to the last timestep, ``-2`` to the next-to-last etc.

