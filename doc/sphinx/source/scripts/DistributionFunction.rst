.. _module-distribution:

@DistributionFunction
*********************
The distribution function module, which defines a distribution function to
weigh results with. A number of different types of distribution functions can be
generated or loaded in SOFT2, including the analytical avalanche distribution [#fulop2006]_,
direct output from the kinetic solvers CODE, LUKE and GO+CODE, as well as regular
SOFT phase-space distribution functions.

.. [#fulop2006] Fülöp et al., 2006 "Destabilization of magnetosonic-whistler waves by a relativistic runaway beam". *Physics of Plasmas* **13** (6), 062506 `doi:10.1063/1.2208327 <https://doi.org/10.1063/1.2208327>`_.

.. toctree::
   :hidden:

   DistributionFunction/Avalanche
   DistributionFunction/CODE
   DistributionFunction/GOCODE
   DistributionFunction/LUKE
   DistributionFunction/Numerical
   DistributionFunction/Pitch
   DistributionFunction/Unit

Types
^^^^^
SOFT provides several different types of distribution functions. Which type is
configured by a particular block is specified by the "secondary type" of the
block, i.e. by giving the type in parentheses after the block name. The available
distribution function types are listed in the table below.

+--------------------------------------+-------------------------------------------------------------+
| **Type**                             | **Description**                                             |
+--------------------------------------+-------------------------------------------------------------+
| :ref:`module-distribution-avalanche` | An analytical distribution function based on [#fulop2006]_. |
+--------------------------------------+-------------------------------------------------------------+
| :ref:`module-distribution-code`      | A momentum-space distribution function directly from CODE.  |
+--------------------------------------+-------------------------------------------------------------+
| :ref:`module-distribution-gocode`    | GO-CODE distribution function.                              |
+--------------------------------------+-------------------------------------------------------------+
| :ref:`module-distribution-luke`      | Distribution function from the LUKE 1D2P kinetic solver.    |
+--------------------------------------+-------------------------------------------------------------+
| :ref:`module-distribution-numerical` | A SOFT numerical distribution function.                     |
+--------------------------------------+-------------------------------------------------------------+
| :ref:`module-distribution-pitch`     | An exponential pitch distribution function.                 |
+--------------------------------------+-------------------------------------------------------------+
| :ref:`module-distribution-unit`      | A distribution function that is one everywhere.             |
+--------------------------------------+-------------------------------------------------------------+

Example configuration
^^^^^^^^^^^^^^^^^^^^^
Please, see the pages for the various distribution function types to view
examples of how to configure a distribution function.

