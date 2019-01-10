.. _module-radiation:

@Radiation
**********
This tool calculates various radiation quantities, including radiation images,
spectra and Green's functions. It is the tool to use for studying bremsstrahlung
and synchrotron radiation.

The basic purpose of the ``@Radiation`` tool is to evaluate various forms of the
radiation diagnostic integral [#hoppe2019lic]_.

.. math::
   :label: raddiagint

   I = \int \Theta\left( \frac{\boldsymbol{r}}{r} \right) \frac{\boldsymbol{r}\cdot\hat{\boldsymbol{n}}}{r^3}
   \frac{\mathrm{d}I(\boldsymbol{x},\boldsymbol{p},\boldsymbol{r})}{\mathrm{d}\Omega}
   f(\boldsymbol{x},\boldsymbol{p})\,\mathrm{d}\boldsymbol{p}\,\mathrm{d} V\,\mathrm{d} A.


where :math:`I` denotes a general radiation quantity (e.g. radiation power,
spectral power, etc.), :math:`\Theta` is the field-of-view step function,
:math:`\boldsymbol{r}` is a vector from the particle to the detector,
:math:`\hat{\boldsymbol{n}}` is the detector surface normal vector,
:math:`\mathrm{d}I/\mathrm{d}\Omega` is the angular distribution of radiation
and :math:`f(\boldsymbol{x},\boldsymbol{p})` is the distribution function.
The integral is taken over all of momentum space (indicated by the differential
:math:`\mathrm{d}\boldsymbol{p}`), real space (indicated by the volume element
:math:`\mathrm{d}V`), and the detector surface (indicated by
:math:`\mathrm{d} A`).

To simplify the computation, SOFT evaluates the above integral in guiding-center
coordinates. Using these coordinates, the radiation diagnostic integral
:eq:`raddiagint` can be written

.. math::

   I = \int\Theta\left( \frac{\boldsymbol{r}}{r} \right) \frac{\boldsymbol{r}\cdot\hat{\boldsymbol{n}}}{r^3}
   \frac{\mathrm{d}I(\boldsymbol{X},\boldsymbol{p},\boldsymbol{r})}{\mathrm{d}\Omega}
   f(\rho,p_\parallel,p_\perp)\,J\,\underbrace{\mathrm{d}p_\parallel\mathrm{d}p_\perp\mathrm{d}\zeta}_{\mathrm{d}\boldsymbol{p}}
   \,\underbrace{\mathrm{d}\rho\mathrm{d}\tau\mathrm{d}\phi}_{\mathrm{d}\boldsymbol{X}}\,\mathrm{d} A.

where now :math:`\boldsymbol{X}` denotes the guiding-center position,
:math:`p_\parallel` and :math:`p_\perp` are the particle momenta in the
directions parallel and perpendicular respectively to the magnetic field,
:math:`\zeta` is the gyro angle, :math:`\rho` is the maximum major radius
visited by a guiding-center along its orbit, :math:`\tau` is the time along the
orbit of the guiding-center and :math:`\phi` is the toroidal angle.

.. [#hoppe2018a] Hoppe et al., 2018, "SOFT: a synthetic synchrotron diagnostic for runaway electrons". *Nuclear Fusion* **58** (2), 026032 `doi:10.1088/1741-4326/aa9abb <https://doi.org/10.1088/1741-4326/aa9abb>`_.

.. [#hoppe2019lic] Hoppe, 2019, "Simulation and analysis of radiation from runaway electrons". *Licentiate thesis* `Available online <http://ft.nephy.chalmers.se/publications/Hoppe_Licentiate_Thesis.pdf>`_.

Summary of options
------------------

+-----------------------------------+----------------------------------------------------------+
| **Option**                        | **Description**                                          |
+-----------------------------------+----------------------------------------------------------+
| :option:`@Radiation detector`     | Specifies which detector configuration to use            |
+-----------------------------------+----------------------------------------------------------+
| :option:`@Radiation model`        | Specifies which radiation model to use                   |
+-----------------------------------+----------------------------------------------------------+
| :option:`@Radiation ntoroidal`    | Toroidal resolution parameter                            |
+-----------------------------------+----------------------------------------------------------+
| :option:`@Radiation output`       | Specifies which output module(s) to use                  |
+-----------------------------------+----------------------------------------------------------+
| :option:`@Radiation torthreshold` | Parameter for :ref:`toroidal-optim`                      |
+-----------------------------------+----------------------------------------------------------+
| :option:`@Radiation torquad`      | Quadrature rule to use when evaluating toroidal integral |
+-----------------------------------+----------------------------------------------------------+
| :option:`@Radiation wall_opacity` | Specifies the wall "opacity"                             |
+-----------------------------------+----------------------------------------------------------+

Example configuration
---------------------
The ``@Radiation`` is merely the parent of a set modules which together produce
the desired simulation output. As such, we must specify both the detector, the
radiation model and output type. An example configuration of a ``@Radiation``
module, along with its required sub-modules, is::

   @Radiation rad {
       detector     = det;
       model        = cone;
       ntoroidal    = 7500;
       output       = image topview;
   }

   @Detector det {
       aperture     = 0.006;
       direction    = 0, 1, 0;
       position     = 0, 1.7, 0;
       vision_angle = 1.25 fov;
       spectrum     = 440e-9, 790e-9, 40;
   }

   @RadiationModel cone (cone) {
       emission = synchrotron;
   }

   @RadiationOutput image (image) {
       pixels = 600;
       output = "myimage.mat";
   }


Available sub-modules
---------------------
There are three types of sub-modules that must be configured for the
``@Radiation`` module. In addition to a :ref:`module-detector`, one radiation
model must specified as well as *at least* one output module.

.. _module-radiation-output:

Output sub-modules
^^^^^^^^^^^^^^^^^^

Radiation output modules are specified with the block type
:ref:`module-radiationoutput`. The secondary type of the block (in parentheses
after the block name) determines which type of output the block configures.
The available secondary types of :ref:`module-radiationoutput` are

+---------------------------------+-------------------------------+
| **Module name**                 | **Output description**        |
+---------------------------------+-------------------------------+
| :ref:`module-ro-greensfunction` | Green's/weight functions      |
+---------------------------------+-------------------------------+
| :ref:`module-ro-image`          | Camera images                 |
+---------------------------------+-------------------------------+
| :ref:`module-ro-space3d`        | 3D maps of radiation          |
+---------------------------------+-------------------------------+
| :ref:`module-ro-spectrum`       | Radiation spectra             |
+---------------------------------+-------------------------------+
| :ref:`module-ro-topview`        | Tokamak topviews of radiation |
+---------------------------------+-------------------------------+


.. _module-radiation-models:

Radiation model sub-modules
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Radiation model modules are specified with the block type
:ref:`module-radiationmodel`. The secondary type of the block (in parentheses
after the block name) determines which type of model the block configures.
The available secondary types of :ref:`module-radiationmodel` are

+--------------------------------------+-------------------------------------------------------+
| **Module name**                      | **Model description**                                 |
+--------------------------------------+-------------------------------------------------------+
| :ref:`module-rm-angulardistribution` | Full angular (and spectral) distribution of radiation |
+--------------------------------------+-------------------------------------------------------+
| :ref:`module-rm-cone`                | Special model for approximating directed radiation    |
+--------------------------------------+-------------------------------------------------------+
| :ref:`module-rm-isotropic`           | Special model for perfectly isotropic radiation       |
+--------------------------------------+-------------------------------------------------------+

.. _toroidal-optim:

Toroidal optimization
---------------------

Options
-------

.. program:: @Radiation

.. option:: detector

   :Default value: Nothing
   :Allowed values: Name of any :ref:`module-detector` configuration block

   Specifies the name of the configuration block to use for setting the
   properties of the detector.

.. option:: model

   :Default value: Nothing
   :Allowed values: Name of any radiation model configuration block

   Specifies the name of the configuration block to use for setting the
   radiation model to use. The radiation model basically specifies how the
   angular distribution of radiation is handled. SOFT can take the full
   angular distribution of radiation into account, but usually, for synchrotron
   radiation, the approximative model known as the "cone model" is often used
   instead. A list of available radiation models can be found above, under the
   section :ref:`module-radiation-models`.

.. option:: ntoroidal

   :Default value: ``3500``
   :Allowed values: Any positive integer

   Number of toroidal sections to divide the tokamak into. This is the
   resolution parameter for the toroidal integral in the radiation diagnostic
   integral evaluated by the ``@Radiation`` tool.

.. option:: output

   :Default value: Nothing
   :Allowed values: List of names of radiation output module configuration blocks

   List of names of configuration blocks setting the properties of the output
   modules to use.

   The ``@Radiation`` tool only facilitates the computation of various radiation
   quantities (such as images and spectra). The actual evaluation of these
   quantities, as well as subsequent generation of output files, are handled by
   the corresponding "radiation output" modules. A full list of available
   radiation output modules can be found above under the section
   :ref:`module-radiation-output`.

.. option:: torthreshold

   :Default value: ``0``
   :Allowed values: Any real value between or equal to ``0`` and ``1``

   Threshold for neglecting the integrand when using the ``maximize`` quadrature
   to evaluate the toroidal integral. The integration stops as soon as the value
   of the integrand is a fraction ``torthreshold`` of the maximum integrand
   value seen so far.
   
   For the cone model, this parameter can safely be set to ``0``. When used
   together with the models that take the full angular distribution into
   account, this parameter should be set to a value greater than ``0`` (yet
   less than ``1``).

.. option:: torquad

   :Default value: ``maximize``
   :Allowed values: ``maximize``, ``trapz``

   Determines which quadrature rule to use for evaluating the toroidal integral.
   The ``trapz`` quadrature is a simple trapezoidal rule. The ``maximize`` rule
   is based on the trapezoidal rule, but uses an optimization algorithm to
   determine which parts of space that will contribute with radiation. The
   ``maximize`` quadrature is often between a factor 25-100 faster than the
   regular trapezoidal rule.

.. option:: wall_opacity

   :Default value: ``semi``
   :Allowed values: ``opaque``, ``semi``, ``transparent``

   Sets the "opacity" level of the wall. If ``opaque``, all walls are fully
   accounted for, and radiation is not allowed to pass the wall. Conversely,
   when set to ``transparent``, walls are not accounted for, and the tokamak
   appears to be transparent, effectively allowing radiation to pass through
   walls unaffected.

   The setting ``semi`` is a middle-ground, where only wall segments located at
   a radius less than the tokamak major radius are accounted for. This means
   that the tokamak central column is correctly accounted for, while the camera
   is allowed to be located outside the tokamak wall without radiation being
   blocked from it. This setting is a way of emulating diagnostic ports in which
   the radiation diagnostic may be somewhat retracted behind the regular tokamak
   wall boundary level.

