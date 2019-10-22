.. _dynamics:

========
Dynamics
========

The details of molecular dynamics (MD) simulations are specified here. Not all
keywords apply to every available method. 
Note that there is a fundamental difference in using dynamics
to either sample free energies,
or as a means to escape local minima in the minima hopping method.
For the former, the goal is an accurate and unbiased sampling,
while for the latter it is the efficient escape from a catchment basin
of the PES.
Hence, the available parameters and their suggested values 
differ significantly between these two classes of dynamics.

dynamics options
==================

**md_method**: (string) Type of MD simulation ensemble.

   default: ``No default value.``

   options: 

        ``nvt_nose``:   NVT ensemble with the Nose-Hoover chain thermostat.

        ``nvt_langev``:  NVT ensemble with the Langevin thermostat.

        ``nve``:   NVE ensemble (under development).


**nmd**: (integer)  Number of MD steps.

   default:  ``300``

**nfreq**: (integer) The frequency at which detailed output is written to disk.

   default: ``0``

**dt**: (real)  Time step size. Units in fs.

   default: ``No default value.``

**temp**: (real)  Target temperature. In units of Kelvin.

   default: ``0.d0``

**init_temp**: (real) Initial temperature. In units of Kelvin.

   default: ``0.d0``

**highest_freq**: (real)  Coupling constant to the extended system. Units in :math:`{\textrm{Ha}}/{\textrm{THz}^{2}}`.

   default: ``1.d1``

**ntherm**: (integer) Number of thermostats used for the Nose-Hoover thermostat chain.

   default: ``2``

**restart**: (logical) Activates a restart from a previously interrupted MD run.

   default: ``False``

Parameters specifically for the MD escape trials in **task** ``minhocao``
------------------------------------------------------------------------------

**algo**: (integer)  Choice of variable cell shape algorithms. 

   default: ``1``

   options:

      ``1``: Parrinello-Rahman

      ``2``: Cleveland

      ``3``: Wentzcovitch

      ``4``: Andersen (variable cell volume only, fixed cell shape)

**integrator**: (integer) Integrator for the variable cell shape
MD equations of motion.

   default: ``3``

   options:

      ``1``: Velocity Verlet

      ``2``: Classic Verlet

      ``3``: Beeman corrector-predictor scheme

**presscomp**: (real) Sets a bias pressure for MD runs in **task** ``minhocao``
simulations. The high temperatures in the ``minhocao`` cycles
can lead to the cell expanding
significantly beyond the equlibrium volume. 
To counteract this thermal expansing,
a bias pressure can be imposed during the MD run.
In units of GPa.

   default: ``0.d0``

**cellmass**: (real) Fictitious cell mass for all variable
cell shape MD tasks. Arbitrary units.

   default: ``1.d0``

**mdmin_init**: (integer) The initial stopping criterion ``mdmin``
for an MD run in **task** ``minhocao``.
The MD run is terminated as soon as ``mdmin`` minima in the 
potential energy/enthalpy along the MD trajectory is encountered.

   default: ``2``

**auto_mdmin**: (logical) The stopping criterion for an MD run 
is adjusted dynamically if ``True``.
``mdmin`` will fluctuate between **mdmin_min** and **mdmin_max**.
Used for the MD escape steps in **task** ``minhocao``.

   default: ``False``

**mdmin_min**: (integer)
The minimal value of ``mdmin``.
Used for the MD escape steps in **task** ``minhocao`` and if **auto_mdmin** is ``True``.

   default: ``2``

**mdmin_max**: (integer)  
The maximal value of ``mdmin``.
Used for the MD escape steps in **task** ``minhocao`` and if **auto_mdmin** is ``True``.

   default: ``2``

**auto_mddt**: (logical) Activates the dynamical adjustment of the MD
timestep ``dt`` in **task** ``minhocao``. For clusters, ``dt`` is adjusted
such that the total energy is conserved to within 1 %,
and has to be accompanied by setting **encon** to ``True``.
For peridic systems, 
``dt`` is adjusted such that the number of sampling points 
per oscillation in the energy/enthalpy along an MD trajectory is close to the
target value **nit_per_min**.

   default: ``False``

**encon**: (logical)  Activates the dynamical adjustment of the MD
timestep ``dt`` in **task** ``minhocao`` for clusters
based on energy conservation.

   default: ``False``

**nit_per_min**: (integer) Target number of MD samples per
energy/enthalpy oscillation. Only used if **auto_mddt** is ``True``.

   default: ``25``

**dt_init**: (real) Initial MD time step ``dt``. In atomic units.                       

   default: ``2.d1``

