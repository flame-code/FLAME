.. _dynamics:

========
Dynamics
========

EXPLAIN Dynamics

dynamics options
==================

**nstep**: Number of molecular dynamics steps

   default: ``1``

**nfreq**: (integer) 

   default: ``0``

**dt**: (real) 

   default: ``No default value.``

**temp**: (real)

   default: 0.0``

**init_temp**: (real)

   default: 0.0``

**highest_freq**: (real)   

   default: ``1.d1``

**ntherm**: (integer)   

   default: ``2``

**md_method**: (string) 

   default: ``No default value.``

   options: 

**print_force**: (logical)

   default: ``False``

**restart**: (logical)

   default: ``False``


Parameters specifically for the MD escape trials in **task** ``minhocao``
------------------------------------------------------------------------------

**nmd**: (integer)  Number of molecular dynamics steps specifically for
the escape trials in the **task** ``minhocao``.

   default:  ``300``

**algo**: (integer)  Variable cell shape algorithm. 

   default: ``1``

   options:

      ``1``: Parrinello-Rahman

      ``2``: Cleveland

      ``3``: Wentzcovitch

      ``4``: Andersen (variable cell volume)

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

   default: ``0.d0``

**cellmass**: (real) Fictitious cell mass for all variable
cell shape MD tasks.

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
such that the totoal energy is conserved to within 1 %,
and has to be accompanied by setting **encon** to ``True``.
For peridic systems, 
``dt`` is aiondjusted such that the number of sampling points 
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

