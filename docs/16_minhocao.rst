.. _minhocao:

=================================================
Minima Hopping for Crystal Structure Optimization
=================================================

The minima hopping method implements a 
global optimization technique to identify the ground state structure of
any chemical system. 

Input/Output files
==================

In addition to the standard *flame_in.yaml* file, further input files are
required to start a MHM run. Also, important information
will be written into output files other than *flame_log.yaml*, e.g., 
*earr.dat* and *global.mon*.



**poscur.ascii/poscur.vasp**: intitial or current configuration of a MHM run
either in the *ascii* or *vasp* format, in Angstrom. It will be 
overwritten continuously as the simulation progresses with the 
current local minimum structure.


**ioput**: contains parameter values which are read and continuously updated throughout an MHM run.
In one line, 3 real numbers must be provided, which correspond to
*ediff* (in Ha/cell), *temperature* (in K and scaled based on the atomic masses ``amass``), 
*maximum temperature*  (in K and scaled based on the atomic masses ``amass``).

   Suggested initial choice is ``1.d-2 5.d2 1.d5`` (for a system with around 16 atoms)

**earr.dat**: contains information about the initial/current state of the MHM run, 
and is populated with structural data during the simulation.
The first two lines must be always provided:

* Line 1: *nmin_cur*,  *nmin_max*
* Line 2: *delta_enthalpy*, *delta_fingerprint*


*nmin_cur* is the currently known number of local minima,
and *nmin_max* is the maximal number of minima to be found in
the next MHM run. *nmin_max* must always be larger than 
*nmin_cur*. Every new MHM run must start with *nmin_cur: 0*
In a restart run, *nmin_cur* corresponds to the number of 
local minima found in the previous run(s).


*delta_enthalpy*  (units of Ha/cell) and *delta_fingerprint* are the 
minimal enthalpies and structural fingerprint tolerances, respectively,
below which two structures will be considered identical.

   Suggested initial choice: ``2.d-3 1.d-1`` (for a system with around 16 atoms and the
   Oganov fingerprint method)

During an MHM run, the value *nmin_cur* will incrementally increase,
and addtional *nmin_cur* lines will be written to **earr.dat**,
each corresponding to a *distinct* local minimum.
This list of local minima is constantly sorted based on the enthalpy.
Every line consists of 6 entries: 
   
   * Line i + 2: *imin*, *enthalpy*, *energy*, *nvisits*, *spg*, *spg_tol*

where *imin* is the index of the local minimum, *enthalpy* and *energy* are
its enthalpy/energy (in Ha/cell),
*nvisits* is the number of times it has been encountered during the run,
*spg* is its space group index (for crystalline systems), and  *spg_tol*
is the tolerance used in SPGLIB to determine the space group.


**poslowXXXXX.ascii/poslowXXXXX.vasp**:
The structures corresponding to the *nmin_cur* local minima are written into 
separate files, where *XXXXX* corresponds to the index *imin* in **earr.dat**.


``minhocao`` options
======================

**nsoften**: (integer) Number of softening iterations to eliminate the
high-frequency modes from the initial MD velocities.

    default: ``7``

**alpha_at**: (real) Stepsize for the softening algorithm on the atomic degrees of freedom. Arbitrary units.

    default: ``5.d-1``

**alpha_lat**: (real) Stepsize for the softening algorithm on the lattice degrees of freedom. Arbitrary units.

    default: ``5.d-1``

**auto_soften**: (logical) Choice of automatically adjusting the stepsize of the softening iterations based on 
gradient feedback. 

    default: ``True``

**eref**: (real) Reference energy, MHM run will stop as soon a minimum
with an energy/enthlapy less than *eref* is found. In units of Ha/cell.

    default: ``-1.d50``

..  **alpha1**: (real) Feedback parameter to adjust *ediff*.
..  Reduces *ediff* whenever a minimum is accepted. alpha1 must be smaller than one.
..  
..      default: ``1.d0/1.02d0``
..  
..  **alpha2**: (real) Feedback parameter to adjust *ediff*.
..  Increases *ediff* whenever a minimum is rejected. alpha1 must be larger than one.
..  
..      default: ``1.02d0``
..  
..  **beta1**: (real) Feedback parameter to adjust the kinetic energy of the MD escape trials.
..  Increases the kinetic energy whenever an MD escape trial fails.
..  beta1 must be larger than one.
..  
..      default: ``1.05d0``
..  
..  **beta2**: (real) Feedback parameter to adjust the kinetic energy of the MD escape trials.
..  Increases the kintetic energy whenever an MD escape 
..  step leads to a known structure.
..  beta2 must be larger than one.
..  
..      default: ``1.05d0``
..  
..  **beta3**: (real) Feedback parameter to adjust the kinetic energy of the MD escape trials.
..  Reduces the kintetic energy whenever an MD escape succeeds
..  and ends up into a new minimum.
..  beta3 must be smaller than one.
..  
..      default: ``1.d0/1.05d0``



