.. _minhocao:

=================================================
Minima hopping for crystal structure optimization
=================================================

The minima hopping method implements a 
a global optimization technique to identify the ground state structure of
any chemical system. 

Input/Output files
==================
**poscur.ascii/poscur.vasp**: starting or current configuration of a MHM run
either in the *ascii* or *vasp* format, in Angstrom. Will be 
overwritten continuously as the simulation progresses with the 
current local minimum structure.


**ioput**: contains parameter values which are read and continuously updated throughout an MHM run.
In one line, 3 real numbers must be provided, which correspond to
*ediff* (in Ha/cell), *temperature* (in K and scaled based on the atomic masses), 
*maximum temperature*  (in K and scaled based on the atomic masses).

   Suggested initial choice is ``1.d-2 5.d2 1.d5`` (for a system with around 16 atoms)

**earr.dat**: contains informations about the initial/current state of the MHM run, 
and is populated with structural data during the simulation.
The first two lines must be always provided:

* Line 1: *nmin_cur*,  *nmin_max*
* Line 2: *delta_enthalpy*, *delta_fingerprint*


*nmin_cur* is the currently known number of local minima,
and *nmin_max* is the maximal number of minima to be found in
the next MHM run. *nmin_max* must always be larger than 
*nmin_max*. Every new MHM run must start with *nmin_cur: 0*

*delta_enthalpy*  (units of Ha/cell) and *delta_enthalpy* are the 
minimal enthalpy and structural fingerprint tolerances, respectively,
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
The structures correspoding to the *nmin_cur* local minima are written into 
separate files, where *XXXXX* corresponds to the index *imin* in **earr.dat**.


minhocao options
=================

**nsoften**: Number of softening iterations to eliminate the
high-frequency modes from the initial MD velocities (integer).

    default: ``7``

**alpha_at**: Stepsize for the softening algorithm on the atomic degrees of freedom (real). Arbitrary units.

    default: ``5.d-1``

**alpha_lat**: Stepsize for the softening algorithm on the lattice degrees of freedom (real). Arbitrary units.

    default: ``5.d-1``

**auto_soften**: Choice of automatically adjusting the stepsize of the softening iterations based on 
a gradient feedback (boolean). 

    default: ``True``

**eref**: Reference energy, minhopp will stop as soon as finds a minimum
as low as *eref* (real). In units of Ha/cell.

    default: ``-1.d50``

**alpha1**: Feedback parameter to adjust *ediff* (real).
Reduces *ediff* whenever a minimum is accepted. alpha1 must be smaller than one.

    default: ``1.0/1.02``

**alpha2**: Feedback parameter to adjust *ediff* (real).
Increases *ediff* whenever a minimum is rejected. alpha1 must be larger than one.

    default: ``1.02d0``

**beta1**: Feedback parameter to adjust the kinetic energy of the MD escape trials (real).
Increases the kinetic energy whenever an MD escape trial fails.
beta1 must be larger than one.

    default: ``1.05``

**beta2**: Feedback parameter to adjust the kinetic energy of the MD escape trials (real).
Increases the kintetic energy whenever an MD escape 
step leads to a known structure.
beta2 must be larger than one.

    default: ``1.05``

**beta3**: Feedback parameter to adjust the kinetic energy of the MD escape trials (real).
Reduces the kintetic energy whenever an MD escape succeeds
and ends up into a new mimimum.
beta3 must be smaller than one.

    default: ``1.0/1.05``

