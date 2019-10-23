

Dummy title
========================

Welcome to FLAME
========================




FLAME -- a library of atomistic modeling environments --
is a software package to perform a wide range of atomistic simulations
for exploring the potential energy surfaces (PES) of complex condensed matter systems.
Available methods includes molecular dynamics simulations to sample free energy landscapes, 
saddle point searches to identify transition states, and gradient relaxations
to find dynamically stable geometries.
In addition to such common tasks, FLAME implements a structure prediction algorithm
based on the minima hopping method (MHM) to identify the ground state
structure of any system given solely the chemical composition, and a
framework to train a neural network potential to
reproduce the PES from *ab initio* calculations.
The combination of neural network potentials
with the MHM in FLAME allows highly
efficient and reliable identification of the ground state
as well as metastable structures  of molecules and crystals, 
as well as of nano structures, including surfaces, interfaces, 
and two-dimensional materials.




The FLAME Team
===============

The development of methods implemented in FLAME started back in 2007 in Basel, Switzerland, with the
implementation of various atomistic simulation tasks in 
different software packages. These tasks were later merged into two
main codes (Alborz and Minhocao), which were independently maintained and
developed by S. Alireza Ghasemi and Maximilian Amsler, respectively. 
In 2018, the FLAME code evolved as a streamlined integration 
of these two codes into a single package to 
facilitate atomistic simulation workflows.


The core developer team consists of:

*       S Alireza Ghasemi
*       Maximilian Amsler
*       Samare Rostami
*       Hossein Tahmasbi
*       Ehsan Rahmatizad
*       Somayeh Faraji
*       Robabe Rasoulkhani

Following contributors are greatefully acknowledged:

* Luigi Genovese
* Thomas Lenosky
* Li Zhu
* Miguel Marques

Licensing and Citation
========================

FLAME is distributed under the GPLv3 license. 

Please cite reference :cite:`amsler_flame_2019` when publishing results based on using FLAME,
which describes the main functionalities of the code.


Additionally, consider citing the following publications when using
specific portions of the FLAME code:

* Neural network potential (CENT): :cite:`ghasemi_interatomic_2015`
* Minima hopping method: :cite:`goedecker_minima_2004,amsler_2010_crystal,amsler_minima_2018`
* SQNM optimizer: :cite:`schaefer_stabilized_2015`
* Saddle point searches: :cite:`ghasemi_an_2011,schaefer_minima_2014`
* Structural fingerprints: :cite:`sadeghi_metrics_2013,li_fingerprint_2016`
* Electrostatic particle-particle, particle-density method: :cite:`ghasemi_particle_2007,Rostami2016`
* Interatimic potentials for silicon: :cite:`ghasemi_energy_2010`



It is the responsibility of the user to
follow the respective licensing terms
when FLAME is used in conjunction with external (quantum) engines.
