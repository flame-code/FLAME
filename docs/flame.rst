


Welcome to FLAME
========================




FLAME is a software package to perform a wide range of atomistic simulations
for exploring the potential energy surfaces (PES) of complex condensed matter systems.
The range of methods include molecular dynamics simulations to sample free energy landscapes, 
saddle point searches to identify transition states, and gradient relaxations
to find dynamically stable geometries.
In addition to such common tasks, FLAME implements a structure prediction algorithm
based on the minima hopping method (MHM) to identify the ground state
structure of any system given solely the chemical composition, and a
framework to train a neural network potential to
reproduce the PES from *ab initio* calculations.
The combination of neural network potentials
with the MHM in FLAME allows a highly
efficient and reliable identification of the ground state
as well as metastable structures  of molecules and crystals, 
as well as of nano structures, including surfaces, interfaces, 
and two-dimensional materials.


The FLAME team
===============

The FLAME software is available on the flame website http://flame-code.org, and the
latest development version can be found on GitHub https://github.com/flame-code.


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


Please cite reference :cite:`amsler_flame_2019` when publishing results based on using FLAME.


Additionally, consider citing the following publications when using
specific portions of the FLAME code:

* Neural network potential (CENT): :cite:`ghasemi_interatomic_2015`
* Minima hopping method: :cite:`goedecker_minima_2004,amsler_2010_crystal,amsler_minima_2018`.
* SQNM optimizer: :cite:`schaefer_stabilized_2015`.
* Saddle point searches: :cite:`ghasemi_an_2011,schaefer_minima_2014`.

