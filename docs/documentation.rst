
Documentation
=============

Running FLAME
---------------

FLAME consists of a single executable
named ``flame```, which will reside in the ``src``
directory upon successful compilation.
The main input file :ref:`flame_in.yaml <flame_in>` must be provided
in all FLAME runs and contains the simulation directives.
Most simulation tasks require additonal input files, 
like the input structure (commonly :ref:`posinp.yaml <structure_yaml>`), parameter files 
for atomic potentials, or run scripts to
call external codes.
