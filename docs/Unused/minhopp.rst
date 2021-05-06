.. _minhopp:

==================================
Minima hopping global optimization
==================================

The minima hopping (MH) method is a global optimization technique.
It employs a feedback mechanism which is based on the history
of visited minima. The feedback mechanism controls the
temperature (namely kinetic energy) in MD part.
The Monte Carlo algorithm used in the minima hopping uses information of
the history as well and a different feedback mechanism controls the
acceptance/rejection of the metropolis algorithm.

List of available options in minhopp

minhopp options
=================

**nstep**: Number of Monte Carlo Steps

    default: ``0``

**nsoften**: Initial velocity vector in MD part is chosen randomly.
The performance of the MH method can be enhanced by eliminating the
component of the velocity along hard mode of the Hessian of the potential.
This is done in an algorithm similar to the dimer method used for saddle point
search methods. The process is called softening and it is performed iteratively.
nstep determines the number of iterations in softening.

    default: ``7``

**eref**: Reference energy, minhopp will stop as soon as finds a minimum
as low as *eref*.

    default: ``-10**50``

**alpha1**: The feedback parameter in the minima hopping method for adjusting ediff.
It reduces ediff whenever a minimum is accepted. alpha1 must be smaller than one.

    default: ``1.0/1.02``

**alpha2**: The feedback parameter in the minima hopping method for adjusting ediff.
It increases ediff whenever a minimum is rejected. alpha1 must be larger than one.

    default: ``1.02d0``

**beta1**: The feedback parameter in the minima hopping method for adjusting
kinetic energy. It increases the kinetic energy whenever an MD escape fails.
beta1 must be larger than one.

    default: ``1.05``

**beta2**:  The feedback parameter in the minima hopping method for adjusting
kinetic energy. It increases the kinetic energy whenever an MD escape succeeds
but ends up into a minimum which is already visited.
beta2 must be larger than one.

    default: ``1.05``

**beta3**: The feedback parameter in the minima hopping method for adjusting
kinetic energy. It reduces the kinetic energy whenever an MD escape succeeds
and ends up into a new minimum.
beta3 must be smaller than one.

    default: ``1.0/1.05``

**trajectory**: is a logical variable. If it is set to .true., minhopp prints
info files all configurations in different parts of the algorithm like
molecular dynamics and optimization, softening.

    default: ``.false.``

**print_force**: If .true., atomic forces are printed for
every configuration written into files.

    default: ``.false.``

**etoler**: The tolerance criterion to check similarity and dissimilarity
between two configurations used in the minhopp.
Since there is a similarity check in FLAME, this option
and the depending codes will be removed and those of
:ref:`conf_comp <conf_comp>` will be employed.

    default: ``1.d-2``

**nrandoff**: An integer based on which minhopp will conducts different
MPI processes, i.e. minima hopping threads, to proceed with distinct pathways.
This is done by letting different random number to be used.

    default: ``0``

**minter**: It defines how frequent minhop results are written into files.
These files can be used for a restart minhopp run. For MH runs in which
the potential is obtained from fast methods such as force fields, it is
recommended to use a much larger value than 1, otherwise FLAME will be
slow due to the printing intermediate results many times.

    default: ``1``

**mdmin**: This parameter determines the number of minima along the
MD pathway before quitting MD escape.
These minima should not be confused with minima of the energy landscape.

    default: ``3``

