.. _potential:

==================================
Potential
==================================

EXPLAIN 

List available potentials:

potential options
====================

**potential**:

    default: None, must be set

    options:

        ``mpmd``: Lennard-Jones

        ``lj``: Lennard-Jones

        ``ltb``: EXPLAIN

        ``ann``: EXPLAIN

        ``bigdft``: EXPLAIN

        ``vasp``: EXPLAIN

        ``coulomb``: EXPLAIN

        ``siesta``: EXPLAIN

**potential_sec**:

    default: None

    options: Identical to ***potential***.



confine options
--------------------
**confinement**: Determines if one or more 2D confinement potentials will be imposed based on polynomic
functions (boolean). The general form of the potential is :math:`P = A(|e-\textbf{r}_i^\alpha|-r_c)^n`.
Where :math:`A` is the amplitude, :math:`e` is the eqilibrium position along the
dimension :math:`\alpha`, :math:`r_c` is the cutoff distance, 
and :math:`i` runs over all atoms that interact with the potential :math:`P`.

   default: ``False``


**nconfine**: Number of confinement potentials (integer). 

   default: ``0``

**cartred**: Choice of Cartesian or reduced coordinates for setting up the confinement potential (string).
Given as a list of length **nconfine** if more than one confinement potential is imposed.

   default: ``C``

   options: 
      ``C``: Cartesian coordinates
      ``R``: Reduced coordinates

**dim**: Axis along which the confinement potential is applied (integer).
Given as a list of length **nconfine** if more than one confinement potential is imposed.

   default: ``1``

   options: ``1``, ``2``, ``3`` for the x, y and z directions, respectively.

**exp**: Exponent *n* of the potential (integer).
Given as a list of length **nconfine** if more than one confinement potential is imposed.

   default: ``4``

**prefac**: Prefactor or the amplitude *A* of the potential, in units of eV (real).
Given as a list of length **nconfine** if more than one confinement potential is imposed.

   default: ``1.d-2``


**cut**: Cutoff distance :math:`r_c` of the potential, in units of Angstrom (real).
Given as a list of length **nconfine** if more than one confinement potential is imposed.

   default: ``1.d0``

**av**: Method of defining the equilibrium position of the potential, :math:`r_c` (integer).
Given as a list of length **nconfine** if more than one confinement potential is imposed.

   default: ``2``
   
   options: 
      
      ``1``: with respect to a predetermined value along the dimension :math:`\alpha` set in **dim**

      ``2``: with respect to the average value of all involved atoms along the dimension :math:`\alpha` set in **dim**

**eq**: Equilibrium position :math:`e_i` of the potential (real). Only relevant if **av** is 1.
The unit depends on the choice of **cartred**: Angstrom for ``C``, in reduced units if ``R``
Given as a list of length **nconfine** if more than one confinement potential is imposed.

   default: ``0.d0``

**nat**: Number of atoms that are subjected to the potential (integer).
Given as a list of length **nconfine** if more than one confinement potential is imposed.

   default: ``0``

**nat**: The indices of the atomix that are subjected to the potential (list of integers, or string).
If all atoms are affected by the potential, the string "all" can be used instead of listing all atomic indices.
Given as a list of length **nconfine** (list of lists) if more than one confinement potential is imposed.

   default: ``all``

   options: 

      ``all``: all atoms aresubjected to the potential 

      ``[...]``: list of atomic indices
