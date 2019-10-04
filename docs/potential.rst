.. _potential:

==================================
Potential
==================================

The potential parameters are set in this block.
Some parameters are used by all methods, while others
apply only to specific interatomic potentials.


potential parameters
=========================


general ``potential`` parameters
------------------------------------------

**potential**: (string) Method to evaluate the atomic interaction.

    default: ``No default value.``

    options:
        
        Methods implemented in FLAME:

            ``mpmd``: Lennard-Jones.
    
            ``lj``: Lennard-Jones.
    
            ``blj``: Binary Lennard-Jones.
    
            ``mlj``: Multicomponent Lennard-Jones with arbiratry number of species.
    
            ``ltb``: Lenosky's Tight-Binding method for silicon.
    
            ``edip``: Bazant's  Environment-Dependent Interatomic Potential for silicon.
    
            ``tersoff``: Tersoff's potential for silicon and sp3 carbon.
    
            ``ann``: Artificial neural network potential.
    
            ``coulomb``: Coulomb potential (electrostatic interations).

            ``msock``: Network socket interface. Uses the i-Pi protocol to interact with extrenal codes.

         External codes that can be linked to FLAME:
            
            ``tinker``: TINKER molecular mechanics code.
    
            ``lammps``: LAMMPS molecular mechanics code.
            
            ``mopac``: MOPAC molecular mechanics code.
    
            ``dftb``: Density functional tight binding as implemented in dftb+.
    
            ``bigdft``: BigDFT wavelet DFT code.
    
            ``vasp``: Plane wave VASP code.
    
            ``abinit``: Plane wave Abinit code.
    
            ``espresso``: Plane wave Quantum Espresso code as implemented in pw.x.
    
            ``siesta``: Siesta DFT code.
    
            ``cp2k``: CP2K DFT/QM/MM code.
    
**potential_sec**: (string) Secondary interatomic potential, usually used to perform a preliminary relaxations
or to evaluate an approximate Hessian matrix.

    default: None

    options: Identical to **potential**.

**core_rep**: (logical) If selected, adds a repulsive potential to the **potential** method based on the
atomic radii. Some potentials (like PAW DFT) tend to cause issues if atoms get too close to each other
and the core regions start to overlap. To avoid atoms from getting too close, a repulsive
:math:`\frac{1}{r^{12}}` term is added.

    default: ``False``

**kptmesh**: (list of three integers)
Desired k-points mesh. Will be overruled if **auto_kpt** is ``True``.
Only relevant for periodic electronic structure codes. 


    default: ``[1, 1, 1]``


**auto_kpt**: (logical) 
Activates a scheme to automatically compute the k-points mesh given a predefined
density. Only relevant for periodic electronic structure codes. 

    default: ``True``

**kptden**: (list of two reals)
Desired k-points density along every dimension for the fine and the coarse potential settings. 
In units of the reciprocal lattice vectors, :math:`2\pi/Bohr`.  Recommended value are 
in the range of ``0.015`` and ``0.04`` for metals and insulators, respectively
Only relevant for periodic electronic structure codes. 

    default: ``[4.d-2, 6.d-2]``

``msock`` parameters
--------------------

**sockinet**: (integer) Selects Unix socket if 0, and internet (TCP) socket otherwise.

   default: ``0``

**sockport**: (integer) Socket port number.

   default: ``3141``

**sockhost**: (string) Socket address. If **sockinet** is 0, a string with the sockhost name will be
created in a temporary directory. Otherwise, a valid IP addres must be provided (127.0.0.1 for localhost).
    
    default: ``mh-driver``

**sockcutwf**: (list of reals) Plane wave cutoff energies for the fine and coarse settings sent along 
with the i-PI protocol. Only relevant for plane wave DFT codes that support this feature (like Quantum Espresso).

    default: ``[1.d0, 1.d0]``


``confine`` parameters
--------------------------
**confinement**: (logical) Determines if one or more 2D confinement potentials will be imposed based on polynomic
functions. The general form of the potential is :math:`P = A(|e-\textbf{r}_i^\alpha|-r_c)^n`.
Where :math:`A` is the amplitude, :math:`e` is the eqilibrium position along the
dimension :math:`\alpha`, :math:`r_c` is the cutoff distance, 
and :math:`i` runs over all atoms that interact with the potential :math:`P`.

   default: ``False``


**nconfine**: (integer) Number of confinement potentials.

   default: ``0``

**cartred**: (string) Choice of Cartesian or reduced coordinates for setting up the confinement potential.
Given as a list of length **nconfine** if more than one confinement potential is imposed.

   default: ``C``

   options: 
      ``C``: Cartesian coordinates
      ``R``: Reduced coordinates

**dim**: (integer) Axis along which the confinement potential is applied.
Given as a list of length **nconfine** if more than one confinement potential is imposed.

   default: ``1``

   options: ``1``, ``2``, ``3`` for the x, y and z directions, respectively.

**exp**: (integer) Exponent *n* of the potential.
Given as a list of length **nconfine** if more than one confinement potential is imposed.

   default: ``4``

**prefac**: (real) Prefactor or the amplitude *A* of the potential, in units of eV.
Given as a list of length **nconfine** if more than one confinement potential is imposed.

   default: ``1.d-2``


**cut**: (real) Cutoff distance :math:`r_c` of the potential, in units of Angstrom.
Given as a list of length **nconfine** if more than one confinement potential is imposed.

   default: ``1.d0``

**av**: (integer) Method of defining the equilibrium position of the potential, :math:`r_c`.
Given as a list of length **nconfine** if more than one confinement potential is imposed.

   default: ``2``
   
   options: 
      
      ``1``: with respect to a predetermined value along the dimension :math:`\alpha` set in **dim**

      ``2``: with respect to the average value of all involved atoms along the dimension :math:`\alpha` set in **dim**

**eq**: (real) Equilibrium position :math:`e_i` of the potential. 
Only relevant if **av** is 1.
The unit depends on the choice of **cartred**: Angstrom for ``C``, in reduced units if ``R``
Given as a list of length **nconfine** if more than one confinement potential is imposed.

   default: ``0.d0``

**nat**: (integer) Number of atoms that are subjected to the potential.
Given as a list of length **nconfine** if more than one confinement potential is imposed.

   default: ``0``

**nat**: (list of integers, or strings) The indices of the atomix that are subjected to the potential.
If all atoms are affected by the potential, the string "all" can be used instead of listing all atomic indices.
Given as a list of length **nconfine** (list of lists) if more than one confinement potential is imposed.

   default: ``all``

   options: 

      ``all``: all atoms aresubjected to the potential 

      ``[...]``: list of atomic indices
