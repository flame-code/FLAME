.. _potential:

==================================
Potential
==================================

The interatomic potential parameters are set in this block.
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


        External codes that can be linked to FLAME:
            
            ``tinker``: TINKER molecular mechanics code :cite:`Tinker`.
    
            ``lammps``: LAMMPS molecular mechanics code :cite:`Lammps`.
            
            ``mopac``: MOPAC molecular mechanics code :cite:`Mopac`.
    
            ``dftb``: Density functional tight binding as implemented in dftb+ :cite:`Dftb`.
    
            ``bigdft``: BigDFT wavelet DFT code :cite:`Bigdft`.
    
            ``vasp``: Plane wave VASP code :cite:`Vasp`.
    
            ``abinit``: Plane wave Abinit code :cite:`Abinit`.
    
            ``espresso``: Plane wave Quantum Espresso code, pw.x :cite:`Espresso`.
    
            ``siesta``: Siesta DFT code :cite:`Siesta`.
    
            ``cp2k``: CP2K DFT/QM/MM code :cite:`Cp2k`.

            ``msock``: Network socket interface. Uses the i-Pi protocol to interact with extrenal codes :cite:`ceriotti_i-pi:_2014`.

..  warning:: mpmd seems identical to lj, which one to select?



**potential_sec**: (string) Secondary interatomic potential, usually used to perform a preliminary relaxation
or to evaluate an approximate Hessian matrix. All available options are identical to  **potential**.

    default: None

**core_rep**: (logical) If selected, a repulsive potential to the **potential** method is added based on the
atomic type. Some potentials (like PAW DFT) tend to cause issues if atoms get too close to each other
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
In units of the reciprocal lattice vectors, :math:`2\pi/\textrm{Bohr}`.  Recommended values are 
in the range of ``0.015`` and ``0.040`` for metals and insulators, respectively
Only relevant for periodic electronic structure codes. 

    default: ``[4.d-2, 6.d-2]``

``msock`` parameters
--------------------

**sockinet**: (integer) Selects Unix socket or internet (TCP) socket.

    default: ``0``

    options:
        
        ``0``: Unix socket

        ``1``: internet (TCP) socket

**sockport**: (integer) Socket port number.

   default: ``3141``

**sockhost**: (string) Socket address. If **sockinet** is ``0``, a string with the **sockhost** name will be
created in a temporary directory. Otherwise, a valid IP address must be provided (`127.0.0.1` for localhost).
    
    default: ``mh-driver``

**sockcutwf**: (list of two reals) Plane wave cutoff energies for the fine and coarse settings sent along 
with the i-Pi protocol. Only relevant for plane wave DFT codes that support this feature (like Quantum Espresso).

    default: ``[1.d0, 1.d0]``


``confine`` parameters
--------------------------
**confinement**: (logical) Determines if one or more 2D confinement potentials will be imposed based 
on polynomial functions. The general form of the potential 
is :math:`P = A(|e-\textbf{r}_i^\alpha|-r_c)^n`.
Where :math:`A` is the amplitude, :math:`e` is the equilibrium position along the
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

**av**: (integer) Method of defining the equilibrium position :math:`r_c` of the potential.
Given as a list of length **nconfine** if more than one confinement potential is imposed.

   default: ``2``
   
   options: 
      
      ``1``: The equilibrium position is set once during initiallization with respect to a predetermined value along the dimension :math:`\alpha` set in **dim**

      ``2``: The equilibrium position is set dynamically with respect to the average value of all involved atoms along the dimension :math:`\alpha` set in **dim**

**eq**: (real) Equilibrium position :math:`e_i` of the potential. 
Only relevant if **av** is set to ``1``.
The unit depends on the choice of **cartred**: Angstrom for ``C``, in reduced units if ``R``.
Given as a list of length **nconfine** if more than one confinement potential is imposed.

   default: ``0.d0``

**nat**: (integer) Number of atoms that are subjected to the potential.
Given as a list of length **nconfine** if more than one confinement potential is imposed.

   default: ``0``

**nat**: (list of integers and/or strings) The indices of the atoms that are subjected to the potential.
If all atoms are affected by the potential, the string ``all`` can be used instead of listing all atomic indices.
Given as a list of length **nconfine** (list of lists) if more than one confinement potential is imposed.

   default: ``all``

   options: 

      ``all``: all atoms aresubjected to the potential 

      ``[...]``: list of atomic indices


..  warning:: major part on the ann is still missing

..  warning:: major part on the electrostatic interactions is missing
