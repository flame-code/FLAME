.. _potential:

==================================
Potential
==================================

The interatomic potential parameters are set in this block.
Some parameters are used by all methods, while others
apply only to specific interatomic potentials.



potential parameters
==================================

general ``potential`` parameters
------------------------------------------

**potential**: (string) Method to evaluate the atomic interaction.

    default: ``No default value.``

    options:
        
        Methods implemented in FLAME:

            ``lj``: Lennard-Jones.
    
            ``blj``: Binary Lennard-Jones.
    
            ``mlj``: Multicomponent Lennard-Jones with arbitrary number of species.
    
            ``ltb``: Lenosky's Tight-Binding method for silicon :cite:`ghasemi_energy_2010`.
    
            ``edip``: Bazant's  Environment-Dependent Interatomic Potential for silicon :cite:`ghasemi_energy_2010`.
    
            ``tersoff``: Tersoff's potential for silicon and sp3 carbon :cite:`ghasemi_energy_2010`.
    
            ``ann``: Artificial neural network potential.
    
            ``coulomb``: Coulomb potential (electrostatic interactions).


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

            ``msock``: Network socket interface. Uses the i-Pi protocol to interact with external codes :cite:`ceriotti_i-pi_2014`.

**potential_sec**: (string) Secondary interatomic potential, usually used to perform a preliminary relaxation
or to evaluate an approximate Hessian matrix. All available options are identical to  **potential**.

    default: None

**core_rep**: (logical) If selected, a repulsive potential to the **potential** method is added based on the
atomic type. Some potentials (like PAW DFT) tend to cause issues if atoms get too close to each other
and the core regions start to overlap. To avoid atoms from getting too close, a repulsive
:math:`\frac{1}{r^{12}}` term is added.

    default: ``False``

**kptmesh**: (list of three integers)
Desired k-points mesh. It will be overruled if **auto_kpt** is ``True``.
Only relevant for periodic electronic structure codes. 


    default: ``[1, 1, 1]``


**auto_kpt**: (logical) 
Activates a scheme to automatically compute the k-points mesh given a predefined
density. Only relevant for periodic electronic structure codes. 

    default: ``True``

**kptden**: (list of two reals)
Desired k-points density along every dimension for the fine and the coarse potential settings. 
In units of the reciprocal lattice vectors, :math:`2\pi/\textrm{Bohr}`.  Recommended values are 
in the range of ``0.015`` and ``0.040`` for metals and insulators, respectively.
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

**confine**: One or more 2D confinement potentials can be imposed based 
on polynomial functions. The general form of the potential 
is :math:`P = A(|e-\textbf{r}_i^\alpha|-r_c)^n`.
Where :math:`A` is the amplitude, :math:`e` is the equilibrium position along the
dimension :math:`\alpha`, :math:`r_c` is the cutoff distance, 
and :math:`i` runs over all atoms that interact with the potential :math:`P`.


   **confinement**: (logical) Determines if one or more 2D confinement potentials will be imposed.
   
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
   
   **av**: (integer) Method of defining the equilibrium position :math:`e` of the potential.
   Given as a list of length **nconfine** if more than one confinement potential is imposed.
   
      default: ``2``
      
      options: 
         
         ``1``: The equilibrium position is set once during initialization with respect to a predetermined value along the dimension :math:`\alpha` set in **dim**
   
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
   
         ``all``: all atoms are subjected to the potential 
   
         ``[...]``: list of atomic indices
   
``ewald`` parameters
--------------------------
   
**ewald**: If electrostatics is a part of the interactions in any FLAME potential, e.g., 
in the CENT potential, then the ``ewald`` key can be used to set the relevant parameters.

    **ewald**: (logical) This subkey determines whether the Ewald method is invoked.
    If ``True``, the Ewald approach is used up
    to speed up the calculations, especially when the calculations involve localized charge densities,
    e.g., when the Gaussian width of atomic charge densities in CENT are small.

        default: ``False``

    **psolver**: (string) Determines the method for the Poisson solver.

        default: ``No default value.``
                
        options: 
                
            ``p3d``: The P3D method is used, applicable only for slab boundary conditions.

            ``kwald``: Fourier summation, applicable only in the CENT potential and for bulk boundary condition.

            ``bigdft``: The BigDFT PSolver is invoked if FLAME is linked with the BigDFT PSolver.
            Currently, only applicable for bulk and free boundary conditions.

    **cell_ortho**: (logical) Activates efficient subroutines to place Gaussian
    charge densities on the grid. ``True`` can be used only when the simulation cell
    is orthogonal and the type of simulation does not change the cell variables.
    If ``False``, then generic subroutines are called to put Gaussian charge densities on the grid.

        default: ``False``

    **ecut**: (real) The cutoff energy that specifies how dense the basis set is when solving the Poisson's equation.
    The value is used for every non-pairwise method available in FLAME. There is no default value
    and it must be set. Units in Ha.

        default: ``No default value.``

    **ecutz**: (real) The cutoff energy that specifies how dense the basis set is in the *z*-direction when solving the Poisson's equation.
    The value is used only when the ``p3d`` method is selected in **psolver**. There is no default value and it must be set.
    Units in Ha.

        default: ``No default value.``

    **rgcut**: (real) The cutoff radius beyond which the atomic Gaussian charge densities are assumed to
    vanish. This parameter is  *not* the actual cutoff radius but is a unitless parameter that
    is multiplied by the Gaussian width value. There is no default value and it must be set.
    Typically, ``6.0`` is a reasonable choice, and for very high accuracy one may use values
    up to ``9.0``. Arbitrary units.

        default: ``No default value.``

    **bias_type** (string) Select schemes to treat external fields or special boundary conditions
    if the ``p3d`` method is used.

        default: ``no``
                
        options: 

            ``p3d_bias``:  For modeling conditions in which an ionic material is confined between
                           two parallel metallic plates at two different electric potentials.  
                           For more information see :cite:`Rostami2016`.

            ``fixed_efield``: An external uniform electric field is applied along the non-periodic *z*-direction.

    **plane_voltageu** (real) Voltage of the upper plate used if ``p3d_bias`` is selected for  **bias_type**. Units in Ha/e. 
        default: ``0.d0``

    **plane_voltagel** (real) Voltage of the lower plate used if ``p3d_bias`` is selected for  **bias_type**. Units in Ha/e.
        default: ``0.d0``

    **external_field** (real) Uniform electric field used if ``fixed_efield`` is selected for **bias_type**. Units in Ha/e/Bohr.
        default: ``0.d0``

