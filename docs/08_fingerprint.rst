.. _fingerprint:

===========
Fingerprint
===========

The structural fingerprint parameters are set in this block.  
Some parameters are used by all methods, while others apply only to specific fingerprints.  
Note that the fingerprint methods discussed here are used in
conjunction with the ``minhocao`` task.


fingerprint parameters
=================================

general ``fingerprint`` parameters
------------------------------------------
**method**: (string) Method to evaluate the structural fingerprint.

    default: ``oganov``

    options:

       ``oganov``: Oganov's fingerprint method for crystalline systems based on radial distribution functions :cite:`oganov_how_2009`.

       ``gauss``:  Gaussian orbital overlap matrix method for clusters and molecules :cite:`sadeghi_metrics_2013`.
        
       ``gom``: Gaussian orbital overlap matrix method for crystalline system :cite:`li_fingerprint_2016`.
    
       ``molgom``: Gaussian orbital overlap matrix method for molecular crystals :cite:`li_fingerprint_2016`.

       ``coganov``: Continuous version of  Oganov's fingerprint, without histogram discretization along the radial direction.
    
       ``bcm``: Bond characterization matrix method for crystalline systems based on the Calypso method :cite:`wang_calypso_2012`.
    
       ``atorb``: Atomic orbitals overlap method for crystalline systems.
    
    
**rcut**: (real) Cutoff distance of the fingerprint, in units of Angstrom. Only relevant for periodic, crystalline systems.

    default: ``15.d0``

``oganov`` parameters
----------------------

**dbin**: (real) Bin size for the discretization of the Oganov fingerprint, in units of Angstrom.

   default: ``5.d-2``

``coganov`` parameters
---------------------------

**atnmax**: (integer) Maximum number of neighboring atoms to take into account. Only used to allocate temporary arrays.

   default: ``10000``

``oganov`` & ``coganov`` parameters
-------------------------------------

**sigma**: (real) Broadening of the Gaussian functions placed on the atomic coordinates, in units of Angstrom.

   default: ``2.d-2``

``bcm`` & ``atorb`` parameters
--------------------------------

**nl**: (integer) Maximum degree :math:`l` of the spherical harmonics to include.

   default: ``6``

``molgom`` parameters
----------------------

**principleev**: (integer) Number of principle eigenvectors for the molecular contraction.

   default: ``6``

**molecules**: (integer) Number of molecules in the system.

   default: ``1``

**expa**: (integer) Number of neighboring cells to be considered when computing the fingerprint. 

   default: ``1``

**molsphere**: (integer) Maximum number of molecules taken into account in the cutoff sphere.

   default: ``50``


**widthcut:**: (real) Characteristic length scale of the spherical molecular cutoff function. Arbitrary units, scaling factor.

   default: ``1.d0``


**widthover**: (real) Characteristic width of the molecular overlap orbitals, which will be scaled by the Van der Waals radii.

   default: ``1.d0``



``gom`` & ``molgom`` parameters
--------------------------------

**natx**: (integer) Maximum number of neighboring atoms allowed in the cutoff shell, per atom.

   default: ``75``

**nexcut**: (integer) Exponent of the spherical cutoff function.

   default: ``3``

**orbital**: (string) Degree of Gaussian type orbitals to include.

   default: ``S``

   options:

        ``S``: s-type orbitals

        ``P``: p-type orbitals
