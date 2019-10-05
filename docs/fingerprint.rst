.. _fingerprint:

===========
Fingerprint
===========

The structural fingerprint parameters are set in this block.  Some parameters are used by all methods, while others apply only to specific fingerprints.  


fingerprint parameters
=================================

general ``fingerprint`` parameters
------------------------------------------
**method**: (string) Method to evaluate the structural fingerprint.

    default: ``oganov``

    options:

       ``gauss``:  Gaussion orbital overlap matrix method for clusters and molecules.
        
       ``gom``: Gaussion orbital overlap matrix method for crystalline system.
    
       ``molgom``: Gaussion orbital overlap matrix method for molecular crystals.

       ``oganov``: Oganov's fingerprint method for crystalline system based on radial distribution functions.

       ``coganov``: Continuous version of  Oganov's fingerprint, without histogram discretization along the radial direction.
    
       ``bcm``: Bond characterization matrix method for crystalline system based on the Calypso code.
    
       ``atorb``: Atomic orbitals overlap method for crystalline system.
    
    
**rcut**: (real) Cutoff distance of the fingerprints, in units of Angstrom.  Only relevant for periodic, crystalline systems.

    default: ``15.d0``






``oganov`` parameters
----------------------

**dbin**: (real) Bin size for the discretization of the Oganov fingperprint, in units of Angstrom.

   default: ``5.d-2``

``coganov`` parameters
---------------------------

**atnmax**: (integer) Maximum number of neighboring atoms to consider.

   default: ``10000``

``oganov`` & ``coganov`` parameters
-------------------------------------

**sigma**: (real) Broadening of the Gaussian functions placed on the atomic coordinates, in units of Angstrom.

   default: ``2.d-2``

``bcm`` & ``atorb`` parameters
--------------------------------

**nl**: (integer) Maximum degree :math:`l` of the spherical harmics to include.

   default: ``6``


``molgom`` parameters
----------------------

**principleev**: (integer) Number of principle eigenvectors for the molecular contraction.

   default: ``6``

**molecules**: (integer) Number of molecules in the system.

   default: ``1``

**expa**: (integer) Number of neighoring cells to be considered when computing the fingerprint. 

   default: ``1``

**molsphere**: (integer) Maximum number of molecules taken into account in the cutoff sphere.

   default: ``50``


**widthcut:**: (real) Characteristic length scale of the spherical molecular cutoff function.

   default: ``1.d0``


**widthover**: (real) Characteristic width of the molecular overlap orbitals, which will be scaled by the van der Vaals radii.

   default: ``1.d0``



``gom`` & ``molgom`` parameters
--------------------------------

**natx**: (integer) Maximum number of neighboring atoms allowed in the cutoff shell per atom.

   default: ``75``

**nexcut**: (integer) Exponent of the spherical cutoff function.

   default: ``3``

**orbital**: (string) Degree of Gaussian type orbitals to include.

   default: ``S``

   options:

        ``S``: s-type orbitals

        ``P``: p-type orbitals
