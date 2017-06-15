c
c     #############################################################
c     ##                                                         ##
c     ##  Tools to extract energies and forces from xyz feeds    ##
c     ##  Interface to be used with the Minhocao package         ##
c     ##                                                         ##
c     #############################################################
c
c
c     initialization routines called for the first time
      subroutine mhm_tinker_init(nat,maxbond,bc,ntypat,rxyz,
     &dist_ang,char_type,typat,bonds)
      implicit none
      integer,intent(in) :: nat
      integer,intent(in) :: bc
      character*120 xyzfile
      integer,intent(in) :: ntypat
      integer,intent(in) :: typat(nat)
      real(8),intent(in) :: rxyz(3,nat)
      real(8),intent(in) :: dist_ang(6)
      character(2),intent(in) :: char_type(nat)
      integer,intent(in):: maxbond
      integer, intent(in):: bonds(maxbond,nat)

      include 'sizes.i'
      include 'atoms.i'
      include 'files.i'
      include 'inform.i'
      include 'iounit.i'
      include 'bound.i'
      include 'boxes.i'
      include 'output.i'
      call initial
      xyzfile="tinker.xyz"
      call basefile (xyzfile)
      call suffix (xyzfile,'xyz','old')
      coordtype = 'CARTESIAN'
      call mhm_tinker_readxyz(nat,maxbond,ntypat,rxyz,char_type,typat,
     &bonds)
c
c     set the default values for the unitcell variables
c
      xbox = 0.0d0
      ybox = 0.0d0
      zbox = 0.0d0
      alpha = 0.0d0
      beta = 0.0d0
      gamma = 0.0d0
      orthogonal = .false.
      monoclinic = .false.
      triclinic = .false.
      octahedron = .false.
      spacegrp = '          '
      use_bounds=.false.
c
c     set the correct boundary conditions: 1 for periodic, 2 for free
c
      if(bc==1) then 
         use_bounds=.true.
         triclinic = .true.
         xbox=dist_ang(1)
         ybox=dist_ang(2)
         zbox=dist_ang(3)
         alpha=dist_ang(4)
         beta=dist_ang(5)
         gamma=dist_ang(6)
      elseif(bc==2) then
      else
         stop "Wround boundary condition in mhm_init of TINKER"
      endif
      use_replica = .false.
      call mhm_tinker_mechanic
      end subroutine

c     Update the cell and xyz coordinates, to be called prior to evaluation of energy and gradient
      subroutine mhm_tinker_update(in_xyz,in_xbox,in_ybox,in_zbox,
     &in_alpha,in_beta,in_gamma)
      implicit none
      integer::i
      real(8),intent(in)::in_xyz(3,n)
      real(8),intent(in)::in_xbox,in_ybox,in_zbox
      real(8),intent(in)::in_alpha,in_beta,in_gamma
      include 'sizes.i'
      include 'bound.i'
      include 'boxes.i'
      include 'atoms.i'
      include 'files.i'
      include 'inform.i'
      include 'iounit.i'
      include 'chgpot.i'
      include 'moment.i'
      include 'virial.i'
      include 'action.i'
      include 'cutoff.i'
      include 'warp.i'
      do i=1,n
         x(i)=in_xyz(1,i)
         y(i)=in_xyz(2,i)
         z(i)=in_xyz(3,i)
      enddo
      xbox=in_xbox
      ybox=in_ybox
      zbox=in_zbox
      alpha=in_alpha
      beta=in_beta
      gamma=in_gamma
      call lattice
      end subroutine

c     Computes the energy and forces from tinker, with already assigned number of atoms etc internally in tinker
      subroutine mhm_tinker_energyandforces(energy,derivs,virial)
      implicit none
      include 'sizes.i'
      include 'bound.i'
      include 'boxes.i'
      include 'atoms.i'
      include 'files.i'
      include 'inform.i'
      include 'iounit.i'
      include 'chgpot.i'
      include 'moment.i'
      include 'virial.i'
      include 'action.i'
      include 'cutoff.i'
      include 'warp.i'
      include 'cell.i'
      integer i,j,ixyz
      integer frame
      integer freeunit
      integer trimtext
      integer list(20)
      logical dosystem,doparam
      logical doenergy,doatom
      logical dolarge,dodetail
      logical doprops,doconect
      logical exist
      logical, allocatable :: active(:)
      character*1 letter
      character*120 record
      character*120 string
      character*120 xyzfile
      real*8:: energy,virial(3,3)
      real*8:: derivs(3,n)
      character*20 value
      call gradient (energy,derivs)
      virial=vir
      write(*,*) "TINKER_ENERGY", energy
      return
      end

c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ## mod subroutine getxyz  --  get Cartesian coordinate structure  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "getxyz" asks for a Cartesian coordinate file name,
c     then reads in the coordinates file
c
c
      subroutine mhm_tinker_getxyz
      implicit none
      include 'inform.i'
      include 'iounit.i'
      include 'output.i'
      integer ixyz
      integer freeunit
      logical exist
      character*120 xyzfile
c
c     override reading filename: read from a custom written tinker.xyz file, and use tinker.key and other extensions
c
      xyzfile="tinker.xyz"
      call basefile (xyzfile)
      call suffix (xyzfile,'xyz','old')
      inquire (file=xyzfile,exist=exist)
      if(.not.exist) then 
        write(*,*) "No tinker.xyz file found!!!"
        stop
      endif
c
c     first open and then read the Cartesian coordinates file
c
      coordtype = 'CARTESIAN'
      ixyz = freeunit ()
      open (unit=ixyz,file=xyzfile,status='old')
      rewind (unit=ixyz)
      call readxyz (ixyz)
      close (unit=ixyz)
c
c     quit if the Cartesian coordinates file contains no atoms
c
      if (abort) then
         write(*,*) "tinker.xyz is incorrect"
         call fatal
      end if
      return
      end subroutine
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  mod subroutine mechanic  --  initialize molecular mechanics  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "mechanic" sets up needed parameters for the potential energy
c     calculation and reads in many of the user selectable options
c
c
      subroutine mhm_tinker_mechanic
      implicit none
      include 'cutoff.i'
      include 'inform.i'
      include 'iounit.i'
      include 'potent.i'
      include 'vdwpot.i'
c
c
c     set the bonded connectivity lists and active atoms
c
      call attach
      call active
c
c     find bonds, angles, torsions, bitorsions and small rings
c
      call bonds
      call angles
      call torsions
      call bitors
      call rings
c
c     get the base force field from parameter file and keyfile
c
      call field
c
c     find unit cell type, lattice parameters and cutoff values
c
c      call unitcell
      call lattice
      call polymer
      call cutoffs
c
c     setup needed for potential energy smoothing methods
c
      call flatten
c
c     assign atom types, classes and other atomic information
c
      call katom
c
c     assign atoms to molecules and set the atom groups
c
      call molecule
      call cluster
c
c     find any pisystem atoms, bonds and torsional angles
c
      call orbital
c
c     assign bond, angle and cross term potential parameters
c
      if (use_bond .or. use_strbnd .or. use_strtor
     &    .or. (use_vdw .and. vdwtyp.eq.'MM3-HBOND'))  call kbond
      if (use_angle .or. use_strbnd .or. use_angang)  call kangle
      if (use_strbnd)  call kstrbnd
      if (use_urey)  call kurey
      if (use_angang)  call kangang
c
c     assign out-of-plane deformation potential parameters
c
      if (use_angle .or. use_opbend)  call kopbend
      if (use_angle .or. use_opdist)  call kopdist
      if (use_improp)  call kimprop
      if (use_imptor)  call kimptor
c
c     assign torsion and torsion cross term potential parameters
c
      if (use_tors .or. use_strtor .or. use_tortor)  call ktors
      if (use_pitors)  call kpitors
      if (use_strtor)  call kstrtor
      if (use_tortor)  call ktortor
c
c     assign van der Waals and electrostatic potential parameters
c
      if (use_vdw .or. use_solv)  call kvdw
      if (use_charge .or. use_chgdpl .or. use_solv)  call kcharge
      if (use_dipole .or. use_chgdpl)  call kdipole
      if (use_mpole .or. use_polar .or.
     &    use_solv .or. use_rxnfld)  call kmpole
      if (use_polar .or. use_solv)  call kpolar
      if (use_ewald)  call kewald
c
c     assign solvation, metal, pisystem and restraint parameters
c
      if (use_solv)  call ksolv
      if (use_metal)  call kmetal
      if (use_orbit)  call korbit
      if (use_geom)  call kgeom
      if (use_extra)  call kextra
c
c     set hybrid parameter values for free energy perturbation
c
      call mutate
c
c     quit if essential parameter information is missing
c
      if (abort) then
         write (iout,10)
   10    format (/,' MECHANIC  --  Some Required Potential Energy',
     &              ' Parameters are Undefined')
         call fatal
      end if
      return
      end
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  mod subroutine readxyz  --  input of Cartesian coordinates  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     Allocates and initiallizes all necessary variables passed to the subroutine, 
c     which are usually taken from an xyz-file
c
c
      subroutine mhm_tinker_readxyz(nat,maxbond,ntypat,rxyz,char_type,
     &typat,bonds)
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'couple.i'
      include 'files.i'
      include 'inform.i'
      include 'iounit.i'
      include 'titles.i'
      integer i,j,k,m
      integer ixyz,nmax
      integer next,size
      integer first,last
      integer nexttext
      integer trimtext
      integer, allocatable :: list(:)
      logical exist,opened
      logical quit,reorder
      logical clash
      character*120 xyzfile
      character*120 record
      character*120 string
      integer:: nat,ntypat,typat(nat)
      real(8):: rxyz(3,nat)
      character(2):: char_type(nat)
      integer:: maxbond, bonds(maxbond,nat)
c
c
c     initialize the total number of atoms in the system
c
      n = nat
c
c     extract the title and determine its length
c
         title = ' '
         ltitle = 0
c
c     initialize coordinates and connectivities for each atom
c
      do i = 1, n
         tag(i) = 0
         name(i) = '   '
         x(i) = 0.0d0
         y(i) = 0.0d0
         z(i) = 0.0d0
         type(i) = 0
         do j = 1, maxval
            i12(j,i) = 0
         end do
      end do
c
c     read the coordinates and connectivities for each atom
c
      do i = 1, n
         tag(i) = i
         name(i) = trim(char_type(i))
         x(i) = rxyz(1,i)
         y(i) = rxyz(2,i)
         z(i) = rxyz(3,i)
         type(i) = typat(i)
         i12(1:maxbond,i)=bonds(:,i)
      end do
c
c     for each atom, count and sort its attached atoms
c
      do i = 1, n
         n12(i) = 0
         do j = maxval, 1, -1
            if (i12(j,i) .ne. 0) then
               n12(i) = j
               goto 90
            end if
         end do
   90    continue
         call sort (n12(i),i12(1,i))
      end do
c
c     perform dynamic allocation of some local arrays
c
      nmax = 0
      do i = 1, n
         nmax = max(tag(i),nmax)
         do j = 1, n12(i)
            nmax = max(i12(j,i),nmax)
         end do
      end do
      allocate (list(nmax))
c
c     check for scrambled atom order and attempt to renumber
      reorder = .false.
      do i = 1, n
         list(tag(i)) = i
         if (tag(i) .ne. i)  reorder = .true.
      end do
      if (reorder) then
         write (iout,100)
  100    format (/,' READXYZ  --  Atom Labels not Sequential,',
     &              ' Attempting to Renumber')
         do i = 1, n
            tag(i) = i
            do j = 1, n12(i)
               i12(j,i) = list(i12(j,i))
            end do
            call sort (n12(i),i12(1,i))
         end do
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (list)
c
c     check for atom pairs with identical coordinates
c
      clash = .false.
c      if (n .le. 10000)  call chkxyz (clash)
c
c     make sure that all connectivities are bidirectional
c
      do i = 1, n
         do j = 1, n12(i)
            k = i12(j,i)
            do m = 1, n12(k)
               if (i12(m,k) .eq. i)  goto 120
            end do
            write (iout,110)  k,i
            write(*,*) "nat",n
  110       format (/,' READXYZ  --  Check Connection of Atom',
     &                 i6,' to Atom',i6)
            call fatal
  120       continue
         end do
      end do
      return
      end
