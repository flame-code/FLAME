subroutine init_fp(parini,fp_len,latvec)
!This routine will initiallize the parameters for the fingerprinting
!For 10<fp_method<20: fully periodic systems
!For 20<fp_method<30: molecular systems
use mod_parini, only: typ_parini
use fingerprint
use global, only: units
use defs_basis, only: Bohr_Ang,pi
implicit none
type(typ_parini), intent(in):: parini
integer:: fp_len,iat,nmax
real(8):: convert,latvec(3,3),vol
!Convert cutoff, sigma and dbin into atomic units
convert=1.d0
if(units=="angstroem") then
  convert=1.d0/Bohr_Ang
endif

!Get recomended size of the continuous oganov FP
if(fp_method==15.or.fp_method==16) then
   call getvol(latvec,vol)
   call estimate_nmax_per_atom(vol,parini%nat,parini%ntypat_global,parini%fp_rcut*convert,pi,nmax)
   write(*,*) vol,parini%fp_rcut*convert,nmax,parini%nat,pi,parini%ntypat_global
   write(*,'(a,i10)') " # Estimated fingerprint size for COGANOV and CAOGANOV: ",nmax
   if(nmax.gt.parini%fp_at_nmax) write(*,'(a,i10,i10)') " # WARNING: FPATNMAX too small!", parini%fp_at_nmax, nmax
endif

select case(fp_method)
  case(11)!Oganov method
!     read(56,*) fp_11_rcut,fp_11_sigma,fp_11_dbin
     fp_11_rcut=parini%fp_rcut*convert
     fp_11_sigma=parini%fp_sigma*convert
     fp_11_dbin=parini%fp_dbin*convert
     fp_11_fp_size=ceiling(fp_11_rcut/fp_11_dbin)
     fp_11_fp_dim=parini%ntypat_global*(parini%ntypat_global+1)/2
     fp_len=fp_11_fp_size*fp_11_fp_dim
     if(.not.allocated(fp_11_nkinds_sum)) allocate(fp_11_nkinds_sum(parini%ntypat_global))
     fp_11_nkinds_sum=0
     do iat=1,parini%nat
       fp_11_nkinds_sum(parini%typat_global(iat))=fp_11_nkinds_sum(parini%typat_global(iat))+1
     enddo
  case(12)!Calypso method
!TEMPORARY READ FROM STDINPUT
!     read(*,*) tmpvar,fp_12_nl
!     read(56,*) tmpvar,fp_12_nl
     fp_12_nl=parini%fp_nl
     fp_12_fp_dim=parini%ntypat_global*(parini%ntypat_global+1)/2
     fp_len=fp_12_fp_dim*fp_12_nl
     if(.not.allocated(fp_12_r_cut)) allocate(fp_12_r_cut(fp_12_fp_dim))
     fp_12_r_cut=parini%fp_rcut*convert
  case(13)!Modified Calypso method
!TEMPORARY READ FROM STDINPUT
!     read(*,*) tmpvar,fp_13_nl
!     read(56,*) tmpvar,fp_13_nl
     fp_12_nl=parini%fp_nl
     fp_13_fp_dim=parini%ntypat_global!*(ntypat+1)/2
     fp_len=fp_13_nl*parini%nat*fp_13_fp_dim
     if(.not.allocated(fp_13_r_cut)) allocate(fp_13_r_cut(fp_13_fp_dim))
     fp_13_r_cut=parini%fp_rcut*convert
  case(14)!xyz2sm
     fp_len=3*parini%fp_14_m*parini%nat
  case(15)!Continuous Oganov method
     fp_15_rcut=parini%fp_rcut*convert
     fp_15_sigma=parini%fp_sigma*convert
     fp_15_fp_size=parini%fp_at_nmax
     fp_15_fp_dim=parini%ntypat_global*(parini%ntypat_global+1)/2
     fp_len=3*fp_15_fp_size*fp_15_fp_dim
!Careful: the FP has 3 entries for the gaussians
     if(.not.allocated(fp_15_nkinds_sum)) allocate(fp_15_nkinds_sum(parini%ntypat_global))
     fp_15_nkinds_sum=0
     do iat=1,parini%nat
       fp_15_nkinds_sum(parini%typat_global(iat))=fp_15_nkinds_sum(parini%typat_global(iat))+1
     enddo
  case(16)!Continuous Atomic Oganov method
     fp_16_rcut=parini%fp_rcut*convert
     fp_16_sigma=parini%fp_sigma*convert
     fp_16_fp_size=parini%fp_at_nmax
     fp_16_fp_dim=parini%ntypat_global*(parini%ntypat_global+1)/2
     fp_len=3*fp_16_fp_size*parini%ntypat_global*parini%nat
!Careful: the FP has 3 entries for the gaussians, nmax entries for all possible neighbors, ndim possible AB interaction, nat atomic lists of gaussians
     if(.not.allocated(fp_16_nkinds_sum)) allocate(fp_16_nkinds_sum(parini%ntypat_global))
     fp_16_nkinds_sum=0
     do iat=1,parini%nat
       fp_16_nkinds_sum(parini%typat_global(iat))=fp_16_nkinds_sum(parini%typat_global(iat))+1
     enddo
  case(17)
     fp_len=parini%fp_17_lseg*(parini%ntypat_global+1)*parini%nat
  case(18) !Molecular gaussian orbital fingerprint
     fp_len=parini%fp_18_lseg*parini%fp_18_molecules_sphere*parini%fp_18_principleev*parini%fp_18_molecules
  case(21)!Gaussian molecular overlap
!The method only has a FP of length nat
     fp_len=parini%nat  !If we only have stype orbitals, alse fp_len=4*nat
  case default
     stop "Wrong choice for FP"
end select
close(56)
end subroutine

!**********************************************************************************************

subroutine get_fp(parini,fp_len,pos_red,latvec,fp)
!This routine will initiallize the parameters for the fingerprinting
!For 10<fp_method<20: fully periodic systems
!For 20<fp_method<30: molecular systems
use mod_parini, only: typ_parini
use fingerprint, only: fp_15_fp_size, fp_method, fp_11_rcut, fp_11_sigma, fp_11_dbin
use fingerprint, only: fp_11_fp_size, fp_11_nkinds_sum, fp_11_fp_dim
use fingerprint, only: fp_12_r_cut, fp_12_fp_dim, fp_16_fp_size, fp_12_nl, fp_13_nl
use fingerprint, only: fp_13_r_cut, fp_16_fp_dim
use fingerprint
use yaml_output
implicit none
type(typ_parini), intent(in):: parini
integer:: fp_len,iat,natmol
real(8):: fp(fp_len),pos_red(3,parini%nat),latvec(3,3),rxyz(3,parini%nat),vol,rcov_arr(parini%nat),fp_coganov_atomic(3,fp_15_fp_size,parini%ntypat_global,parini%nat)
real(8):: rvan(parini%nat) !nat*molecules)
character(len=2):: finalchar(parini%nat) ! dimension(nat*molecules)

select case(fp_method)
  case(11)!Oganov method
     call getvol(latvec,vol)
     call yaml_map('Suggested minimal cutoff radius for Oganov FP',vol**(1.d0/3.d0)*2.d0)
     !write(*,'(a,es15.7)') " # Suggested minimal cutoff radius for Oganov FP: ", vol**(1.d0/3.d0)*2.d0
     call rxyz_int2cart(latvec,pos_red,rxyz,parini%nat)
     call get_fp_oganov(parini%nat,rxyz,latvec,fp_11_rcut,fp_11_sigma,fp_11_dbin,&
          &parini%typat_global,parini%ntypat_global,fp_11_nkinds_sum,fp_11_fp_size,fp_11_fp_dim,fp)
  case(12)!Calypso method
     call rxyz_int2cart(latvec,pos_red,rxyz,parini%nat)
     call get_fp_calypso(parini%nat,rxyz,latvec,fp_12_r_cut,parini%typat_global,parini%ntypat_global,fp_12_fp_dim,fp_12_nl,fp)
  case(13)!Modified Calypso method
     call rxyz_int2cart(latvec,pos_red,rxyz,parini%nat)
     call get_fp_malypso(parini%nat,rxyz,parini%rcov,latvec,fp_13_r_cut,parini%typat_global,parini%ntypat_global,fp_13_fp_dim,fp_13_nl,fp)
  case(14)!xyz2sm fingerprint
     do iat=1,parini%nat
        rcov_arr(iat)=parini%rcov(parini%typat_global(iat))
     enddo
     call rxyz_int2cart(latvec,pos_red,rxyz,parini%nat)
     call xyz2sm(parini%nat,latvec,rxyz,rcov_arr,parini%fp_14_w1,parini%fp_14_w2,parini%fp_14_m,fp)
  case(15)!C-Oganov fingerprint
     call rxyz_int2cart(latvec,pos_red,rxyz,parini%nat)
     call get_fp_coganov(parini%nat,rxyz,latvec,fp_15_rcut,fp_15_sigma,parini%rcov,&
          &parini%typat_global,parini%ntypat_global,fp_15_nkinds_sum,fp_15_fp_size,fp_15_fp_dim,fp)
  case(16)!C-Atomic-Oganov fingerprint
     call rxyz_int2cart(latvec,pos_red,rxyz,parini%nat)
     call get_fp_coganov_atomic(parini%nat,rxyz,latvec,fp_16_rcut,fp_16_sigma,parini%rcov,&
          &parini%typat_global,parini%ntypat_global,fp_16_nkinds_sum,fp_16_fp_size,fp_16_fp_dim,fp)
  case(17)!GOM
     do iat = 1, parini%nat
        rcov_arr(iat) = parini%rcov(parini%typat_global(iat))
     end do
     call rxyz_int2cart(latvec,pos_red,rxyz,parini%nat)
     call get_fp_gauss(parini%nat, parini%ntypat_global, parini%fp_17_natx_sphere, parini%typat_global, parini%fp_17_lseg, parini%fp_17_width_cutoff,&
          & parini%fp_17_nex_cutoff, latvec, rxyz, rcov_arr, fp)
  case(18)!MOLGOM
!This fingerprint wants to have the number of atoms per molecule
     natmol=parini%nat/parini%fp_18_molecules     
     if(natmol*parini%fp_18_molecules.ne.parini%nat) stop "Something wrong with the number of molecules"
     call findmolecule(parini,rxyz,latvec,finalchar,pos_red,parini%char_type,parini%typat_global,parini%ntypat_global,natmol)
!
! Assign the Van-der-Waals radii
!
     do iat = 1, parini%nat
         call sym2rvan(finalchar(iat), rvan(iat))
     end do
!
!from this system a fingerprint is taken
!
     call periodic_fingerprint(parini,rxyz,latvec,finalchar,rvan,fp,natmol)

  case(21)!Gaussian molecular overlap
     do iat=1,parini%nat
        rcov_arr(iat)=parini%rcov(parini%typat_global(iat))
     enddo
     call rxyz_int2cart(latvec,pos_red,rxyz,parini%nat)
     call fingerprint_gaussmol(parini%nat,fp_len,rxyz,rcov_arr,fp)
  case default
     stop "Wrong choice for FP"
end select
end subroutine
