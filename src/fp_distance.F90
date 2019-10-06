subroutine get_fp_distance(parini,fp_len,fp1,fp2,fp_dist)
use mod_parini, only: typ_parini
!This routine will initiallize the parameters for the fingerprinting
!For 10<fp_method<20: fully periodic systems
!For 20<fp_method<30: molecular systems
use fingerprint
use defs_basis, only: pi
implicit none
type(typ_parini), intent(in):: parini
integer:: fp_len
real(8):: fp(fp_len),pos_red(3,parini%nat),latvec(3,3),rxyz(3,parini%nat),fp1(fp_len),fp2(fp_len),fp_dist
select case(fp_method)
  case(11)!Oganov method
        call get_cosinedistance(parini,fp1,fp2,fp_11_fp_size,fp_11_fp_dim,parini%ntypat_global,fp_11_nkinds_sum,fp_dist)
  case(12)!Calypso method
        call get_distance_calypso(parini,fp1,fp2,fp_12_fp_dim,fp_12_nl,fp_dist)
  case(13)!Modified Calypso method
        call get_distance_malypso(fp1,fp2,fp_13_fp_dim,parini%nat,parini%typat_global,fp_13_nl,fp_dist)
  case(14)!xyz2sm
        call get_distance_xyz2sm(parini%nat,parini%typat_global,fp_len/parini%nat,fp1,fp2,fp_dist)
  case(15)!Continuous Oganov
        call get_cosinedistance_coganov(fp1,fp2,fp_15_fp_size,fp_15_fp_dim,parini%ntypat_global,fp_15_nkinds_sum,fp_15_rcut,pi,fp_dist)
  case(16)!Continuous Atomic Oganov
        call get_cosinedistance_coganov_atomic(fp1,fp2,parini%nat,fp_16_fp_size,fp_16_fp_dim,parini%typat_global,parini%ntypat_global,fp_16_nkinds_sum,fp_16_rcut,pi,fp_dist)
  case(17)!GOM
        call get_distance_gauss(fp1, fp2, parini%fp_17_lseg, parini%nat, parini%ntypat_global, parini%typat_global, fp_dist)
  case(18)!MOLGOM
        call get_distance_molgom(fp1,fp2,fp_dist,parini%fp_18_lseg,parini%fp_18_molecules,parini%fp_18_molecules_sphere,parini%fp_18_principleev)
  case(21)!Gaussian molecular overlap
        call fpdistance_gaussmol(fp_len,fp1,fp2,fp_dist)
  case default
     stop "Wrong choice for FP"
end select
end subroutine
