module interface_code
  use global
  use void_lj_params
  use interface_core_repulsion
  use defs_basis
  !use cell_utils
  use interface_abinit
  use interface_cp2k
  use interface_mopac
  use interface_siesta
  use interface_vasp
  use interface_dftb
  use interface_lenosky_tb
  use interface_lenosky_meam
  use interface_blj
  use interface_mlj
  use interface_espresso  
#if defined(LAMMPS)
  use interface_lammps
#endif
  use interface_lenosky_tb_lj
#if defined(TINKER)
  use interface_tinker
#endif
!!! #if defined(ALBORZ)
  use interface_alborz
!!! #endif
  use interface_tersoff
  use interface_edip
  use interface_ipi
  use interface_msock
  use interface_lj_voids

  implicit none
  
  private
  public :: &
    get_energyandforces_single, &
    get_dos,                    &
    geopt_external

contains
 
  subroutine geopt_external(parini,latvec,xred,fcart,strten,energy,iprec,ka,kb,kc,counter)
    use mod_parini, only: typ_parini
    type(typ_parini), intent(in):: parini
    real(8), intent(in)  :: latvec(3,3)
    real(8), intent(in)  :: xred(3,nat)
    real(8), intent(out) :: fcart(3,nat)
    real(8), intent(out) :: strten(6)
    real(8), intent(out) :: energy
    real(8), intent(inout):: counter
    integer, intent(in)  :: iprec
    integer:: ka,kb,kc
    if(voids) then
      stop "Cannot run external geometry optimizer when using voids!"
    endif
    if(trim(parini%potential_potential)=="siesta") then
      call siesta_geopt(parini,latvec,xred,fcart,strten,energy,iprec,ka,kb,kc,counter)
    elseif(trim(parini%potential_potential)=="vasp") then
      call vasp_geopt(parini,latvec,xred,fcart,strten,energy,iprec,ka,kb,kc,counter)
    elseif(trim(parini%potential_potential)=="espresso") then
      call espresso_geopt(parini,latvec,xred,fcart,strten,energy,iprec,ka,kb,kc,counter)
    elseif(trim(parini%potential_potential)=="dftb") then
      call dftb_geopt(parini,latvec,xred,fcart,strten,energy,iprec,ka,kb,kc,counter)
    else
      stop "Code interface not yet implemented"
    endif
  end subroutine

  subroutine get_energyandforces_single(parini,latvec, xred, fcart, strten, energy, iprec, getwfk)
    use mod_parini, only: typ_parini
    implicit none
    type(typ_parini), intent(in):: parini
    real(8), intent(in)  :: latvec(3,3)
    real(8), intent(in)  :: xred(3,nat)
    real(8), intent(out) :: fcart(3,nat)
    real(8), intent(out) :: strten(6)
    real(8), intent(out) :: energy
    integer, intent(in)  :: iprec
    logical, intent(in)  :: getwfk


    real(8) :: dkpt,fcart_c(3,nat),energy_c,strten_c

!For the voids
    real(8):: energy_voidlj,strten_voidlj(6),fcart_voidlj(3,nat)
    integer:: nat_all, ntypat_all, l
!For the core repulsion
    real(8):: energy_rep,strten_rep(6),fcart_rep(3,nat)

    fcart=0.d0
    strten_voidlj=0.d0
    fcart_voidlj=0.d0

    if(iprec == 1) then
      dkpt = dkpt1
    else
      dkpt = dkpt2
    endif


!Setting up the k-point mesh
    if((trim(parini%potential_potential)=="siesta".and.siesta_kpt_mode.ne.2).or.&
       (trim(parini%potential_potential)=="vasp".and.vasp_kpt_mode.ne.2).or.&
       (trim(parini%potential_potential)=="abinit".and.abinit_kpt_mode.ne.2)) then
    write(*,'(a)') " # KPT mesh set up with kpt_mode /= 2"
    else
        if(dkpt.ne.0.d0) then
          call find_kpt(ka, kb, kc, latvec, dkpt)
          if(max_kpt) then
!if the ka1, kb1 and kc1 were reset or are complete garbage, set up new kpt mesh
            if(ka1.le.0) ka1=ka
            if(kb1.le.0) kb1=kb
            if(kc1.le.0) kc1=kc
            ka = max(min(ka1+1, ka), ka1)
            kb = max(min(kb1+1, kb), kb1)
            kc = max(min(kc1+1, kc), kc1)
          endif
          if(reuse_kpt.and.ka1.ne.0.and.kb1.ne.0.and.kc1.ne.0) then
            ka = ka1
            kb = kb1
            kc = kc1
          endif
        endif
        ka1 = ka
        kb1 = kb
        kc1 = kc
if(parini%verb.gt.0.and.trim(parini%potential_potential).ne."lammps") write(*,'(a,3(1x,i5))') " # KPT mesh set up as follows: ", ka, kb, kc
    endif
!Trigger exit if reuse and wrong kpt_options are used
    if((trim(parini%potential_potential)=="siesta".and.reuse_kpt.and.siesta_kpt_mode==1).or.&
       (trim(parini%potential_potential)=="vasp".and.reuse_kpt.and.vasp_kpt_mode==1).or.&
       (trim(parini%potential_potential)=="abinit".and.reuse_kpt.and.abinit_kpt_mode==1))&
       stop "Incompatible kpt-modes"

!Copy global array sizes to handle voids for single point energy evaluation
    if(voids) then
       !Backup the number of atoms and atom types
       nat_all=nat
       ntypat_all=ntypat
       !Replace those numbers with the ones used to get the atomic forces
       nat=nat_atoms
       ntypat=ntypat_atoms
    endif  
    
    if(trim(parini%potential_potential)=="abinit") then
      call make_input_abinit(latvec, xred, iprec, ka, kb, kc, getwfk)
    elseif(trim(parini%potential_potential)=="cp2k") then
      call make_input_cp2k(latvec, xred, iprec, ka, kb, kc, getwfk)
    elseif(trim(parini%potential_potential)=="siesta") then
      call make_input_siesta(latvec, xred, iprec, ka, kb, kc, getwfk)
    elseif(trim(parini%potential_potential)=="vasp") then
      call make_input_vasp  (latvec, xred, iprec, ka, kb, kc, getwfk)
    elseif(trim(parini%potential_potential)=="espresso") then
      call make_input_espresso  (latvec, xred, iprec, ka, kb, kc, getwfk)
    elseif(trim(parini%potential_potential)=="mopac") then
      call make_input_mopac (latvec, xred, iprec, ka, kb, kc, getwfk)
    elseif(trim(parini%potential_potential)=="dftb") then
      call make_input_dftb (latvec, xred, iprec, ka, kb, kc, getwfk)
    elseif(trim(parini%potential_potential)=="lenosky_tb") then
    elseif(trim(parini%potential_potential)=="lenosky_meam") then
    elseif(trim(parini%potential_potential)=="lenosky_tb_lj") then
!!! #if defined(ALBORZ)
    elseif(trim(parini%potential_potential)=="alborz") then
        icount_alborz=icount_alborz+1
        if(icount_alborz==1) then
            call call_to_alborz_init(nat)
        endif
!!! #endif
    elseif(trim(parini%potential_potential)=="blj") then
    elseif(trim(parini%potential_potential)=="mlj") then
    elseif(trim(parini%potential_potential)=="tersoff") then
    elseif(trim(parini%potential_potential)=="edip") then
    elseif(trim(parini%potential_potential)=="ipi") then
    elseif(trim(parini%potential_potential)=="msock") then
#if defined(LAMMPS)
    elseif(trim(parini%potential_potential)=="lammps") then
          if(count_lammps==0) call init_lammps(nat)
          count_lammps=count_lammps+1
#endif
#if defined(TINKER)
    elseif(trim(parini%potential_potential)=="tinker") then
          if(count_tinker==0) call init_tinker(nat,xred,latvec)
          count_tinker=count_tinker+1
#endif
    else
      stop "Code interface not yet implemented"
    endif
    
    if(trim(parini%potential_potential)=="lenosky_tb") then
      call lenosky_tb(latvec,xred,iprec,ka,kb,kc,fcart,energy,strten)
    elseif(trim(parini%potential_potential)=="lenosky_tb_lj") then
      call lenosky_tb_lj(latvec,xred,iprec,ka,kb,kc,fcart,energy,strten)
    elseif(trim(parini%potential_potential)=="lenosky_meam") then
      call lenosky_meam(latvec,xred,iprec,ka,kb,kc,fcart,energy,strten)
!!! #if defined(ALBORZ)
    elseif(trim(parini%potential_potential)=="alborz") then
        call call_to_alborz_get('bulk',nat,latvec,xred,fcart,energy,strten)
!!! #endif
    elseif(trim(parini%potential_potential)=="blj") then
      call blj(latvec,xred,fcart,strten,energy)
    elseif(trim(parini%potential_potential)=="mlj") then
      call mlj(latvec,xred,fcart,strten,energy)
    elseif(trim(parini%potential_potential)=="tersoff") then
      call tersoff(latvec,xred,fcart,strten,energy)
    elseif(trim(parini%potential_potential)=="edip") then
      call edip(latvec,xred,fcart,strten,energy)
    elseif(trim(parini%potential_potential)=="ipi") then
      call evaluate_ipi(latvec, xred, fcart, strten, energy, ka, kb, kc, iprec)
    elseif(trim(parini%potential_potential)=="msock") then
      call evaluate_msock(latvec, xred, fcart, strten, energy, ka, kb, kc, iprec)
#if defined(LAMMPS)
    elseif(trim(parini%potential_potential)=="lammps") then
      call call_lammps(latvec,xred,fcart,energy,strten)
#endif
#if defined(TINKER)
    elseif(trim(parini%potential_potential)=="tinker") then
      call tinker(latvec,xred,fcart,energy,strten)
#endif
    else
    !call system("sleep 1")
      call system("./runjob.sh")
    !call system("sleep 1")
    endif

    if(trim(parini%potential_potential)=="abinit") then
      call get_output_abinit(fcart, energy, strten)
    elseif(trim(parini%potential_potential)=="cp2k") then
      call get_output_cp2k(fcart, energy, strten)
    elseif(trim(parini%potential_potential)=="siesta") then
      call get_output_siesta(fcart, energy, strten)
    elseif(trim(parini%potential_potential)=="vasp") then
      call get_output_vasp  (fcart, energy, strten)
    elseif(trim(parini%potential_potential)=="espresso") then
      call get_output_espresso  (fcart, energy, strten)
    elseif(trim(parini%potential_potential)=="mopac") then
      call get_output_mopac (fcart, energy, strten)
    elseif(trim(parini%potential_potential)=="dftb") then
      call get_output_dftb (fcart, energy, strten)
    end if

!Copy back global array sizes and compte/add the LJ forces
    if(voids) then
       !Put back the number of atoms and atom types
       nat=nat_all
       ntypat=ntypat_all
       !Compute the LJ forces and stresses
       call lj_void_addon(latvec,xred,fcart_voidlj,strten_voidlj,energy_voidlj)
       energy=energy+energy_voidlj
       fcart=fcart+fcart_voidlj
       strten=strten+strten_voidlj
    endif  

 
!Add confinement potential
    if(confine.ge.2) then
      write(*,'(a,i3)') " # Running confinement with option :",confine
      call confinement_energy_forces(nat,xred,latvec,energy_c,fcart_c,strten_c)
      energy=energy+energy_c
      fcart=fcart+fcart_c
      strten=strten+strten_c
    endif

!Add core repulsion (no repulsion on LJ particles)
    if(core_rep) then
       call core_repulsion(latvec,xred,fcart_rep,strten_rep,energy_rep)
       energy=energy+energy_rep
       fcart=fcart+fcart_rep
       strten=strten+strten_rep
    endif

  end subroutine get_energyandforces_single


  !************************************************************************************
  subroutine get_dos(parini, latvec, xred, efermi, fdos, iprec, getwfk)
    use mod_parini, only: typ_parini
    implicit none
    type(typ_parini), intent(in):: parini
    real(8), intent(in)  :: latvec(3,3)
    real(8), intent(in)  :: xred(3,nat)
    real(8), intent(out) :: efermi
    real(8), intent(out) :: fdos
    integer, intent(in)  :: iprec
    logical, intent(in)  :: getwfk

    real(8) :: dkpt

    if(iprec == 1) then
      dkpt = dkpt1
    else
      dkpt = dkpt2
    endif


!Setting up the k-point mesh
    if((trim(parini%potential_potential)=="siesta".and.siesta_kpt_mode.ne.2).or.&
       (trim(parini%potential_potential)=="vasp".and.vasp_kpt_mode.ne.2).or.&
       (trim(parini%potential_potential)=="abinit".and.abinit_kpt_mode.ne.2)) then
    write(*,'(a)') " # KPT mesh set up with kpt_mode /= 2"
    else
        if(dkpt.ne.0.d0) then
          call find_kpt(ka, kb, kc, latvec, dkpt)
          if(max_kpt) then
            ka = max(min(ka1+1, ka), ka1)
            kb = max(min(kb1+1, kb), kb1)
            kc = max(min(kc1+1, kc), kc1)
          endif
          if(reuse_kpt.and.ka1.ne.0.and.kb1.ne.0.and.kc1.ne.0) then
            ka = ka1
            kb = kb1
            kc = kc1
          endif
        endif
        ka1 = ka
        kb1 = kb
        kc1 = kc
if(parini%verb.gt.0.and.trim(parini%potential_potential).ne."lammps") write(*,'(a,3(1x,i5))') " # KPT mesh set up as follows: ", ka, kb, kc
    endif
!Trigger exit if reuse and wrong kpt_options are used
    if((trim(parini%potential_potential)=="siesta".and.reuse_kpt.and.siesta_kpt_mode==1).or.&
       (trim(parini%potential_potential)=="vasp".and.reuse_kpt.and.vasp_kpt_mode==1).or.&
       (trim(parini%potential_potential)=="abinit".and.reuse_kpt.and.abinit_kpt_mode==1))&
       stop "Incompatible kpt-modes"


    if(trim(parini%potential_potential)=="abinit") then
      call make_input_abinit(latvec, xred, iprec, ka, kb, kc, getwfk, dos=.true.)
    elseif(trim(parini%potential_potential)=="vasp") then
      call make_input_vasp(latvec,xred,iprec,ka,kb,kc,getwfk, dos=.true.)
    elseif(trim(parini%potential_potential)=="dftb") then
      call make_input_dftb(latvec,xred,iprec,ka,kb,kc,getwfk, dos=.true.)     
    else
      stop "Code not yet implemented"
    endif

!    call system("sleep 1")
    call system("./runjob.sh")
!    call system("sleep 1")

    if(trim(parini%potential_potential)=="abinit") then
      call get_dos_abinit(fdos,efermi)
    elseif(trim(parini%potential_potential)=="vasp") then
      call get_dos_vasp(fdos,efermi)
    elseif(trim(parini%potential_potential)=="dftb") then
      call get_dos_dftb(fdos,efermi)
    endif

  end subroutine get_dos

end module interface_code
