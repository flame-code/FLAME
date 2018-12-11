   module callback
     use mod_parini, only: typ_parini
     use mod_atoms, only: typ_atoms, set_rat
     implicit none
     type(typ_parini):: parini_lammps
     type(typ_atoms):: atoms
#if defined(HAVE_LAMMPS)
     contains
       subroutine copy_parini_for_lammps(parini)
       type(typ_parini):: parini
       parini_lammps=parini
       end subroutine copy_parini_for_lammps
       subroutine fortran_callback(lmp, timestep, nlocal, ids, c_pos, c_fext) & 
     & bind(C, name='f_callback')
       use, intrinsic :: ISO_C_binding
       use LAMMPS
       implicit none
       type (C_ptr), value :: lmp
       integer(C_int64_t), intent(in), value :: timestep
       integer(C_int), intent(in), value :: nlocal
       real (C_double), dimension(:,:), pointer :: x
       type(c_ptr) :: c_pos, c_fext, c_ids
       double precision, pointer :: fext(:,:), pos(:,:)
       integer, intent(in) :: ids(nlocal)
       real(C_double) :: virial(6)
       real (C_double) :: etot
       integer :: natom , i, iat
       double precision lx, ly, lz
       real (C_double), pointer :: boxxlo, boxxhi
       real (C_double), pointer :: boxylo, boxyhi
       real (C_double), pointer :: boxzlo, boxzhi
       real (C_double), pointer :: boxxy, boxxz, boxyz
       double precision :: volume
       type (C_ptr) :: Cptr
       type (C_ptr), pointer, dimension(:) :: Catom

       ! nlocal is the number of atoms
       ! pos array is the position of atoms
       call c_f_pointer(c_pos, pos, [3,nlocal])
       call c_f_pointer(c_fext, fext, [3,nlocal])
       call lammps_extract_global(boxxlo, lmp, 'boxxlo')
       call lammps_extract_global(boxxhi, lmp, 'boxxhi')
       call lammps_extract_global(boxylo, lmp, 'boxylo')
       call lammps_extract_global(boxyhi, lmp, 'boxyhi')
       call lammps_extract_global(boxzlo, lmp, 'boxzlo')
       call lammps_extract_global(boxzhi, lmp, 'boxzhi')
       call lammps_extract_global(boxxy, lmp, 'xy')
       call lammps_extract_global(boxxz, lmp, 'xz')
       call lammps_extract_global(boxyz, lmp, 'yz')
       if(boxxlo/=0.d0) stop 'ERROR: boxxlo/=0.d0'
       if(boxylo/=0.d0) stop 'ERROR: boxylo/=0.d0'
       if(boxzlo/=0.d0) stop 'ERROR: boxzlo/=0.d0'
       lx = boxxhi - boxxlo
       ly = boxyhi - boxylo
       lz = boxzhi - boxzlo
       !volume = lx*ly*lz
       !call LJ(pos,fext,nlocal,etot)
       if(atoms%nat/=nlocal) then
           write(*,'(a,2i7)') 'atoms%nat/=nlocal',atoms%nat,nlocal
           stop
       endif
       !back origin of box to (0,0,0)
       !atoms%cellvec(1,1)= boxxhi-boxxlo   !ax
       !atoms%cellvec(2,1)=-boxylo          !ay
       !atoms%cellvec(3,1)=-boxzlo          !az
       !atoms%cellvec(1,2)=-boxxlo          !bx
       !atoms%cellvec(2,2)= boxyhi-boxylo   !by
       !atoms%cellvec(3,2)=-boxzlo          !bz
       !atoms%cellvec(1,3)=-boxxlo          !cx
       !atoms%cellvec(2,3)=-boxylo          !cy
       !atoms%cellvec(3,3)= boxzhi-boxzlo   !cz
       !do iat=1,atoms%nat
       !    atoms%rat(1,iat)=pos(1,iat)-boxxlo
       !    atoms%rat(2,iat)=pos(2,iat)-boxylo
       !    atoms%rat(3,iat)=pos(3,iat)-boxzlo
       !enddo
       atoms%cellvec(1,1) = boxxhi
       atoms%cellvec(1,2) = boxxy 
       atoms%cellvec(2,2) = boxyhi
       atoms%cellvec(1,3) = boxxz 
       atoms%cellvec(2,3) = boxyz 
       atoms%cellvec(3,3) = boxzhi
       call getvol_alborz(atoms%cellvec,volume)
       call set_rat(atoms,pos,setall=.true.)
       call cal_potential_forces(parini_lammps,atoms)
       etot=atoms%epot
       do iat=1,atoms%nat
           fext(1,iat)=atoms%fat(1,iat)
           fext(2,iat)=atoms%fat(2,iat)
           fext(3,iat)=atoms%fat(3,iat)
       enddo
       !The unit of stress in FLAME is Ha/bohr^3
       ! 1Ha/bohr^3 = 29421.02648438959 GPa 
       ! virial_LAMMPS = -stress_FLAME*convert_to_Pascals
       virial(1) = atoms%stress(1,1)*-1.d0/volume !*29421.02648438959d9
       virial(2) = atoms%stress(2,2)*-1.d0/volume !*29421.02648438959d9
       virial(3) = atoms%stress(3,3)*-1.d0/volume !*29421.02648438959d9
       virial(4) = atoms%stress(1,2)*-1.d0/volume !*29421.02648438959d9
       virial(5) = atoms%stress(1,3)*-1.d0/volume !*29421.02648438959d9
       virial(6) = atoms%stress(2,3)*-1.d0/volume !*29421.02648438959d9
       !write(*,'(a,f20.10)') 'VOL',volume
       !write(*,*) 'STRESSS',virial(1)
       !write(*,'(a,3es24.15)') 'STRESS ',atoms%stress(1,1),atoms%stress(2,1),atoms%stress(3,1)
       !write(*,'(a,3es24.15)') 'STRESS ',atoms%stress(1,2),atoms%stress(2,2),atoms%stress(3,2)
       !write(*,'(a,3es24.15)') 'STRESS ',atoms%stress(1,3),atoms%stress(2,3),atoms%stress(3,3)
       call lammps_set_user_energy (lmp, etot)
       call lammps_set_user_virial (lmp, virial)
       end  subroutine
#endif
    end module callback
