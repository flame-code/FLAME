   module callback
     use mod_parini, only: typ_parini
     use mod_atoms, only: typ_atoms
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
       !real(C_double) :: virial(6)
       real (C_double) :: etot
       real(C_double), pointer :: ts_lmp
       double precision :: stress(3,3), ts_dftb
       integer :: natom , i, iat
       !real (C_double), parameter :: econv = 627.4947284155114 ! converts from Ha to
       !double precision, parameter :: fconv = 1185.793095983065 ! converts from Ha/bohr to
       !double precision, parameter :: autoatm = 2.9037166638E8
       double precision lx, ly, lz
       real (C_double), pointer :: boxxlo, boxxhi
       real (C_double), pointer :: boxylo, boxyhi
       real (C_double), pointer :: boxzlo, boxzhi
       real (C_double), pointer :: boxxy, boxxz, boxyz
       !double precision, parameter :: nktv2p = 68568.4149999999935972
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
       volume = lx*ly*lz
       !open (unit = 10, status = 'replace', action = 'write', file='lammps.gen')
       !write(10,*)nlocal,"S"
       !write(10,*) "C"
       !do i = 1, nlocal
       !  write(10,'(2I,3F15.6)')i,1,pos(:,ids(i))
       !enddo
       !write(10,*)"0.0 0.0 0.0"
       !write(10,*)lx,0,0
       !write(10,*)0,ly,0
       !write(10,*)0,0,lz
       !close(10)
       !call system("./dftb+ > dftb.out")
       !open (unit = 10, status = 'old', file = 'results.out')
       !read(10,*)etot
       !read(10,*)ts_dftb
       !do i = 1, 3
       !  read(10,*)stress(i,:)
       !enddo
       !stress (:,:) = stress(:,:)*autoatm
       !virial(1) = stress(1,1)/(nktv2p/volume)
       !virial(2) = stress(2,2)/(nktv2p/volume)
       !virial(3) = stress(3,3)/(nktv2p/volume)
       !virial(4) = stress(1,2)/(nktv2p/volume)
       !virial(5) = stress(1,3)/(nktv2p/volume)
       !virial(6) = stress(2,3)/(nktv2p/volume)
       !etot = etot*econv
       !call lammps_set_external_vector(lmp,1,ts_dftb*econv)
       !call LJ(pos,fext,nlocal,etot)
       if(atoms%nat/=nlocal) stop 'ERROR: atoms%nat/=nlocal'
       write(*,'(a,6es24.15)') 'CELLVEC ',atoms%cellvec(1,1),boxxhi
       write(*,'(a,6es24.15)') 'CELLVEC ',atoms%cellvec(1,2),boxxy
       write(*,'(a,6es24.15)') 'CELLVEC ',atoms%cellvec(2,2),boxyhi
       write(*,'(a,6es24.15)') 'CELLVEC ',atoms%cellvec(1,3),boxxz
       write(*,'(a,6es24.15)') 'CELLVEC ',atoms%cellvec(2,3),boxyz
       write(*,'(a,6es24.15)') 'CELLVEC ',atoms%cellvec(3,3),boxzhi
       !stop 'STOPPED AT CELLVEC'
       if(abs(atoms%cellvec(1,1)-boxxhi)>1.d-15*boxxhi) stop 'ERROR: inconsistency in boxxhi'
       if(abs(atoms%cellvec(1,2)-boxxy )>1.d-12       ) stop 'ERROR: inconsistency in boxxy '
       if(abs(atoms%cellvec(2,2)-boxyhi)>1.d-15*boxyhi) stop 'ERROR: inconsistency in boxyhi'
       if(abs(atoms%cellvec(1,3)-boxxz )>1.d-12       ) stop 'ERROR: inconsistency in boxxz '
       if(abs(atoms%cellvec(2,3)-boxyz )>1.d-12       ) stop 'ERROR: inconsistency in boxyz '
       if(abs(atoms%cellvec(3,3)-boxzhi)>1.d-15*boxzhi) stop 'ERROR: inconsistency in boxzhi'
       do iat=1,atoms%nat
           atoms%rat(1,iat)=pos(1,iat)
           atoms%rat(2,iat)=pos(2,iat)
           atoms%rat(3,iat)=pos(3,iat)
           !write(31,*) pos(1,iat),pos(2,iat),pos(3,iat)
           !atoms%rat(1,iat)=pos(0,iat-1)
           !atoms%rat(2,iat)=pos(1,iat-1)
           !atoms%rat(3,iat)=pos(2,iat-1)
       enddo
       call cal_potential_forces(parini_lammps,atoms)
       etot=atoms%epot
       do iat=1,atoms%nat
           fext(1,iat)=atoms%fat(1,iat)
           fext(2,iat)=atoms%fat(2,iat)
           fext(3,iat)=atoms%fat(3,iat)
           !write(41,'(3f20.14)') fext(1,iat),fext(2,iat),fext(3,iat)
           !fext(0,iat-1)=atoms%fat(1,iat)
           !fext(1,iat-1)=atoms%fat(2,iat)
           !fext(2,iat-1)=atoms%fat(3,iat)
       enddo
       !do i = 1, nlocal
         !read(10,*)fext(:,ids(i))
       !  fext(:,ids(i)) = fext(:,ids(i))*fconv
       !enddo
       !close(10)
       call lammps_set_user_energy (lmp, etot)
       !call lammps_set_user_virial (lmp, virial)

       end  subroutine
#endif
    end module callback

