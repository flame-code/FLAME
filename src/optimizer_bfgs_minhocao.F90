module mbfgs_interface
!Module used for the geometry optimizer BFGS with ABINIT LINESEARCH
interface
   !subroutine unit_matrix(mat,ndim)
   !  implicit none
   !  real(8),DIMENSION(:,:), INTENT(INOUT) :: mat
   !  integer:: ndim
   !end subroutine unit_matrix

   function vabs(v) result(res)
     implicit none
     real(8),dimension(:):: v
     real(8):: res
   end function vabs

   function outerprod(a,b)
     real(8),dimension(:),intent(in)::a,b
     real(8),dimension(size(a),size(b))::outerprod
   end function outerprod
end interface
end module mbfgs_interface


!>   Geometry optimization, parametrisation routine.
subroutine geopt_init()
  use minpar
  implicit none

  parmin_bfgs%approach  = 'unknown'
  parmin_bfgs%iter      = 0
  parmin_bfgs%iflag     = 0
  parmin_bfgs%verbosity = 1
  parmin_bfgs%MSAVE=7
  parmin_bfgs%MP=16
  parmin_bfgs%LP=16
  parmin_bfgs%MAXFEV=10
  parmin_bfgs%GTOL=9.d-1
  parmin_bfgs%XTOL=1.d-15
  parmin_bfgs%FTOL=1.d-6
  parmin_bfgs%STPMIN=1.d-20
  parmin_bfgs%STPMAX=20.d0
  parmin_bfgs%DIAGCO=.FALSE.
  parmin_bfgs%IWRITE=.FALSE.

END SUBROUTINE geopt_init




!> @file
!!  Routines to do BFGS geometry optimisation
!! @author
!!    Copyright (C) 2007-2011 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS

subroutine GEOPT_RBFGS_MHM(parini,parres,latvec_in,xred_in,fcart_in,strten_in,etot_in,iprec,counter,folder)
 use minpar
 use mod_parini, only: typ_parini
implicit none
type(typ_parini), intent(in):: parini
type(typ_parini), intent(inout):: parres
integer:: iprec,iat,i,istr
real(8):: latvec_in(3,3),xred_in(3,parini%nat),fcart_in(3,parini%nat),etot_in,strten_in(6),counter,flat(9)
real(8):: fmax_tol,fmax_at,fmax_lat,fmax,pressure,strtarget(6),dstr(6),enthalpy
logical:: fail,getwfk
character(40):: folder
parmin_bfgs%maxiter_lat=20
fmax_tol=1.d10
counter=0.d0
getwfk=.false.
call get_BFGS_forces_lattice(parini,parres,xred_in,flat,latvec_in,enthalpy,getwfk,iprec,latvec_in,xred_in,etot_in,fcart_in,strten_in)
do 
!!Compute maximal component of forces, EXCLUDING any fixed components
 call get_fmax(parini,fcart_in,strten_in,fmax,fmax_at,fmax_lat)
 fmax_tol=min((fmax_at+fmax_lat)/2.d0*0.3d0,fmax_tol*0.5d0)
 fmax_tol=max(fmax_tol,parini%paropt_geopt%fmaxtol*1.d-2)
 write(*,*) "# New tolarance",fmax_tol
 if(fmax.lt.parini%paropt_geopt%fmaxtol.or.int(counter).gt.parini%paropt_geopt%nit) exit
 call bfgs_driver_atoms(parini,parres,latvec_in,xred_in,fcart_in,strten_in,etot_in,iprec,counter,fmax_tol,folder)
 call get_fmax(parini,fcart_in,strten_in,fmax,fmax_at,fmax_lat)
 if(fmax.lt.parini%paropt_geopt%fmaxtol.or.int(counter).gt.parini%paropt_geopt%nit) exit
 fail=.false.
 call lbfgs_driver_lattice(parini,parres,latvec_in,xred_in,fcart_in,strten_in,etot_in,iprec,counter,fail,fmax_tol,folder)
 if(fail) write(*,*) "# LBFGS GEOPT FAILED, switching to backup routine by ALIREZA"
 if(fmax.lt.parini%paropt_geopt%fmaxtol.or.int(counter).gt.parini%paropt_geopt%nit) exit
 if(fail) call bfgs_driver_lattice(parini,parres,latvec_in,xred_in,fcart_in,strten_in,etot_in,iprec,counter,fmax_tol,folder)
enddo

call get_enthalpy(latvec_in,etot_in,parini%target_pressure_habohr,enthalpy)
write(*,'(a,i5,a,es15.7,a,es12.4)') " # Combined BFGS  exited in iterations: ", int(counter), " Enthalpy=",enthalpy," fmax=",fmax

end subroutine GEOPT_RBFGS_MHM
!contains

subroutine bfgs_driver_atoms(parini,parres,latvec_in,xred_in,fcart_in,strten_in,etot_in,iprec,counter,fmax_tol,folder)
 use global, only: units
 use defs_basis
!subroutine bfgsdriver(nat,nproc,iproc,rxyz,fxyz,epot,ncount_bigdft)!nproc,iproc,rxyz,fxyz,epot,at,rst,in,ncount_bigdft)
!    use module_base
!    use module_types
!    use module_interfaces
    use minpar
    use mod_parini, only: typ_parini
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_parini), intent(inout):: parres
    real(8) :: latvec_in(3,3),xred_in(3,parini%nat),fcart_in(3,parini%nat),etot_in,strten_in(6),enthalpy
    real(8), intent(inout) :: counter
!    type(atoms_data), intent(inout) :: at
!    type(input_variables), intent(inout) :: in
!    type(restart_objects), intent(inout) :: rst
    real(8) :: epot,tmplat(3,3),str_matrix(3,3),flat(3,3),pressure_mat(3,3),sigma(3,3),dstr(6),pressure,ent_pos_0,strtarget(6)
    real(8), dimension(3*parini%nat) :: rxyz
    real(8), dimension(3*parini%nat) :: fxyz
    real(8) :: fmax,fmax_at,fmax_lat,fmax_tol,en0000,betax
    integer :: infocode,i,ixyz,iat,istat,icall,icheck,istr,iexit
    character(len=4) :: fn4
    character(len=40) :: comment,filename,coord,folder
    logical ::  getwfk
    integer ::  nr
    integer ::  nwork,iprec
    real(8),allocatable:: x(:),f(:),work(:)
!Output BETAX and BETAX_LAT
write(*,'(a,es15.7,es15.7)') " # BFGS BETAX, BETAX_LAT: ", parmin_bfgs%betax, parmin_bfgs%betax_lat
coord="atoms"
!Reset counter
!counter=0.d0
pressure=parini%target_pressure_habohr
!Generate a set of variables containing all degrees of freedome
    call rxyz_int2cart(latvec_in,xred_in,rxyz,parini%nat)
!    counter=counter+1
!Here is the real driver of BFGS from reza
    icheck=0
    nr=3*parini%nat
    parmin_bfgs%iflag=0
    nwork=nr*nr+3*nr+3*nr*nr+3*nr
    allocate(work(nwork),stat=istat)
    if(istat/=0) stop 'ERROR: failure allocating work.'
    allocate(x(nr),stat=istat)
    if(istat/=0) stop 'ERROR: failure allocating x.'
    allocate(f(nr),stat=istat)
    if(istat/=0) stop 'ERROR: failure allocating f.'
    icall=0
    do 
!Here we perform the force call
       if(parini%usewf_geopt) then
           getwfk=.true.
       else
           getwfk=.false.
       endif
!       call get_BFGS_forces_PR(parini,rxyz,fxyz,epot,getwfk,iprec,latvec_in,xred_in,etot_in,fcart_in,strten_in)
!       call get_BFGS_forces_max(parini,rxyz,fxyz,epot,getwfk,iprec,latvec_in,xred_in,etot_in,fcart_in,strten_in)
       if(icall.ne.0) call get_BFGS_forces_atom(parini,parres,rxyz,fxyz,latvec_in,enthalpy,getwfk,iprec,latvec_in,&
                           &xred_in,etot_in,fcart_in,strten_in)
       if(icall.ne.0) counter=counter+1
       if(icall==0)   then
         do iat=1,parini%nat
           fxyz((iat-1)*3+1:iat*3)=fcart_in(:,iat)
         enddo
       endif
       call get_enthalpy(latvec_in,etot_in,pressure,enthalpy)
       epot=enthalpy
       if(icall==0) ent_pos_0=epot
!      call call_bigdft()!nproc,iproc,at,rxyz,in,epot,fxyz,fnoise,rst,infocode)

        !endif
        call atomic_copymoving_forward(parini%nat,3*parini%nat,fxyz,nr,f)
        call atomic_copymoving_forward(parini%nat,3*parini%nat,rxyz,nr,x)
!FIRE: check for convergence
!!Compute maximal component of forces, EXCLUDING any fixed components
 call get_fmax(parini,fcart_in,strten_in,fmax,fmax_at,fmax_lat)
         iexit=0
         if(fmax_at.lt.fmax_tol.or.fmax.lt.parini%paropt_geopt%fmaxtol) iexit=1
         if(iexit==1) parmin_bfgs%converged=.true.
          
!        call fnrmandforcemax(fxyz,fnrm,fmax,nat)
!        if(fmax<3.d-1) call updatefluctsum(nat,fnoise,fluct)
!        call convcheck()!fnrm,fmax,fluct*in%frac_fluct,in%forcemax,icheck)
!        if(iproc==0) write(*,*) 'ICHECK ',icheck
!        if(icheck>5) parmin_bfgs%converged=.true.
        !call calmaxforcecomponentanchors(atoms,np,f(1,1),fnrm,fspmax)
        !call checkconvergence(parmin_bfgs,fspmax)
        !if(ncount_bigdft>in%ncount_cluster_x-1)
        !if(iproc==0) write(*,*) 'nr=',nr,f(1)
!        if (iproc == 0) then
!           write(fn4,'(i4.4)') ncount_bigdft
!           write(comment,'(a,1pe10.3)')'BFGS:fnrm= ',sqrt(fnrm)
!           call  write_atomic_file()!trim(in%dir_output)//'posout_'//fn4,epot,rxyz,at,trim(comment),forces=fxyz)
!        endif

        call bfgs_reza(parini%nat,nr,x,epot,f,nwork,work,parmin_bfgs%betax,parmin_bfgs%betax_lat,fmax,fmax_at,fmax_lat,int(counter),coord)
!        call sd_minhocao(nat,nr,x,epot,f,parmin_bfgs%betax,parmin_bfgs%betax_lat,fmax,int(counter))



!        call bfgs_reza(iproc,nr,x,epot,f,nwork,work,parmin_bfgs%betax,sqrt(fnrm),fmax, &
!            ncount_bigdft,fluct*parmin_bfgs%frac_fluct,fluct)


        !x(1:nr)=x(1:nr)+1.d-2*f(1:nr)
        call atomic_copymoving_backward(parini%nat,nr,x,3*parini%nat,rxyz)
!        if(iexit==1) then !if(parmin_bfgs%converged) then
!        write(*,'(a,i0,a)') " #  BFGS converged in ",icall," iterations"
!MHM: Write output to file in every step***********************************
       call get_enthalpy(latvec_in,etot_in,pressure,enthalpy)
       en0000=enthalpy-ent_pos_0
       if(icall.ne.0.or.int(counter)==0) then
       write(fn4,'(i4.4)') int(counter)
       filename=trim(folder)//"posgeopt."//fn4//".ascii"
       units=units
       write(*,'(a,a)') " # Writing the positions in BFGS ATOMS  : ",filename
       call write_atomic_file_ascii(parini,filename,parini%nat,units,xred_in,latvec_in,fcart_in,strten_in,&
            &parini%char_type(1:parini%ntypat_global),parini%ntypat_global,parini%typat_global,parini%fixat,parini%fixlat,etot_in,pressure,enthalpy,en0000)
        write(*,'(a,i4,2(1x,es17.8))') " # GEOPT ",int(counter),enthalpy, fmax 
        if(iexit==1) then
          write(*,'(a,i4,2(1x,es25.15))') " #GEOPT converged", icall,enthalpy,fmax
          exit
        endif
       endif
!*********************************************************************
!              write(fn4,'(i4.4)') icall
!              write(comment,'(a,1pe10.3)')'BFGS:fnrm= ',sqrt(fnrm)
!              call  write_atomic_file()!trim(in%dir_output)//'posout_'//fn4,epot,rxyz,at,trim(comment),forces=fxyz)
!        endif
        !if(ncount_bigdft>in%ncount_cluster_x-1)
        !do ip=1,np-1
        !    call atomic_copymoving_backward(atoms,nr,xa(1,ip),n,x(1,ip))
        !enddo
!        if(parmin_bfgs%converged) exit
!        if(parmin_bfgs%iflag<=0) exit
        icall=icall+1
        if(int(counter)>parini%paropt_geopt%nit) exit
    enddo
    deallocate(work,stat=istat);if(istat/=0) stop 'ERROR: failure deallocating work.'
    deallocate(x,stat=istat);if(istat/=0) stop 'ERROR: failure deallocating x.'
    deallocate(f,stat=istat);if(istat/=0) stop 'ERROR: failure deallocating f.'
!contains
END SUBROUTINE

subroutine bfgs_driver_lattice(parini,parres,latvec_in,xred_in,fcart_in,strten_in,etot_in,iprec,counter,fmax_tol,folder)
 use global, only: units,reuse_kpt,ka1,kb1,kc1
 use defs_basis

!subroutine bfgsdriver(nat,nproc,iproc,rxyz,fxyz,epot,ncount_bigdft)!nproc,iproc,rxyz,fxyz,epot,at,rst,in,ncount_bigdft)
!    use module_base
!    use module_types
!    use module_interfaces
    use minpar
    use mod_parini, only: typ_parini
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_parini), intent(inout):: parres
    real(8) :: latvec_in(3,3),xred_in(3,parini%nat),fcart_in(3,parini%nat),etot_in,strten_in(6),enthalpy,latvec(9)
    real(8), intent(inout) :: counter
!    type(atoms_data), intent(inout) :: at
!    type(input_variables), intent(inout) :: in
!    type(restart_objects), intent(inout) :: rst
    real(8) :: epot,tmplat(3,3),str_matrix(3,3),flat(3,3),pressure_mat(3,3),sigma(3,3),dstr(6),pressure,ent_pos_0,strtarget(6)
    real(8) :: fmax,fmax_at,fmax_lat,fmax_tol,en0000
    integer :: infocode,i,ixyz,iat,istat,icall,icheck,istr,iexit
    character(len=4) :: fn4
    character(len=40) :: comment, filename,coord,folder
    logical ::  getwfk
    integer ::  nr
    integer ::  nwork,iprec
    real(8),allocatable:: x(:),f(:),work(:)
    real(8) :: transformed(3,3),transformed_inv(3,3),stressvol(3,3),vol
!We swich on the kpt history and maximization, and here we also initialize
    reuse_kpt=.true.

!Output BETAX and BETAX_LAT
write(*,'(a,es15.7,es15.7)') " # BFGS BETAX, BETAX_LAT: ", parmin_bfgs%betax, parmin_bfgs%betax_lat
write(*,'(a,i5)') " # MAX_LAT_ITER: ", parmin_bfgs%maxiter_lat
coord="lattice"
!Reset counter
!counter=0.d0
pressure=parini%target_pressure_habohr
latvec(1:3)=latvec_in(:,1)
latvec(4:6)=latvec_in(:,2)
latvec(7:9)=latvec_in(:,3)

!Here is the real driver of BFGS from reza
    icheck=0
    nr=9
    parmin_bfgs%iflag=0
    nwork=nr*nr+3*nr+3*nr*nr+3*nr
    allocate(work(nwork),stat=istat)
    if(istat/=0) stop 'ERROR: failure allocating work.'
    allocate(x(nr),stat=istat)
    if(istat/=0) stop 'ERROR: failure allocating x.'
    allocate(f(nr),stat=istat)
    if(istat/=0) stop 'ERROR: failure allocating f.'
    icall=0

       x=latvec
!Convert strten into "stress", which is actually the forces on the cell vectors
       str_matrix(1,1)=strten_in(1)
       str_matrix(2,2)=strten_in(2)
       str_matrix(3,3)=strten_in(3)
       str_matrix(1,2)=strten_in(6)
       str_matrix(2,1)=strten_in(6)
       str_matrix(1,3)=strten_in(5)
       str_matrix(3,1)=strten_in(5)
       str_matrix(2,3)=strten_in(4)
       str_matrix(3,2)=strten_in(4)
       call getvol(latvec_in,vol)
       transformed(1,:)=latvec_in(:,1)
       transformed(2,:)=latvec_in(:,2)
       transformed(3,:)=latvec_in(:,3)
       call invertmat(transformed,transformed_inv,3)
       flat=(-vol*matmul(str_matrix,transformed_inv))
       call stress_volume(latvec_in,vol,pressure,stressvol)
       flat=flat+stressvol
!Finally, write those values into F
        do iat=1,3
          f((iat-1)*3+1:iat*3)=flat(:,iat)
        enddo
    

    do 
!Here we perform the force call
       if(parini%usewf_geopt) then
           getwfk=.true.
       else
           getwfk=.false.
       endif
!       call get_BFGS_forces_PR(parini,rxyz,fxyz,epot,getwfk,iprec,latvec_in,xred_in,etot_in,fcart_in,strten_in)
!       if(icall.ne.0) call get_BFGS_forces_max(parini,xred_in,fxyz,enthalpy,getwfk,iprec,latvec_in,xred_in,etot_in,fcart_in,strten_in)
!       if(icall.ne.0) call get_BFGS_forces_atom(parini,rxyz,fxyz,latvec_in,enthalpy,getwfk,iprec,latvec_in,xred_in,etot_in,fcart_in,strten_in)
       if(ICALL.ne.0) call get_BFGS_forces_lattice(parini,parres,xred_in,f,x,enthalpy,getwfk,iprec,latvec_in,xred_in,etot_in,fcart_in,strten_in)
       if(icall.ne.0) counter=counter+1
       call get_enthalpy(latvec_in,etot_in,pressure,enthalpy)
       epot=enthalpy
       if(icall==0) ent_pos_0=epot
!      call call_bigdft()!nproc,iproc,at,rxyz,in,epot,fxyz,fnoise,rst,infocode)

        !endif
!FIRE: check for convergence
!!Compute maximal component of forces, EXCLUDING any fixed components
         call get_fmax(parini,fcart_in,strten_in,fmax,fmax_at,fmax_lat)
         iexit=0
         if(fmax_lat.lt.fmax_tol.or.fmax.lt.parini%paropt_geopt%fmaxtol) iexit=1
         if(iexit==1) parmin_bfgs%converged=.true.
          
!        call fnrmandforcemax(fxyz,fnrm,fmax,nat)
!        if(fmax<3.d-1) call updatefluctsum(nat,fnoise,fluct)
!        call convcheck()!fnrm,fmax,fluct*in%frac_fluct,in%forcemax,icheck)
!        if(iproc==0) write(*,*) 'ICHECK ',icheck
!        if(icheck>5) parmin_bfgs%converged=.true.
        !call calmaxforcecomponentanchors(atoms,np,f(1,1),fnrm,fspmax)
        !call checkconvergence(parmin_bfgs,fspmax)
        !if(ncount_bigdft>in%ncount_cluster_x-1)
        !if(iproc==0) write(*,*) 'nr=',nr,f(1)
!        if (iproc == 0) then
!           write(fn4,'(i4.4)') ncount_bigdft
!           write(comment,'(a,1pe10.3)')'BFGS:fnrm= ',sqrt(fnrm)
!           call  write_atomic_file()!trim(in%dir_output)//'posout_'//fn4,epot,rxyz,at,trim(comment),forces=fxyz)
!        endif

        call bfgs_reza(parini%nat,nr,x,epot,f,nwork,work,parmin_bfgs%betax,parmin_bfgs%betax_lat,fmax,fmax_at,fmax_lat,int(counter),coord)
!        call sd_minhocao(nat,nr,x,epot,f,parmin_bfgs%betax,parmin_bfgs%betax_lat,fmax,int(counter))



!        call bfgs_reza(iproc,nr,x,epot,f,nwork,work,parmin_bfgs%betax,sqrt(fnrm),fmax, &
!            ncount_bigdft,fluct*parmin_bfgs%frac_fluct,fluct)


        !x(1:nr)=x(1:nr)+1.d-2*f(1:nr)
!        if(iexit==1) then !if(parmin_bfgs%converged) then
!        write(*,'(a,i0,a)') " #  BFGS converged in ",icall," iterations"
!MHM: Write output to file in every step***********************************
       call get_enthalpy(latvec_in,etot_in,pressure,enthalpy)
       en0000=enthalpy-ent_pos_0
     if(icall.ne.0.or.int(counter)==0) then
       write(fn4,'(i4.4)') int(counter)
       filename=trim(folder)//"posgeopt."//fn4//".ascii"
       units=units
       write(*,'(a,a)') " # Writing the positions in BFGS2LATTICE: ",filename
       call write_atomic_file_ascii(parini,filename,parini%nat,units,xred_in,latvec_in,fcart_in,strten_in,&
            &parini%char_type(1:parini%ntypat_global),parini%ntypat_global,parini%typat_global,parini%fixat,parini%fixlat,etot_in,pressure,enthalpy,en0000)
        write(*,'(a,i4,2(1x,es17.8))') " # GEOPT ",int(counter),enthalpy, fmax 
        if(iexit==1) then
          write(*,'(a,i4,2(1x,es25.15))') " #GEOPT converged", icall,enthalpy,fmax
          exit
        endif
      endif
!*********************************************************************
!              write(fn4,'(i4.4)') icall
!              write(comment,'(a,1pe10.3)')'BFGS:fnrm= ',sqrt(fnrm)
!              call  write_atomic_file()!trim(in%dir_output)//'posout_'//fn4,epot,rxyz,at,trim(comment),forces=fxyz)
!        endif
        !if(ncount_bigdft>in%ncount_cluster_x-1)
        !do ip=1,np-1
        !    call atomic_copymoving_backward(atoms,nr,xa(1,ip),n,x(1,ip))
        !enddo
!        if(parmin_bfgs%converged) exit
!        if(parmin_bfgs%iflag<=0) exit
        icall=icall+1
        if(int(counter)>parini%paropt_geopt%nit.or.icall>parmin_bfgs%maxiter_lat) exit
    enddo
    deallocate(work,stat=istat);if(istat/=0) stop 'ERROR: failure deallocating work.'
    deallocate(x,stat=istat);if(istat/=0) stop 'ERROR: failure deallocating x.'
    deallocate(f,stat=istat);if(istat/=0) stop 'ERROR: failure deallocating f.'

!Switch back to default
    reuse_kpt=.false.
!contains
END SUBROUTINE
!subroutine inithess(iproc,nr,nat,rat,hess)!iproc,nr,nat,rat,atoms,hess)
!!    use module_types
!    implicit none
!    integer :: iproc,nr,nat,iat,jat,nsb,nrsqtwo,i,j,k,info
!    real(kind=8) :: rat(3,nat),hess(nr,nr),r0types(4,4),fctypes(4,4),soft,hard
!!    type(atoms_data), intent(inout) :: atoms
!    integer, allocatable::ita(:),isb(:,:)
!    real(8), allocatable::r0bonds(:),fcbonds(:),evec(:,:),eval(:),wa(:)
!    real(kind=8) :: dx,dy,dz,r,tt
!    nrsqtwo=2*nr**2
!    if(nr/=3*nat) then
!        stop 'ERROR: This subroutine works only for systems without fixed atoms.'
!    endif
!    allocate(ita(nat),isb(10*nat,2),r0bonds(10*nat),fcbonds(10*nat))
!    allocate(evec(nr,nr),eval(nr),wa(nrsqtwo))
!    do iat=1,nat
!!        if(trim(atoms%atomnames(atoms%iatype(iat)))=='H') then
!!            ita(iat)=1
!!        elseif(trim(atoms%atomnames(atoms%iatype(iat)))=='C') then
!!            ita(iat)=2
!!        elseif(trim(atoms%atomnames(atoms%iatype(iat)))=='N') then
!!            ita(iat)=3
!!        elseif(trim(atoms%atomnames(atoms%iatype(iat)))=='O') then
!!            ita(iat)=4
!!        else
!!            if(iproc==0) then
!!                write(*,'(a)') 'ERROR: This PBFGS is only implemented for systems which '
!!                write(*,'(a)') '       contain only organic elements, namely H,C,N,O.'
!!                write(*,'(a)') '       so use BFGS instead.'
!!            endif
!!            stop
!!        endif
!    enddo
!    call init_parameters(r0types,fctypes)
!    !r0types(1:4,1:4)=2.d0 ; fctypes(1:4,1:4)=5.d2
!    nsb=0
!    do iat=1,nat
!        do jat=iat+1,nat
!            dx=rat(1,jat)-rat(1,iat)
!            dy=rat(2,jat)-rat(2,iat)
!            dz=rat(3,jat)-rat(3,iat)
!            r=sqrt(dx**2+dy**2+dz**2)
!            !if(iat==21 .and. jat==27 .and. iproc==0) then
!            !    write(*,*) 'REZA ',r,1.35d0*r0types(ita(iat),ita(jat))
!            !endif
!            if(r<1.35d0*r0types(ita(iat),ita(jat))) then
!                nsb=nsb+1
!                if(nsb>10*nat) stop 'ERROR: too many stretching bonds, is everything OK?'
!                isb(nsb,1)=iat
!                isb(nsb,2)=jat
!                r0bonds(nsb)=r0types(ita(iat),ita(jat)) !CAUTION: equil. bond length from amber
!                !r0bonds(nsb)=r !CAUTION: current bond length assumed as equil. 
!                fcbonds(nsb)=fctypes(ita(iat),ita(jat))
!            endif
!        enddo
!    enddo
!    if(iproc==0) write(*,*) 'NSB ',nsb
!    !if(iproc==0) then
!    !    do i=1,nsb
!    !        write(*,'(a,i5,2f20.10,2i4,2(x,a))') 'PAR ', &
!    !            i,r0bonds(i),fcbonds(i),isb(i,1),isb(i,2), &
!    !            trim(atoms%atomnames(atoms%iatype(isb(i,1)))),trim(atoms%atomnames(atoms%iatype(isb(i,2))))
!    !    enddo
!    !endif
!    call pseudohess(nat,rat,nsb,isb(1,1),isb(1,2),fcbonds,r0bonds,hess)
!    evec(1:nr,1:nr)=hess(1:nr,1:nr)
!    !if(iproc==0) write(*,*) 'HESS ',hess(:,:)
!    call DSYEV('V','L',nr,evec,nr,eval,wa,nrsqtwo,info)
!    if(info/=0) stop 'ERROR: DSYEV in inithess failed.'
!    if(iproc==0) then
!        do i=1,nr
!            write(*,'(i5,es20.10)') i,eval(i)
!        enddo
!    endif
!    !stop
!    hard=eval(nr)
!    soft=eval(nr-nsb+1)
!    do k=1,nr
!        if(eval(k)<soft) then
!            eval(k)=soft
!        endif
!        eval(k)=1.d0/sqrt(eval(k)**2+soft**2)
!    enddo
!    do i=1,nr
!    do j=i,nr
!        tt=0.d0
!        do k=1,nr
!            !ep=1.d0/max(1.d-5,eval(k))
!            !ep=sqrt(ep**2+(20.d0/eval(nr))**2)
!            !if(eval(k
!            !ep=sqrt(eval(k)**2+constant**2)
!            tt=tt+eval(k)*evec(i,k)*evec(j,k)
!        enddo
!        hess(i,j)=tt
!    enddo
!    enddo
!    do i=1,nr
!    do j=1,i-1
!        hess(i,j)=hess(j,i)
!    enddo
!    enddo
!    deallocate(ita,isb,r0bonds,fcbonds,evec,eval,wa)
!end subroutine inithess


subroutine init_parameters(r0,fc)
    implicit none
    integer :: i,j
    real(kind=8) :: r0(4,4),fc(4,4)
    !((0.0104 / 0.239) / 27.2114) * (0.529177^2) = 0.000447802457
    r0(1,1)=0.80d0/0.529d0
    r0(2,1)=1.09d0/0.529d0 ; r0(2,2)=1.51d0/0.529d0
    r0(3,1)=1.01d0/0.529d0 ; r0(3,2)=1.39d0/0.529d0 ; r0(3,3)=1.10d0/0.529d0
    r0(4,1)=0.96d0/0.529d0 ; r0(4,2)=1.26d0/0.529d0 ; r0(4,3)=1.10d0/0.529d0 ; r0(4,4)=1.10/0.529d0
    do i=1,4
        do j=i+1,4
            r0(i,j)=r0(j,i)
        enddo
    enddo
    fc(1,1)=1.00d3*4.48d-4
    fc(2,1)=3.40d2*4.48d-4 ; fc(2,2)=3.31d2*4.48d-4
    fc(3,1)=4.34d2*4.48d-4 ; fc(3,2)=4.13d2*4.48d-4 ; fc(3,3)=4.56d3*4.48d-4
    fc(4,1)=5.53d2*4.48d-4 ; fc(4,2)=5.43d2*4.48d-4 ; fc(4,3)=4.56d3*4.48d-4 ; fc(4,4)=4.56d3*4.48d-4
    do i=1,4
        do j=i+1,4
            fc(i,j)=fc(j,i)
        enddo
    enddo
end subroutine init_parameters


subroutine pseudohess(nat,rat,nbond,indbond1,indbond2,sprcons,xl0,hess)
    implicit none
    integer :: nat,nbond,indbond1(nbond),indbond2(nbond)
    real(kind=8) :: rat(3,nat),sprcons(nbond),xl0(nbond),hess(3*nat,3*nat)
    integer :: iat,jat,i,j,ibond
    real(kind=8) :: dx,dy,dz,r2,r,rinv,r3inv !,rinv2,rinv4,rinv8,rinv10,rinv14,rinv16
    real(kind=8) :: dxsq,dysq,dzsq,dxdy,dxdz,dydz,tt1,tt2,tt3
    real(kind=8) :: h11,h22,h33,h12,h13,h23
    do j=1,3*nat
        do i=1,3*nat
            hess(i,j)=0.d0
        enddo
    enddo
    do ibond=1,nbond
        iat=indbond1(ibond)
        jat=indbond2(ibond)
        dx=rat(1,iat)-rat(1,jat)
        dy=rat(2,iat)-rat(2,jat)
        dz=rat(3,iat)-rat(3,jat)
        r2=dx**2+dy**2+dz**2
        r=sqrt(r2) ; rinv=1.d0/r ; r3inv=rinv**3
        !rinv2=1.d0/r2
        !rinv4=rinv2*rinv2
        !rinv8=rinv4*rinv4
        !rinv10=rinv8*rinv2
        !rinv14=rinv10*rinv4
        !rinv16=rinv8*rinv8
        dxsq=dx*dx ; dysq=dy*dy ; dzsq=dz*dz
        dxdy=dx*dy ; dxdz=dx*dz ; dydz=dy*dz
        !tt1=672.d0*rinv16
        !tt2=48.d0*rinv14
        !tt3=192.d0*rinv10
        !tt4=24.d0*rinv8
        tt1=sprcons(ibond)
        tt2=xl0(ibond)*rinv
        tt3=tt2*rinv**2
        !calculating the six distinct elements of 6 by 6 block
        !h11=dxsq*tt1-tt2-dxsq*tt3+tt4
        !h22=dysq*tt1-tt2-dysq*tt3+tt4
        !h33=dzsq*tt1-tt2-dzsq*tt3+tt4
        !h12=dxdy*tt1-dxdy*tt3
        !h13=dxdz*tt1-dxdz*tt3
        !h23=dydz*tt1-dydz*tt3

        !k_b*(1-l0/l+l0*(x_i-x_j)^2/l^3)
        h11=tt1*(1.d0-tt2+dxsq*tt3)
        h22=tt1*(1.d0-tt2+dysq*tt3)
        h33=tt1*(1.d0-tt2+dzsq*tt3)
        h12=tt1*dxdy*tt3
        h13=tt1*dxdz*tt3
        h23=tt1*dydz*tt3
        i=3*(iat-1)+1 ; j=3*(jat-1)+1
        !filling upper-left traingle (summing-up is necessary)
        hess(i+0,i+0)=hess(i+0,i+0)+h11
        hess(i+0,i+1)=hess(i+0,i+1)+h12
        hess(i+1,i+1)=hess(i+1,i+1)+h22
        hess(i+0,i+2)=hess(i+0,i+2)+h13
        hess(i+1,i+2)=hess(i+1,i+2)+h23
        hess(i+2,i+2)=hess(i+2,i+2)+h33
        !filling lower-right traingle (summing-up is necessary)
        hess(j+0,j+0)=hess(j+0,j+0)+h11
        hess(j+0,j+1)=hess(j+0,j+1)+h12
        hess(j+1,j+1)=hess(j+1,j+1)+h22
        hess(j+0,j+2)=hess(j+0,j+2)+h13
        hess(j+1,j+2)=hess(j+1,j+2)+h23
        hess(j+2,j+2)=hess(j+2,j+2)+h33
        !filling 3 by 3 block
        !summing-up is not needed but it may be necessary for PBC
        hess(i+0,j+0)=-h11 ; hess(i+1,j+0)=-h12 ; hess(i+2,j+0)=-h13
        hess(i+0,j+1)=-h12 ; hess(i+1,j+1)=-h22 ; hess(i+2,j+1)=-h23
        hess(i+0,j+2)=-h13 ; hess(i+1,j+2)=-h23 ; hess(i+2,j+2)=-h33
        !write(*,'(i3,5es20.10)') ibond,hess(i+0,i+0),tt1,tt2,tt3,xl0(ibond)
    enddo
    !filling the lower triangle of 3Nx3N Hessian matrix
    do i=1,3*nat-1
        do j =i+1,3*nat
            hess(j,i)=hess(i,j)
        enddo
    enddo
    
end subroutine pseudohess


subroutine bfgs_reza(nat,nr,x,epot,f,nwork,work,alphax_at,alphax_lat,fmax,fmax_at,fmax_lat,counter,coord)
    use minpar, only:parmin_bfgs
!    use module_types
    implicit none
    integer :: nat,nr,nwork,mf,my,ms,nrsqtwo,iw1,iw2,iw3,iw4,info,i,j,l,mx
    integer :: counter
    real(kind=8) :: x(nr),f(nr),epot,work(nwork),alphax_at,alphax_lat,alphax
!    type(atoms_data), intent(inout) :: atoms
    !real(8), allocatable::eval(:),umat(:)
    !type(parameterminimization)::parmin_bfgs
    real(kind=8) :: DDOT,tt1,tt2,de,fnrm,fmax,beta,beta_lat,fmax_at,fmax_lat
    real(kind=8),save :: tt1_sum
    real(kind=8) :: tt3,tt4,tt5,tt6
    real(8), save::epotold,alpha,alphamax,zeta
    logical, save::reset
    integer, save::isatur
    character(40):: coord

    if(trim(coord)=="lattice") then
      alphax=alphax_lat
    elseif(trim(coord)=="atoms") then
      alphax=alphax_at
    else
      stop "Wrong coordinates"
    endif

    if(nwork/=nr*nr+3*nr+3*nr*nr+3*nr) then
        stop 'ERROR: size of work array is insufficient.'
    endif
    nrsqtwo=nr*nr*2
    mf=nr*nr+1       !for force of previous iteration in wiki notation
    my=mf+nr         !for y_k in wiki notation
    ms=my+nr         !for s_k in wiki notation
    iw1=ms+nr        !work array to keep the hessian untouched
    iw2=iw1+nr*nr    !for work array of DSYTRF
    iw3=iw2+nrsqtwo  !for p_k in wiki notation
    mx =iw3+nr       !for position of previous iteration
    iw4=mx+nr        !for eigenvalues of inverse og hessian
    if(parmin_bfgs%iflag==0) then
        parmin_bfgs%iflag=1
        parmin_bfgs%converged=.false.   !! STEFAN Stefan stefan
        parmin_bfgs%iter=0
        epotold=epot
        alpha=8.d-1
        reset=.false.
        alphamax=0.9d0
        zeta=1.d0
        isatur=0
        tt1_sum=0.d0
!        open(unit=1390,file='bfgs_eigenvalues.dat',status='replace')
        open(unit=1390,file='bfgs_eigenvalues.dat',status='unknown',position='APPEND')
          write(1390,*) "New relaxation started"
        close(1390)
    else
        parmin_bfgs%iter=parmin_bfgs%iter+1
    endif
!    if(fnrm<min(6.d-2,max(1.d-2,2.d-3*sqrt(real(nr,8))))) then
!    if(fmax<min(1.5d-2,max(0.25d-2,0.5d-3*sqrt(real(nr,8))))) then
    if((trim(coord)=="atoms".and.fmax_at<2.d-3) .OR. (trim(coord)=="lattice".and.fmax_lat<5.d-3)) then
        if(isatur<99) isatur=isatur+1
    else
        isatur=0
    endif
    de=epot-epotold
    !fnrm=calnorm(nr,f);fmax=calmaxforcecomponent(nr,f)
!    if(iproc==0) then
    !write(*,'(a10,i5,es23.15,es11.3,2es12.5,2es12.4,i3)') &
    !    'GEOPT_BFGS',parmin_bfgs%iter,epot,de,fnrm,fmax,zeta,alpha,isatur
    !       '(I5,1x,I5,2x,a10,2x,1pe21.14,2x,e9.2,1(1pe11.3),3(1pe10.2),2x,a,I3,2x,a,1pe8.2E1)'
    if(parmin_bfgs%iter.ne.0.or.counter==0) then
if(trim(coord)=="atoms")    write(*,'(a,i5,2x,2x,1es21.14,2x,3(es15.7),es11.3,2x,a7,i3)')  " # GEOPT BFGS ATOMS  ",&
        counter,epot,fmax,fmax_at,fmax_lat,de,'isatur=',isatur
if(trim(coord)=="lattice")  write(*,'(a,i5,2x,2x,1es21.14,2x,3(es15.7),es11.3,2x,a7,i3)')  " # GEOPT BFGS2LATTICE",&
        counter,epot,fmax,fmax_at,fmax_lat,de,'isatur=',isatur
    endif
!    write(16,'(i5,1x,i5,2x,a10,2x,1es21.14,2x,es9.2,es11.3,2x,a7,i3)') &
!        ncount_bigdft,parmin_bfgs%iter,'GEOPT_BFGS',epot,de,fmax,'isatur=',isatur
!    endif
    close(16)
    open(unit=16,file='geopt.mon',status='unknown',position='APPEND')
    !if(parmin_bfgs%iter==602) then
    !    do i=1,nr/3
    !        write(31,*) x(i*3-2),x(i*3-1),x(i*3-0)
    !    enddo
    !    stop
    !endif
    !if(fmax<parmin_bfgs%fmaxtol) then
    if(parmin_bfgs%converged) then
        !parmin_bfgs%converged=.true.
        parmin_bfgs%iflag=0
!        if(iproc==0) then
        write(*,'(a,i4,es23.15,es12.5)') &
            ' # BFGS FINISHED: iter,epot,fmax ',parmin_bfgs%iter,epot,fmax
!        endif
        return
    endif

    !if(de>0.d0 .and. zeta>1.d-1) then
!    if(de>5.d-2) then
    if(de>8.d-3) then
        write(*,*) "# GEOPT RESET OCCURED, OPT 1"
        epot=epotold
        x(1:nr)=work(mx:mx-1+nr)
        f(1:nr)=work(mf:mf-1+nr)
        reset=.true.
        !alpha=max(alpha*0.5d0/1.1d0,1.d-2)
        zeta=max(zeta*2.d-1,1.d-3)
        isatur=0
    elseif(de>1.d-8) then
        write(*,*) "# GEOPT RESET OCCURED, OPT 2"
        reset=.true.
        zeta=max(zeta*2.d-1,1.d-3)
        isatur=0
    elseif(parmin_bfgs%iter.ne.0 .and. tt1_sum.lt.-1.d2) then
        write(*,*) "# GEOPT RESET OCCURED, OPT 3",tt1_sum
        tt1_sum=0.d0
        reset=.true.
        zeta=max(zeta*2.d-1,1.d-3)
        isatur=0
    else
        !zeta=1.d0
        !if(zeta>1.d-1) zeta=min(zeta*1.1d0,1.d0)
        zeta=min(zeta*1.1d0,1.d0)
        !isatur=isatur+1
    endif
    if(parmin_bfgs%iter==0 .or. reset) then
        reset=.false.
        !if(isatur>=10) then
        !    reset=.false.
        !    !alpha=5.d-1
        !endif

        if(trim(parmin_bfgs%approach)=='PBFGS') then
          !  call inithess(iproc,nr,atoms%nat,x,atoms,work(1))
        else
            work(1:nr*nr)=0.d0
            do i=1,nr
                work(i+(i-1)*nr)=zeta*alphax
!                if(i.le.3*nat) work(i+(i-1)*nr)=zeta*alphax
!                if(i.gt.3*nat) work(i+(i-1)*nr)=zeta*alphax_lat
!                if(i.le.3*nat) work(i+(i-1)*nr)=1.d0*alphax
!                if(i.gt.3*nat) work(i+(i-1)*nr)=1.d0*alphax_lat
            enddo
        endif
        work(iw3:iw3-1+nr)=zeta*alphax*f(1:nr)
!        work(iw3:iw3-1+nr-9)=zeta*alphax*f(1:nr-9)
!        work(iw3+nr-9:iw3-1+nr)=zeta*alphax_lat*f(nr-8:nr)
    else
        work(ms:ms-1+nr)=x(1:nr)-work(mx:mx-1+nr)
        work(my:my-1+nr)=work(mf:mf-1+nr)-f(1:nr)
        tt1=DDOT(nr,work(my),1,work(ms),1)
        do i=1,nr
            tt2=0.d0
            do j=1,nr
                tt2=tt2+work(i+(j-1)*nr)*work(my-1+j)
            enddo
            work(iw2-1+i)=tt2
        enddo
        tt2=DDOT(nr,work(my),1,work(iw2),1)
        !write(21,*) parmin_bfgs%iter,tt1,tt2
        !tt1=max(tt1,1.d-2)
        do i=1,nr
            do j=i,nr
                l=i+(j-1)*nr
                work(l)=work(l)+(tt1+tt2)*work(ms-1+i)*work(ms-1+j)/tt1**2- &
                    (work(iw2-1+i)*work(ms-1+j)+work(iw2-1+j)*work(ms-1+i))/tt1
                work(j+(i-1)*nr)=work(l)
            enddo
        enddo
        !do i=1,nr
        !    tt2=0.d0
        !    do j=1,nr
        !        tt2=tt2+work(j+(i-1)*nr)*f(j)
        !    enddo
        !    work(iw3-1+i)=tt2
        !enddo
        !write(31,*) zeta
        work(iw1:iw1-1+nr*nr)=work(1:nr*nr)
        call DSYEV('V','L',nr,work(iw1),nr,work(iw4),work(iw2),nrsqtwo,info)
        if(info/=0) stop 'ERROR: DSYEV in bfgs_reza failed.'
        tt1=work(iw4+0)    ; tt2=work(iw4+1)    ; tt3=work(iw4+2)
        tt4=work(iw4+nr-3) ; tt5=work(iw4+nr-2) ; tt6=work(iw4+nr-1)
!        if(iproc==0) then
        open(unit=1390,file='bfgs_eigenvalues.dat',status='old',position='append')
!        write(1390,'(i5,6es15.5)') parmin_bfgs%iter,tt1,tt2,tt3,tt4,tt5,tt6
        write(1390,'(i5,6es15.5)') counter,tt1,tt2,tt3,tt4,tt5,tt6
        close(1390)
        if(tt1.lt.0.d0) tt1_sum=tt1_sum+tt1
!        endif
        work(iw3:iw3-1+nr)=0.d0
        if(isatur<3) then
            beta=1.d-1/alphax
!            beta_lat=1.d-1/alphax_lat
        elseif(isatur<6) then
            beta=3.d-2/alphax
!            beta_lat=3.d-2/alphax_lat
        elseif(isatur<10) then
            beta=1.d-2/alphax
!            beta_lat=1.d-2/alphax_lat
        else
            beta=1.d-2/alphax
!            beta_lat=1.d-2/alphax_lat
        endif
!if(isatur<3) then
!beta=2.d-1/alphax
!beta_lat=2.d-1/alphax_lat
!elseif(isatur<6) then
!beta=5.d-2/alphax
!beta_lat=5.d-2/alphax_lat
!else
!beta=2.d-2/alphax
!beta_lat=2.d-2/alphax_lat
!endif
!        if(isatur<3) then
!            beta=1.d-1/alphax
!            beta_lat=1.d-1/alphax_lat
!        elseif(isatur<6) then
!            beta=4.d-2/alphax
!            beta_lat=4.d-2/alphax_lat
!        elseif(isatur<10) then
!            beta=4.d-2/alphax
!            beta_lat=4.d-2/alphax_lat
!        else
!            beta=4.d-2/alphax
!            beta_lat=4.d-2/alphax_lat
!        endif
        !do j=1,nr
        !    if(work(iw4-1+j)>0.d0) then
        !        tt3=work(iw4-1+j)
        !        exit
        !    enddo
        !enddo
        tt3=alphax*0.5d0
        do j=1,nr
            tt1=DDOT(nr,work(iw1+nr*(j-1)),1,f,1)
            if(work(iw4-1+j)<tt3) then
                tt4=tt3
            else
                tt4=work(iw4-1+j)
            endif
              tt2=1.d0/sqrt(1.d0/tt4**2+beta**2)
!            if(j.le.3*nat) then
!              tt2=1.d0/sqrt(1.d0/tt4**2+beta**2)
!            else
!              tt2=1.d0/sqrt(1.d0/tt4**2+beta_lat**2)
!            endif
            do i=1,nr
                work(iw3-1+i)=work(iw3-1+i)+tt1*work(iw1-1+i+nr*(j-1))*tt2
            enddo
        enddo
    endif
    epotold=epot
    work(mf:mf-1+nr)=f(1:nr)
    work(mx:mx-1+nr)=x(1:nr)
    alpha=min(alphamax,alpha*1.1d0)
    x(1:nr)=x(1:nr)+alpha*work(iw3:iw3-1+nr)
end subroutine bfgs_reza


!!> Driver for the LBFGS routine found on the Nocedal Homepage
!!! The subroutines have only been modified slightly, so a VIMDIFF will show all modifications!
!!! This is helpfull when we are looking for the source of problems during BFGS runs
subroutine lbfgs_driver_lattice(parini,parres,latvec_in,xred_in,fcart_in,strten_in,etot_in,iprec,counter,fail,fmax_tol,folder)
!This routine expects to receive "good" forces and energies initially
 use global, only: units,reuse_kpt,ka1,kb1,kc1
 use defs_basis

!subroutine lbfgsdriver(rxyz,fxyz,etot,at,rst,in,ncount_bigdft,fail) 
 use minpar
 use mod_parini, only: typ_parini
 implicit none
  type(typ_parini), intent(in):: parini
  type(typ_parini), intent(inout):: parres
  real(8) :: latvec_in(3*3),xred_in(3,parini%nat),fcart_in(3,parini%nat),etot_in,strten_in(6),enthalpy,enthalpyprev
  real(8), intent(inout) :: counter
!  real(8), dimension(3*at%nat), intent(inout) :: rxyz
  logical, intent(out) :: fail
!  real(8), dimension(3*at%nat), intent(out) :: fxyz
!  real(8), dimension(3*at%nat):: txyz, sxyz
!  real(8) :: fluct,fnrm, fnoise
  real(8) :: fmax,fmax_lat,fmax_at,fmax_tol,latvec_write(3*3),pressure,strtarget(6),dstr(6),de,str_matrix(3,3),vol
  real(8) :: flat(3,3),transformed(3,3),transformed_inv(3,3),stressvol(3,3)
  real(8) :: strten_write(6),fcart_write(3,parini%nat),etot_write
!  logical :: check
  integer :: check,istr,iexit,iprec
  integer :: ixyz,iat,i,iproc
  character(len=4) :: fn4
  character(len=40) :: comment
  logical :: getwfk

  integer ::  n,nr,ndim
  integer ::  NWORK
  real(8),allocatable:: X(:),G(:),DIAG(:),W(:)
  real(8):: F,EPS!,XTOL,GTOL,,STPMIN,STPMAX
!  real(8), dimension(3*at%nat) :: rxyz0,rxyzwrite
  integer ::  IPRINT(2),IFLAG,ICALL,M
!  character(len=*), parameter :: subname='bfgs'
  integer :: i_stat,i_all
  character(40):: filename,folder
write(*,'(a,i5)') " # MAX_LAT_ITER: ", parmin_bfgs%maxiter_lat
!Initialize parameters
  latvec_write=latvec_in
  call geopt_init()
  pressure=parini%target_pressure_habohr
  check=0
  iproc=0
  latvec_write=latvec_in
  etot_write=etot_in
  fcart_write=fcart_in
  strten_write=strten_in

!  call init_driver(par)     !Initialize the parameters
  parmin_bfgs%finstep=0
  parmin_bfgs%alpha=1.d0
  fail=.false.


  
  enthalpyprev=0.d0
  !check if the convergence is reached after entering routine
  call get_fmax(parini,fcart_in,strten_in,fmax,fmax_at,fmax_lat)
         iexit=0
         if(fmax_lat.lt.fmax_tol.or.fmax.lt.parini%paropt_geopt%fmaxtol) iexit=1
         if(iexit==1) parmin_bfgs%converged=.true.
         if(iexit==1) then
               write(*,*) "# Lattice L-BFGS converged before entering routine"
               return
         endif 

!Make a list of all degrees of freedom that should be passed to bfgs
  nr=9
     NDIM=nr
     NWORK=NDIM*(2*parmin_bfgs%MSAVE +1)+2*parmin_bfgs%MSAVE
      
     allocate(X(NDIM),stat=i_stat)
     allocate(G(NDIM),stat=i_stat)
     allocate(DIAG(NDIM),stat=i_stat)
     allocate(W(NWORK),stat=i_stat)

!We swich on the kpt history and maximization, and here we also initialize
    reuse_kpt=.true.

     X=latvec_in
!Convert strten into "stress", which is actually the forces on the cell vectors
       str_matrix(1,1)=strten_in(1)
       str_matrix(2,2)=strten_in(2)
       str_matrix(3,3)=strten_in(3)
       str_matrix(1,2)=strten_in(6)
       str_matrix(2,1)=strten_in(6)
       str_matrix(1,3)=strten_in(5)
       str_matrix(3,1)=strten_in(5)
       str_matrix(2,3)=strten_in(4)
       str_matrix(3,2)=strten_in(4)
       call getvol(latvec_in,vol)
       transformed(1,:)=latvec_in(1:3)
       transformed(2,:)=latvec_in(4:6)
       transformed(3,:)=latvec_in(7:9)


       call invertmat(transformed,transformed_inv,3)
       flat=(-vol*matmul(str_matrix,transformed_inv))
       call stress_volume(latvec_in,vol,pressure,stressvol)
       flat=flat+stressvol
!Finally, write those values into G
        do iat=1,3
          G((iat-1)*3+1:iat*3)=flat(:,iat)
        enddo

     N=nr
     M=parmin_bfgs%MSAVE
     IPRINT(1)= 1
     IPRINT(2)= 0
     call get_enthalpy(latvec_in,etot_in,pressure,enthalpy)
     F=enthalpy
!     We do not wish to provide the diagonal matrices Hk0, and 
!     therefore set DIAGCO to FALSE.a
!     Lets try true for once
     parmin_bfgs%DIAGCO=.true.;DIAG=parmin_bfgs%BETAX_LAT
     EPS=0.0d0
     ICALL=0
     IFLAG=0

 20   CONTINUE
!        if (parmin_bfgs%IWRITE) then
         if(icall.ne.0) then
!MHM: Write output to file in every step***********************************
           call get_enthalpy(latvec_in,etot_in,pressure,enthalpy)
           de=enthalpy-enthalpyprev
         if(icall.ge.2.or.counter==0) then
           write(fn4,'(i4.4)') int(counter)
           filename=trim(folder)//"posgeopt."//fn4//".ascii"
           units=units
           write(*,'(a,a)') " # Writing the positions in BFGS LATTICE: ",filename
           call write_atomic_file_ascii(parini,filename,parini%nat,units,xred_in,latvec_in,fcart_in,strten_in,&
                &parini%char_type(1:parini%ntypat_global),parini%ntypat_global,parini%typat_global,parini%fixat,parini%fixlat,etot_in,pressure,enthalpy,real(icall,8))
           write(*,'(a,i4,4(1x,es17.8),1x,i5)') " # GEOPT ",int(counter),enthalpy,fmax,fmax_at,fmax_lat,icall
           write(*,'(a,i5,2x,2x,1es21.14,2x,3(es15.7),es11.3,2x,a8,i3,1x,a6,1x,1pe8.2E1)')  " # GEOPT BFGS LATTICE",&
           int(counter),enthalpy,fmax,fmax_at,fmax_lat,de,"BFGS-it=",parmin_bfgs%finstep,"alpha=",parmin_bfgs%alpha
            if(iexit==1) then
              write(*,'(a,i4,2(1x,es25.15),1x,i5)') " # GEOPT LBFGS converged", int(counter),enthalpy,fmax_lat,icall
              goto 100
            endif
!*********************************************************************
           if(parmin_bfgs%iwrite) then
             parmin_bfgs%IWRITE=.false.
             latvec_write=latvec_in
             etot_write=etot_in
             fcart_write=fcart_in
             strten_write=strten_in
            endif
        endif
           enthalpyprev=enthalpy
        endif
!              write(*,'(I5,1x,I5,2x,a11,1x,1pe21.14,2x,e9.2,1(1pe11.3),3(1pe10.2),2x,a,I3,2x,a,1pe8.2E1)')&
!              &ncount_bigdft,ICALL,"GEOPT_LBFGS",enthalpy,enthalpy-enthalpyprev,fmax,sqrt(fnrm),fluct*in%frac_fluct,fluct&
!              &,"BFGS-it=",parmin_bfgs%finstep,"alpha=",parmin_bfgs%alpha
!   if (iproc==0.and.ICALL.ne.0.and.parmin_bfgs%verbosity > 0) & 
!              & write(* ,'(I5,1x,I5,2x,a11,1x,1pe21.14,2x,e9.2,1(1pe11.3),3(1pe10.2),2x,a,I3,2x,a,1pe8.2E1)')&
!              &ncount_bigdft,ICALL,"GEOPT_LBFGS",enthalpy,enthalpy-enthalpyprev,fmax,sqrt(fnrm),fluct*in%frac_fluct,fluct&
!              &,"BFGS-it=",parmin_bfgs%finstep,"alpha=",parmin_bfgs%alpha
!             if (iproc==0.and.ICALL.ne.0.and.parmin_bfgs%verbosity > 0) write(*,'(1x,a,1pe14.5,2(1x,a,1pe14.5))')&
!                          'FORCES norm(Ha/Bohr): maxval=',fmax,'fnrm2=',fnrm,'fluct=', fluct
!              call convcheck(fnrm,fmax,fluct*in%frac_fluct, in%forcemax,check)
!              if (ncount_bigdft >= in%ncount_cluster_x) goto 50
!         if(iexit==1) then
!         write(*,'(a,i0,a)') " # LATTICE BFGS converged in ",ICALL," iterations"
!         if (iproc == 0) then
!            write(fn4,'(i4.4)') ncount_bigdft
!            write(comment,'(a,1pe10.3)')'BFGS:fnrm= ',sqrt(fnrm)
!            call  write_atomic_file(trim(in%dir_output)//'posout_'//fn4,etot,rxyz,at,trim(comment),forces=fxyz)
!         endif
!         goto 100
!         endif

      
!      rxyz=rxyz0
!      call atomic_copymoving_backward(at,nr,X,n,rxyz)
!      in%inputPsiId=1
!!      if(ICALL.ne.0) call call_bigdft(nproc,iproc,at,rxyz,in,F,fxyz,rst,infocode)
!      if(ICALL.ne.0) call call_bigdft()!nproc,iproc,at,rxyz,in,F,fxyz,fnoise,rst,infocode)
!Here we perform the force call
       if(parini%usewf_geopt) then
           getwfk=.true.
       else
           getwfk=.false.
       endif
      if(ICALL.ne.0) call get_BFGS_forces_lattice(parini,parres,xred_in,G,X,enthalpy,getwfk,iprec,latvec_in,xred_in,etot_in,fcart_in,strten_in)
      if(ICALL.ne.0) counter=counter+1
      F=enthalpy
      G=-G
!      call fnrmandforcemax(fxyz,fnrm,fmax,at%nat)
!      call fnrmandforcemax(fxyz,fnrm,fmax,at)
!check if the convergence is reached
  call get_fmax(parini,fcart_in,strten_in,fmax,fmax_at,fmax_lat)
         iexit=0
         if(fmax_lat.lt.fmax_tol.or.fmax.lt.parini%paropt_geopt%fmaxtol) iexit=1
         if(iexit==1) parmin_bfgs%converged=.true.
         if(int(counter)>parini%paropt_geopt%nit) goto 50
         if(icall>parmin_bfgs%maxiter_lat) goto 50

      CALL LBFGS(IPROC,parmin_bfgs,N,M,X,F,G,DIAG,IPRINT,EPS,W,IFLAG)
!Switch of the DIAGCO after the first call
      parmin_bfgs%diagco=.false.      
      IF(IFLAG.LE.0) GO TO 50
      ICALL=ICALL + 1
!     We allow at most the given number of evaluations of F and G
!      if(ncount_bigdft>in%ncount_cluster_x-1)  then
      if(int(counter).gt.parini%paropt_geopt%nit) then 
        goto 100
      endif
      close(16)
      open(unit=16,file='geopt.mon',status='unknown',position='append')
      GO TO 20
  50  CONTINUE
        write(*,*) "# Error or limit reached in GEOPT LATTICE, exiting"
        latvec_in=latvec_write
        etot_in=etot_write
        fcart_in=fcart_write
        strten_in=strten_write
        fail=.true.
!MHM: Write output to file in every step***********************************
           call get_enthalpy(latvec_in,etot_in,pressure,enthalpy)
           de=enthalpy-enthalpyprev
           write(fn4,'(i4.4)') int(counter)
           filename=trim(folder)//"posgeopt."//fn4//".ascii"
           units=units
           write(*,'(a,a)') " # Writing the positions in BFGS LATTICE: ",filename
           call write_atomic_file_ascii(parini,filename,parini%nat,units,xred_in,latvec_in,fcart_in,strten_in,&
                &parini%char_type(1:parini%ntypat_global),parini%ntypat_global,parini%typat_global,parini%fixat,parini%fixlat,etot_in,pressure,enthalpy,real(icall,8))
           write(*,'(a,i4,4(1x,es17.8),1x,i5)') " # GEOPT ",int(counter),enthalpy,fmax,fmax_at,fmax_lat,icall
           write(*,'(a,i5,2x,2x,1es21.14,2x,3(es15.7),es11.3,2x,a8,i3,1x,a6,1x,1pe8.2E1)')  " # GEOPT BFGS LATTICE",&
           int(counter),enthalpy,fmax,fmax_at,fmax_lat,de,"BFGS-it=",parmin_bfgs%finstep,"alpha=",parmin_bfgs%alpha
!*********************************************************************
 100  CONTINUE
!Switch back to default
    reuse_kpt=.false.
      deallocate(X,stat=i_stat)
      deallocate(G,stat=i_stat)
      deallocate(DIAG,stat=i_stat)
      deallocate(W,stat=i_stat)
END SUBROUTINE 

subroutine atomic_copymoving_forward(nat,n,x,nr,xa)
!    use module_types
    implicit none
    integer :: n,nr,i,iat,ixyz,ir,nat
    real(kind=8) :: x(n),xa(nr)
    logical :: move_this_coordinate
    ir=0
    do i=1,3*nat
        iat=(i-1)/3+1
        ixyz=mod(i-1,3)+1
!        if(move_this_coordinate(atoms%ifrztyp(iat),ixyz)) then
            ir=ir+1
            xa(ir)=x(i)
!        endif
    enddo
    if(ir/=nr) stop 'ERROR: inconsistent number of relaxing DOF'
END SUBROUTINE atomic_copymoving_forward


subroutine atomic_copymoving_backward(nat,nr,xa,n,x)
!    use module_types
    implicit none
    integer :: n,nr,i,iat,ixyz,ir,nat
    real(kind=8) :: x(n),xa(nr)
    logical :: move_this_coordinate
    ir=0
    do i=1,3*nat
        iat=(i-1)/3+1
        ixyz=mod(i-1,3)+1
!        if(move_this_coordinate(atoms%ifrztyp(iat),ixyz)) then
            ir=ir+1
            x(i)=xa(ir)
!        endif
    enddo
    if(ir/=nr) stop 'ERROR: inconsistent number of relaxing DOF'
END SUBROUTINE atomic_copymoving_backward



subroutine get_BFGS_forces_PR(parini,parres,pos_all,force_all,enthalpy,getwfk,iprec,latvec_in,xred_in,etot_in,fcart_in,strten_in)
!This routine hides away all cumbersome conversion of arrays in lattice and positions and forces and stresses
!such that they can be directly passed on to bfgs. It also outputs the enthalpy instead of the energy
use interface_code
use mod_parini, only: typ_parini
implicit none
type(typ_parini), intent(in):: parini
type(typ_parini), intent(inout):: parres
integer:: iprec,iat
real(8):: pos_all(3*parini%nat+9)
real(8):: force_all(3*parini%nat+9)
real(8):: enthalpy,pressure
real(8):: xred_in(3,parini%nat),fcart_in(3,parini%nat),strten_in(6),etot_in,latvec_in(3,3)
real(8):: str_matrix(3,3),flat(3,3),pressure_mat(3,3),tmplat(3,3),sigma(3,3),crossp(3)
logical:: getwfk

!Here we perform the force call
!****************************************************************************************************************        
!****************************************************************************************************************        
!Get the current acell_in and rprim_in, and also the atomic positions
!Generate a set of variables containing all degrees of freedome
    do iat=1,parini%nat
      xred_in(:,iat)=pos_all((iat-1)*3+1:iat*3)
    enddo
    do iat=1,3
      latvec_in(:,iat)=pos_all(3*parini%nat+(iat-1)*3+1:3*parini%nat+iat*3)
    enddo
       call get_energyandforces_single(parini,parres,latvec_in,xred_in,fcart_in,strten_in,etot_in,iprec,getwfk)
!****************************************************************************************************************   
!Conversion of forces is more complicate: 
!start with atomic forces
        call invertmat(latvec_in,tmplat,3)
        do iat=1,parini%nat
          force_all((iat-1)*3+1:iat*3)=matmul(tmplat,fcart_in(:,iat))
        enddo
!now the stresses   
!Setup pressure matrix
pressure=parini%target_pressure_habohr
pressure_mat=0.d0
pressure_mat(1,1)=1.d0;pressure_mat(2,2)=1.d0;pressure_mat(3,3)=1.d0
pressure_mat=pressure_mat*pressure  !Here the pressure is not passed to the energyandforces, so we move on the ENERGY surface
!Transform the forces on the lattice vectors into the stress tensore according to the paper Tomas Bucko and Jurg Hafner
!This procedure is no longer necessary, since we have the correct stress tensor now!
           str_matrix(1,1)=strten_in(1)
           str_matrix(2,2)=strten_in(2)
           str_matrix(3,3)=strten_in(3)
           str_matrix(1,2)=strten_in(6)
           str_matrix(2,1)=strten_in(6)
           str_matrix(1,3)=strten_in(5)
           str_matrix(3,1)=strten_in(5)
           str_matrix(2,3)=strten_in(4)
           str_matrix(3,2)=strten_in(4)
           flat=-str_matrix
!Here the pressure is applied
           flat=flat-pressure_mat
!We need to multiply with sigma (or eta in my paper)
           call cross_product(latvec_in(:,2),latvec_in(:,3),crossp); sigma(:,1)=crossp
           call cross_product(latvec_in(:,3),latvec_in(:,1),crossp); sigma(:,2)=crossp
           call cross_product(latvec_in(:,1),latvec_in(:,2),crossp); sigma(:,3)=crossp
!Compute the F
           flat=matmul(flat(:,:),sigma(:,:))
!Finally, write those values into fxyz
        do iat=1,3
          force_all(3*parini%nat+(iat-1)*3+1:3*parini%nat+iat*3)=flat(:,iat)
        enddo
!Get the enthalpy
        call get_enthalpy(latvec_in,etot_in,pressure,enthalpy)
end subroutine

subroutine getvol_strain(strain,latvec0,vol)
!Compute the volume based on strain and relative cell
implicit none
real(8), dimension(3,3):: latvec0,latvec,strain,unitmat
real(8):: vol
unitmat=0.d0
unitmat(1,1)=1.d0
unitmat(2,2)=1.d0
unitmat(3,3)=1.d0
    latvec=matmul(unitmat+strain,latvec0) 
    call getvol(latvec,vol)
end subroutine

subroutine  get_BFGS_forces_strainlatt(parini,parres,pos_all,force_all,enthalpy,getwfk,iprec,latvec0,&
           &lattdeg,latvec_in,xred_in,etot_in,fcart_in,strten_in)
!This routine hides away all cumbersome conversion of arrays in lattice and positions and forces and stresses
!such that they can be directly passed on to bfgs. It also outputs the enthalpy instead of the energy
use interface_code

use mod_parini, only: typ_parini
implicit none
type(typ_parini), intent(in):: parini
type(typ_parini), intent(inout):: parres
integer:: iprec,iat,lattdeg
real(8):: pos_all(3*parini%nat+9)
real(8):: force_all(3*parini%nat+9)
real(8):: enthalpy,pressure,vol,unitmat(3,3)
real(8):: xred_in(3,parini%nat),fcart_in(3,parini%nat),strten_in(6),etot_in,latvec_in(3,3),transformed(3,3),transformed_inv(3,3)
real(8):: str_matrix(3,3),flat(3,3),pressure_mat(3,3),tmplat(3,3),sigma(3,3),crossp(3),stressvol(3,3),latvec0(3,3)
real(8):: g(3,3),lattrans(3,3),strainall(3,3),strainalltrans(3,3)
logical:: getwfk

unitmat=0.d0
unitmat(1,1)=1.d0
unitmat(2,2)=1.d0
unitmat(3,3)=1.d0
transformed(:,1)=pos_all(3*parini%nat+(1-1)*3+1:3*parini%nat+1*3) !!latvec_in(1,:)
transformed(:,2)=pos_all(3*parini%nat+(2-1)*3+1:3*parini%nat+2*3) !!latvec_in(2,:)
transformed(:,3)=pos_all(3*parini%nat+(3-1)*3+1:3*parini%nat+3*3) !!latvec_in(3,:)
transformed=transpose(transformed)
if(enthalpy==-1.d10) goto 1000
!Here we perform the force call
!****************************************************************************************************************        
!****************************************************************************************************************        
!Get the current acell_in and rprim_in, and also the atomic positions
!Generate a set of variables containing all degrees of freedome
    do iat=1,parini%nat
      xred_in(:,iat)=pos_all((iat-1)*3+1:iat*3)
    enddo
if(lattdeg==1) then
    do iat=1,3
      latvec_in(:,iat)=pos_all(3*parini%nat+(iat-1)*3+1:3*parini%nat+iat*3)
    enddo
elseif(lattdeg==2) then 
    do iat=1,3
      strainall(:,iat)=pos_all(3*parini%nat+(iat-1)*3+1:3*parini%nat+iat*3)
    enddo
    latvec_in=matmul(unitmat+strainall,latvec0) 
endif
       call get_energyandforces_single(parini,parres,latvec_in,xred_in,fcart_in,strten_in,etot_in,iprec,getwfk)
!****************************************************************************************************************   
1000 continue
!Conversion of forces is more complicate: 
!start with atomic forces
!     if(lattdeg==1) then
       call fxyz_cart2int(parini%nat,fcart_in,force_all(1:3*parini%nat),latvec_in)
!     elseif(lattdeg==2) then
!       do iat=1,nat
!           force_all((iat-1)*3+1:iat*3)=matmul(matmul(transformed,latvec_in),fcart_in(:,iat))
!       enddo
!     endif
!now the stresses   
!Convert strten into "stress", which is actually the forces on the cell vectors
       str_matrix(1,1)=strten_in(1)
       str_matrix(2,2)=strten_in(2)
       str_matrix(3,3)=strten_in(3)
       str_matrix(1,2)=strten_in(6)
       str_matrix(2,1)=strten_in(6)
       str_matrix(1,3)=strten_in(5)
       str_matrix(3,1)=strten_in(5)
       str_matrix(2,3)=strten_in(4)
       str_matrix(3,2)=strten_in(4)
       call getvol(latvec_in,vol)
       pressure=parini%target_pressure_habohr
if(lattdeg==1) then
       call invertmat(transpose(latvec_in),transformed_inv,3)
       flat=(-vol*matmul(str_matrix,transformed_inv))
       call stress_volume(latvec_in,vol,pressure,stressvol)
       flat=flat+stressvol
elseif(lattdeg==2) then
       flat=(str_matrix+unitmat*pressure)*vol
       strainalltrans(1,:)=strainall(:,1)
       strainalltrans(2,:)=strainall(:,2)
       strainalltrans(3,:)=strainall(:,3)
       transformed=unitmat+strainalltrans
       call invertmat(transformed,transformed_inv,3)
       flat=-matmul(flat,transformed_inv)
endif
!Finally, write those values into fxyz
        do iat=1,3
          force_all(3*parini%nat+(iat-1)*3+1:3*parini%nat+iat*3)=flat(:,iat)
        enddo
!Get the enthalpy
        call get_enthalpy(latvec_in,etot_in,pressure,enthalpy)
end subroutine

subroutine correct_hessin(hess,hessin,latvec,ndim,hessupdate,lattdeg)
use minpar
!Depending on hessupdate, hessin will be updated or additionally hess
implicit none
integer:: ndim,LWORK,info,i,j,hessupdate,lattdeg
real(8):: hessin(ndim,ndim),hess(ndim,ndim),hess_tmp(ndim,ndim),dmat(ndim,ndim),latvec(3,3)
real(8), allocatable:: WORK(:),eval(:)
real(8):: minev
minev=0.01d0 !1.d-6
hess_tmp=hessin
LWORK=-1
allocate(WORK(1),eval(ndim))
        call DSYEV('V','L',ndim,hess_tmp,ndim,eval,WORK,LWORK,INFO)
        if (info.ne.0) stop 'DSYEV in correct_hessin'
LWORK=WORK(1)
deallocate(WORK)
allocate(WORK(LWORK))
        call DSYEV('V','L',ndim,hess_tmp,ndim,eval,WORK,LWORK,INFO)
        if (info.ne.0) stop 'DSYEV in correct_hessin'
!!!!        write(*,*) '---   App. eigenvalues in a.u. -------------'
!!!!        do j=1,3*nat+9
!!!!        if(j.ne.3*nat+9) write(*,'(1x,es8.1)',advance="no") eval(j)
!!!!        if(j==3*nat+9) write(*,'(1x,es8.1)') eval(j)
!!!!        enddo
!If the eigenvalues of the inverse hessian becomes negative, we need to do something...
!We have to shift up the negative eigenvalues to a small positive value, and slightly decrease the
!EV of the hard modes as well
dmat=0.d0
if(eval(1).lt.0.d0) then
write(*,*) "# Negative eigenvalue encountered, correcting the hessian..."
  do i=1,ndim
!    if(eval(i).lt.0.d0) then
      eval(i)=1.d0/sqrt(1.d0/eval(i)**2+minev**2)
    dmat(i,i)=eval(i)
  enddo
  hessin=matmul(matmul(hess_tmp,dmat),transpose(hess_tmp))
!Eliminate all off-diagonals
  hess_tmp=hessin
  hessin=0.d0
  do i=1,ndim
    hessin(i,i)=hess_tmp(i,i)
  enddo
  if(hessupdate==2)  then
  call invertmat(hessin,hess,ndim)
!     do i=1,ndim
!     dmat(i,i)=1.d0/dmat(i,i)
!     enddo   
!  hess=matmul(matmul(hess_tmp,dmat),transpose(hess_tmp))  
!!Eliminate all off-diagonals
!  hess_tmp=hess
!  hess=0.d0
!  do i=1,ndim
!    hess(i,i)=hess_tmp(i,i)
!  enddo
  endif
!!Instead, setup completely new hessian
!  call init_hessinv(hessin,latvec,parmin_bfgs%betax,parmin_bfgs%betax_lat,lattdeg) 
!  if(hessupdate==2) call invertmat(hessin,hess,ndim)
endif

!write(*,*) "ALPHA_1",hessin
!write(*,*) "ALPHA_2",hessin
end subroutine

!************************************************************************************

SUBROUTINE unit_matrix(mat,ndim)
implicit none
real(8),DIMENSION(ndim,ndim), INTENT(INOUT) :: mat
integer:: ndim
INTEGER :: i,n
!Action:
!Sets the diagonal components of mat to unity, all other components to zero.
!When mat is square, this will be the unit matrix; otherwise, a unit matrix
!with appended rows or columns of zeros.
mat(:,:)=0.0d0
n=min(size(mat,1),size(mat,2))
do i=1,n
mat(i,i)=1.0d0
end do
END SUBROUTINE unit_matrix

!************************************************************************************

subroutine get_BFGS_forces_max(parini,parres,pos_all,force_all,enthalpy,getwfk,iprec,latvec_in,xred_in,etot_in,fcart_in,strten_in)
!This routine hides away all cumbersome conversion of arrays in lattice and positions and forces and stresses
!such that they can be directly passed on to bfgs. It also outputs the enthalpy instead of the energy
use interface_code

use mod_parini, only: typ_parini
implicit none
type(typ_parini), intent(in):: parini
type(typ_parini), intent(inout):: parres
integer:: iprec,iat
real(8):: pos_all(3*parini%nat+9)
real(8):: force_all(3*parini%nat+9)
real(8):: enthalpy,pressure,vol
real(8):: xred_in(3,parini%nat),fcart_in(3,parini%nat),strten_in(6),etot_in,latvec_in(3,3),transformed(3,3),transformed_inv(3,3)
real(8):: str_matrix(3,3),flat(3,3),pressure_mat(3,3),tmplat(3,3),sigma(3,3),crossp(3),stressvol(3,3)
logical:: getwfk

!Here we perform the force call
!****************************************************************************************************************        
!****************************************************************************************************************        
!Get the current acell_in and rprim_in, and also the atomic positions
!Generate a set of variables containing all degrees of freedome
    do iat=1,parini%nat
      xred_in(:,iat)=pos_all((iat-1)*3+1:iat*3)
    enddo
    do iat=1,3
      latvec_in(:,iat)=pos_all(3*parini%nat+(iat-1)*3+1:3*parini%nat+iat*3)
    enddo
       call get_energyandforces_single(parini,parres,latvec_in,xred_in,fcart_in,strten_in,etot_in,iprec,getwfk)
!****************************************************************************************************************   
!Conversion of forces is more complicate: 
!start with atomic forces
        call fxyz_cart2int(parini%nat,fcart_in,force_all(1:3*parini%nat),latvec_in)
!now the stresses   
!Convert strten into "stress", which is actually the forces on the cell vectors
       str_matrix(1,1)=strten_in(1)
       str_matrix(2,2)=strten_in(2)
       str_matrix(3,3)=strten_in(3)
       str_matrix(1,2)=strten_in(6)
       str_matrix(2,1)=strten_in(6)
       str_matrix(1,3)=strten_in(5)
       str_matrix(3,1)=strten_in(5)
       str_matrix(2,3)=strten_in(4)
       str_matrix(3,2)=strten_in(4)
       call getvol(latvec_in,vol)
       transformed(:,1)=latvec_in(1,:)
       transformed(:,2)=latvec_in(2,:)
       transformed(:,3)=latvec_in(3,:)
       call invertmat(transformed,transformed_inv,3)
       flat=(-vol*matmul(str_matrix,transformed_inv))
       pressure=parini%target_pressure_habohr
       call stress_volume(latvec_in,vol,pressure,stressvol)
       flat=flat+stressvol
!Finally, write those values into fxyz
        do iat=1,3
          force_all(3*parini%nat+(iat-1)*3+1:3*parini%nat+iat*3)=flat(:,iat)
        enddo
!Get the enthalpy
        call get_enthalpy(latvec_in,etot_in,pressure,enthalpy)
end subroutine

subroutine get_BFGS_forces_atom(parini,parres,pos,force,latvec,enthalpy,getwfk,iprec,latvec_in,xred_in,etot_in,fcart_in,strten_in)
!This routine hides away all cumbersome conversion of arrays in lattice and positions and forces and stresses
!such that they can be directly passed on to bfgs. It also outputs the enthalpy instead of the energy
use interface_code

use mod_parini, only: typ_parini
implicit none
type(typ_parini), intent(in):: parini
type(typ_parini), intent(inout):: parres
integer:: iprec,iat
real(8):: pos(3*parini%nat),latvec(3,3)
real(8):: force(3*parini%nat)
real(8):: enthalpy,pressure,vol
real(8):: xred_in(3,parini%nat),fcart_in(3*parini%nat),strten_in(6),etot_in,latvec_in(3,3),transformed(3,3),transformed_inv(3,3)
real(8):: str_matrix(3,3),flat(3,3),pressure_mat(3,3),tmplat(3,3),sigma(3,3),crossp(3),stressvol(3,3)
logical:: getwfk

!Here we perform the force call
!****************************************************************************************************************        
!****************************************************************************************************************        
!Get the current acell_in and rprim_in, and also the atomic positions
!Generate a set of variables containing all degrees of freedome
    call rxyz_cart2int(latvec,xred_in,pos,parini%nat)
    latvec_in=latvec
       call get_energyandforces_single(parini,parres,latvec_in,xred_in,fcart_in,strten_in,etot_in,iprec,getwfk)
!****************************************************************************************************************   
!Conversion of forces is more complicate: 
!start with atomic forces
       force=fcart_in
!Get the enthalpy
        call get_enthalpy(latvec_in,etot_in,pressure,enthalpy)
end subroutine

subroutine get_BFGS_forces_lattice(parini,parres,pos,force,latvec,enthalpy,getwfk,iprec,latvec_in,xred_in,etot_in,fcart_in,strten_in)
!This routine hides away all cumbersome conversion of arrays in lattice and positions and forces and stresses
!such that they can be directly passed on to bfgs. It also outputs the enthalpy instead of the energy
 use interface_code

 use mod_parini, only: typ_parini
implicit none
type(typ_parini), intent(in):: parini
type(typ_parini), intent(inout):: parres
integer:: iprec,iat
real(8):: pos(3,parini%nat),latvec(3,3)
real(8):: force(9)
real(8):: enthalpy,pressure,vol
real(8):: xred_in(3,parini%nat),fcart_in(3*parini%nat),strten_in(6),etot_in,latvec_in(3,3),transformed(3,3),transformed_inv(3,3)
real(8):: str_matrix(3,3),flat(3,3),pressure_mat(3,3),tmplat(3,3),sigma(3,3),crossp(3),stressvol(3,3)
logical:: getwfk

!Here we perform the force call
!****************************************************************************************************************        
!****************************************************************************************************************        
!Get the current acell_in and rprim_in, and also the atomic positions
!Generate a set of variables containing all degrees of freedome
    xred_in=pos
    latvec_in=latvec
       call get_energyandforces_single(parini,parres,latvec_in,xred_in,fcart_in,strten_in,etot_in,iprec,getwfk)
!****************************************************************************************************************   
!Conversion of forces is more complicate: 
!Convert strten into "stress", which is actually the forces on the cell vectors
       str_matrix(1,1)=strten_in(1)
       str_matrix(2,2)=strten_in(2)
       str_matrix(3,3)=strten_in(3)
       str_matrix(1,2)=strten_in(6)
       str_matrix(2,1)=strten_in(6)
       str_matrix(1,3)=strten_in(5)
       str_matrix(3,1)=strten_in(5)
       str_matrix(2,3)=strten_in(4)
       str_matrix(3,2)=strten_in(4)
       call getvol(latvec_in,vol)
       transformed(:,1)=latvec_in(1,:)
       transformed(:,2)=latvec_in(2,:)
       transformed(:,3)=latvec_in(3,:)
       call invertmat(transformed,transformed_inv,3)
       flat=(-vol*matmul(str_matrix,transformed_inv))
       pressure=parini%target_pressure_habohr
       call stress_volume(latvec_in,vol,pressure,stressvol)
       flat=flat+stressvol
!Finally, write those values into fxyz
        do iat=1,3
          force((iat-1)*3+1:iat*3)=flat(:,iat)
        enddo
!Get the enthalpy
        call get_enthalpy(latvec_in,etot_in,pressure,enthalpy)
end subroutine


subroutine sd_minhocao(nat,nr,x,epot,f,betax,betax_lat,fmax,iter)
!Simple SD for testing purpose
implicit none
integer:: nr,i,iter,nat
real(8):: x(nr),f(nr),epot,betax,betax_lat,fmax

do i=1,nr-9
  x(i)=x(i)+betax*f(i)
enddo
do i=nr-8,nr
  x(i)=x(i)+betax_lat*f(i)
enddo

end subroutine

        subroutine stress_volume(latvec,vol,pressure,stressvol)
        !This subroutine will compute the additional components to the negative
        !derivative of the enthalpy with respect to the cell variables
        !For the derivatives with respect to hij, the lattive vector components of the
        !lattece-matrix h, the relation on http://en.wikipedia.org/wiki/Determinant is used:
        !\frac{\partial \det(h)}{\partial h_{ij}}= \det(h)(h^{-1})_{ji}. 
        implicit none
        integer:: nat,i,j
        real(8):: latvec(3,3),vol,stressvol(3,3),inv_latvec(3,3),pressure
        stressvol=0.d0
        call invertmat(latvec,inv_latvec,3)
        do i=1,3
           stressvol(i,1)=-vol*inv_latvec(1,i)*pressure
           stressvol(i,2)=-vol*inv_latvec(2,i)*pressure
           stressvol(i,3)=-vol*inv_latvec(3,i)*pressure
        enddo
        end subroutine


subroutine get_fmax(parini,fcart_in,strten_in,fmax,fmax_at,fmax_lat)
use mod_parini, only: typ_parini
implicit none
type(typ_parini), intent(in):: parini
integer:: iat,i,istr
real(8):: fcart_in(3,parini%nat),strten_in(6),fmax,fmax_at,fmax_lat
real(8):: dstr(6), strtarget(6)

!!Compute maximal component of forces, EXCLUDING any fixed components
 fmax=0.d0
 fmax_at=0.d0
 fmax_lat=0.d0
 do iat=1,parini%nat
   do i=1,3
!     if (dtsets(1)%iatfix(i,iat) /= 1) then
       if( abs(fcart_in(i,iat)) >= fmax_at ) fmax_at=abs(fcart_in(i,iat))
!     end if
   end do
 end do
 strtarget=0.d0
 strtarget(1:3)=-parini%target_pressure_habohr
 dstr(:)=strten_in(:)-strtarget(:)
!Eventually take into account the stress
 do istr=1,6
     if(abs(dstr(istr))*parini%paropt_geopt%strfact >= fmax_lat ) fmax_lat=abs(dstr(istr))*parini%paropt_geopt%strfact
 end do
 fmax=max(fmax_at,fmax_lat)
end subroutine
!******************************************
!                                         *
!GEOMETRY OPTIMIZERS MAX-BFGS             *                   
!                                         *
!******************************************



!************************************************************************************
subroutine init_hessinv(parini,hessin,latvec,omega,b0,lattdeg) 
!This routine will setup an inverse hessian accoprding to Pfrommer et al, J. Comp. Phys 131, 233 1997
!The hessin is in atomic units, taking omega in THZ and B0 in GPa (phonon frequency and bulk modulus) in
use mod_parini, only: typ_parini
use defs_basis
use mbfgs_interface
implicit none
type(typ_parini), intent(in):: parini
integer:: itype,iat,i,j,k,lattdeg
real(8):: omega,b0,hessin(3*parini%nat+9,3*parini%nat+9),diagat,avmass,diaglat
real(8):: amass(parini%nat),rcov,amass_u(parini%ntypat_global),vol
real(8),dimension(3,3):: diagat_lat,diagat_lat_inv,latvec,latvectrans
character(2):: tmp_ch
 write(*,'(a)') " # BFGS: initiallizing hessian"
call unit_matrix(hessin,3*parini%nat+9) !Initialize inverse Hessian to the unit matrix.
!Get the correct atomic masses and atomic character
 do itype=1,parini%ntypat_global
   call atmdata(amass_u(itype),rcov,tmp_ch,parini%znucl(itype))
 enddo
!Assign masses to each atom (for MD)
 do iat=1,parini%nat
   amass(iat)=amu_emass*amass_u(parini%typat_global(iat))
   write(*,'(a,i5,2(1x,es15.7))') " # BFGS: iat, AMU, EM: ", iat, amass_u(parini%typat_global(iat)),amass(iat)
 enddo
!Average mass
 avmass=sum(amass)/real(parini%nat,8)
 write(*,'(a,(1x,es15.7))') " # BFGS: average atomic mass: ",avmass
!Initialize Hessian diagonal elements of the cell part
 call getvol(latvec,vol)
if(lattdeg==1) diaglat=1.d0/(b0)*HaBohr3_GPa/vol**(1.d0/3.d0)  !Appropriate if one operates on lattice vectors
if(lattdeg==2) diaglat=1.d0/(3.d0*vol*b0)*HaBohr3_GPa          !Appropriate if one operates on strain
 do i=3*parini%nat+1,3*parini%nat+9
    hessin(i,i)=diaglat
 enddo 
!Intiallize the atomic part of the hessian
 diagat=1.d0/(avmass*(omega*Time_Sec*1.d12)**2)
 latvectrans(:,1)=latvec(1,:)
 latvectrans(:,2)=latvec(2,:)
 latvectrans(:,3)=latvec(3,:)
 diagat_lat=matmul(latvectrans,latvec)
 call invertmat(diagat_lat,diagat_lat_inv,3)
 diagat_lat=diagat_lat_inv*diagat
 do i=0,3*parini%nat-1,3
   do j=1,3
   do k=1,3
    hessin(i+j,i+k)=diagat_lat(j,k)
   enddo 
   enddo 
 enddo
! write(*,*) diagat_lat(1,:)
! write(*,*) diagat_lat(2,:)
! write(*,*) diagat_lat(3,:)
! do i=1,3*nat+9
! write(*,*) hessin(i,i)
! enddo
end subroutine init_hessinv
!************************************************************************************
subroutine GEOPT_MBFGS_MHM(parini,parres,latvec_in,xred_in,fcart_in,strten_in,etot_in,iprec,counter,folder)
!subroutine bfgs_driver_atoms(latvec_in,xred_in,fcart_in,strten_in,etot_in,iprec,counter,fmax_tol)
 use global, only: units
 use defs_basis
 use minpar
 use mod_fire,   only:dtmin, dtmax
!SUBROUTINE dfpmin_pos(nat,latvec,rxyz,fxyz,stress,pressure,etot,fnrmtol,iter,count)
use mbfgs_interface
use mod_parini, only: typ_parini
IMPLICIT NONE
type(typ_parini), intent(in):: parini
type(typ_parini), intent(inout):: parres
!REAL(8), INTENT(IN) :: fnrmtol!gtol 
REAL(8) :: fret, counter
REAL(8), INTENT(INOUT) :: xred_in(3*parini%nat),latvec_in(9),fcart_in(3*parini%nat),strten_in(6),etot_in
INTEGER, PARAMETER :: ITMAX=4000
REAL(8), PARAMETER :: STPMX=1.0d0,EPS=epsilon(xred_in),TOLX=4.0d0*EPS
!   Given a starting point p that is a vector of length N , the Broyden-Fletcher-Goldfarb-Shanno
!   variant of Davidon-Fletcher-Powell minimization is performed on a function func, using its
!   gradient as calculated by a routine dfunc. The convergence requirement on zeroing the
!   gradient is input as gtol. Returned quantities are p (the location of the minimum), iter
!   (the number of iterations that were performed), and fret (the minimum value of the
!   function). The routine lnsrch is called to perform approximate line minimizations.
!   Parameters: ITMAX is the maximum allowed number of iterations; STPMX is the scaled
!   maximum step length allowed in line searches; EPS is the machine precision; TOLX is the
!   convergence criterion on x values.
INTEGER :: its,i,j,info,LWORK
LOGICAL :: check
REAL(8) :: den,fac,fad,fae,fp,stpmax,sumdg,sumxi
REAL(8):: dg(3*parini%nat+9),g(3*parini%nat+9),hdg(3*parini%nat+9),pnew(3*parini%nat+9),xi(3*parini%nat+9),p(3*parini%nat+9)
REAL(8):: latvec0(9)
REAL(8):: tp(3*parini%nat+9),tg(3*parini%nat+9),dvin(3*parini%nat+9),vout(3*parini%nat+9),vout_prev(3*parini%nat+9)
REAL(8):: vin_min(3*parini%nat+9),vin(3*parini%nat+9),vel_in(3*parini%nat),vel_lat_in(9)
REAL(8):: vout_min(3*parini%nat+9),dedv_min(3*parini%nat+9)
REAL(8), DIMENSION(3*parini%nat+9,3*parini%nat+9) :: hessin,hessin0,hess_tmp,hess,hessin_dsyev
REAL(8) :: alpha_pl,dlatvec(9),dxred(3*parini%nat)
REAL(8) :: gamma0,gammax,lambda_1,lambda_2,tfp,dedv_1,dedv_2,etotal_1,etotal_2,dedv_predict
REAL(8) :: d2edv2_1,d2edv2_2,d2edv2_predict,etotal_predict,lambda_predict,lambda_predict_prev
INTEGER :: choice,status,sumstatus,iprec,iexit,lattdeg,hessupdate
REAL(8) :: rxyz0(3*parini%nat),eval(3*parini%nat+9),fmax,fmax_at,fmax_lat,pressure
LOGICAL :: getwfk
REAL(8) :: ent_pos_0,enthalpy,en0000,vvol_in
character(40)::filename,folder
character(4) ::fn4
logical:: multiprec
real(8):: tolmxf_switch,lambda,vol0,vol1,vol2 !,tolmxf0
real(8), allocatable:: WORK(:)
!real(8):: vabs
!real(8):: outerprod
type(typ_parini):: parini_tmp
!multiprec is hardcoded and, if true, starts a geopt with iprec==2, and then switches 
!to iprec==1 when the fmax==tolmxf_switch. The switch only occurs once
 multiprec=.true.
 tolmxf_switch=10.d0*parini%paropt_geopt%fmaxtol
 counter=0.d0

!lattdeg: this variable defines if the lattice degrees of freedom are treated directly or through the strain
 lattdeg=2    !1 for direct lattice coordinates, 2 for strain
 if(any(parini%fixlat)) lattdeg=1
!hessupdate: option for updating the approximate hessian: either the inverse hessian is updated, or the hessian itselfe (more costly, requires inversion of matrix)
 hessupdate=1 !1 for inverse hessian update, 2 for direct hessian update

!If the initial forces are too large, FIRE should be run initially to get down the high energy components
iprec=2
parini_tmp=parini
!tolmxf0=parini%paropt_geopt%fmaxtol
parini_tmp%paropt_geopt%fmaxtol=5.d-2
vel_in=0.d0;vel_lat_in=0.d0;vvol_in=0.d0
!Some conservative time steps for FIRE
parini_tmp%paropt_geopt%dt_start=10.d0;dtmin=1.d0;dtmax=50.d0
call GEOPT_FIRE_MHM(parini_tmp,parres,latvec_in,xred_in,fcart_in,strten_in,vel_in,vel_lat_in,vvol_in,etot_in,iprec,counter,folder)
!parini_tmp%paropt_geopt%fmaxtol=tolmxf0

!Now start the real BFGS
iexit=0
write(*,'(a,es15.7,es15.7)') " # BFGS OPTICAL FREQUENCY IN THZ, BULK MODULUS IN GPA: ", parmin_bfgs%betax, parmin_bfgs%betax_lat
latvec0=latvec_in
call getvol(latvec0,vol0)

lambda_predict_prev=1.d0
 
pressure=parini%target_pressure_habohr

open(unit=16,file="geopt.mon")
alpha_pl=1.d-0
gammax=1.d0
fret=0.d0
check=.true.
!call rxyz_cart2int(latvec,p(1:3*nat),rxyz,nat)
p(1:3*parini%nat)=xred_in(:)
if(lattdeg==1) then 
p(3*parini%nat+1:3*parini%nat+9)=latvec_in(:)
elseif(lattdeg==2) then
if(any(parini%fixlat)) stop "Fixed cell parameters for BFGS with lattdeg=2 not yet implemented"
p(3*parini%nat+1:3*parini%nat+9)=0.d0
endif

!Setup initial inverse hessian 
!Assume that that betax is the average optical phonon frequency in THz, and betax_lat is the expected Bulk modulus in GPa
   call init_hessinv(parini,hessin,latvec_in,parmin_bfgs%betax,parmin_bfgs%betax_lat,lattdeg) 
if(hessupdate==2) then
   call invertmat(hessin,hess,3*parini%nat+9)
endif

!!Calculate starting function value and gradient.
getwfk=.false.
!iprec=1
!This call is only to map all variables correctly
fp=-1.d10
call  get_BFGS_forces_strainlatt(parini,parres,p,g,fp,getwfk,iprec,latvec0,lattdeg,latvec_in,xred_in,etot_in,fcart_in,strten_in)
call convcheck(parini,parini%nat,latvec_in,fcart_in,strten_in,parini%target_pressure_habohr,parini%paropt_geopt%strfact,fmax,fmax_at,fmax_lat,parini%paropt_geopt%fmaxtol,iexit)
   !Eliminate components not to be changed
   if(any(parini%fixlat)) call elim_fixed_lat(parini,p(3*parini%nat+1:3*parini%nat+9),g(3*parini%nat+1:3*parini%nat+9))
   if(any(parini%fixat))  call elim_fixed_at(parini,parini%nat,g(1:3*parini%nat))
!call get_fmax(fcart_in,strten_in,fmax,fmax_at,fmax_lat)
if(counter==0.d0) then
!MHM: Write output to file in every step***********************************
!INITIAL STEP, STILL THE SAME STRUCTURE AS INPUT
       write(*,*) "Pressure, Energy",pressure,etot_in
       ent_pos_0=fp
       en0000=fp-ent_pos_0
       write(fn4,'(i4.4)') 0
       filename=trim(folder)//"posgeopt."//fn4//".ascii"
       units=units
       write(*,*) "# Writing the positions in :",filename
       call write_atomic_file_ascii(parini,filename,parini%nat,units,xred_in,latvec_in,fcart_in,strten_in,&
            &parini%char_type(1:parini%ntypat_global),parini%ntypat_global,parini%typat_global,parini%fixat,parini%fixlat,etot_in,pressure,fp,en0000)
       write(*,'(a,i4,4(1x,es17.8),1x,es9.2,1x,i4)') " # GEOPT BFGS AC ",int(counter),fp,fmax,fmax_lat,fmax_at,0.d0,iprec
!*********************************************************************
endif
!   if(fmax.lt.parini%paropt_geopt%fmaxtol) iexit=1
   if(iexit==1) then
   write(*,'(a)') " # BFGS converged before entering optimization"
!   call wtpos_inter(nat,rxyz,latvec,555)
   RETURN 
   endif



!Initial iprec after running the first force call
 if(multiprec) iprec=2

!call energyandforces(nat,latvec,rxyz,fxyz,stress,pressure,fp,count)
write(16,*) "Initial energy",fp
!call wtpos_inter(nat,rxyz,latvec,500)
!call fxyz_cart2int(nat,fxyz,g(1:3*nat),latvec)
g(3*parini%nat+1:3*parini%nat+9)=g(3*parini%nat+1:3*parini%nat+9)*alpha_pl
g=-g
xi=-matmul(hessin,g)

!Main loop over the iterations.
! call wtpos_inter(nat,rxyz,latvec,500)
 sumstatus=0
LWORK=100*(3*parini%nat+9)
allocate(WORK(LWORK))
hessin0=hessin
do its=1,ITMAX
!call unit_matrix(hessin)
!hessin=hessin*1.d-2


!!!write(*,*) "Hessin"
!!!do i=1,3*nat+9
!!!  do j=1,3*nat+9
!!!  if(j.ne.3*nat+9) write(*,'(1x,es8.1)',advance="no") hessin(i,j)
!!!  if(j==3*nat+9) write(*,'(1x,es8.1)') hessin(i,j)
!!!  enddo
!!!enddo


!        hessin_dsyev=hessin
!        call DSYEV('V','L',3*nat+9,hessin_dsyev,3*nat+9,eval,WORK,LWORK,INFO)
!        if (info.ne.0) stop 'DSYEV'
!        write(*,*) '---   App. eigenvalues in a.u. -------------'
!        do j=1,3*nat+9
!        if(j.ne.3*nat+9) write(*,'(1x,es8.1)',advance="no") eval(j)
!        if(j==3*nat+9) write(*,'(1x,es8.1)') eval(j)
!        enddo

 lambda=1.0d0
 1001 continue
 do i=1,3*parini%nat+9
    tp(i)=p(i)+lambda*xi(i)
 enddo
!write(*,*) "tp",tp
!write(*,*) "xi",xi
!write(*,*) "p",p

 !Check if displacement was too violent
  
if(lattdeg==1) then
 call getvol(p(3*parini%nat+1:3*parini%nat+9),vol1)
 call getvol(tp(3*parini%nat+1:3*parini%nat+9),vol2)
elseif(lattdeg==2) then
 call getvol_strain(p(3*parini%nat+1:3*parini%nat+9),latvec0,vol1)
 call getvol_strain(tp(3*parini%nat+1:3*parini%nat+9),latvec0,vol2)
endif

 if(vabs(lambda*xi(1:3*parini%nat))*vabs(latvec_in)/real(parini%nat,8).gt.0.1d0.or.&
   &abs(vol1-vol2)/vol1.gt.0.1d0) then
   lambda=lambda*0.5d0
   goto 1001
 endif

 if(parini%usewf_geopt) then
     getwfk=.true.
 else
     getwfk=.false.
 endif
 if(iprec==2.and.multiprec.and.fmax.lt.tolmxf_switch) then
     getwfk=.false.
     iprec=1
 endif
 counter=counter+1.d0
 call  get_BFGS_forces_strainlatt(parini,parres,tp,tg,tfp,getwfk,iprec,latvec0,lattdeg,latvec_in,xred_in,etot_in,fcart_in,strten_in)
 call convcheck(parini,parini%nat,latvec_in,fcart_in,strten_in,parini%target_pressure_habohr,parini%paropt_geopt%strfact,fmax,fmax_at,fmax_lat,parini%paropt_geopt%fmaxtol,iexit)
   !Eliminate components not to be changed
   if(any(parini%fixlat)) call elim_fixed_lat(parini,tp(3*parini%nat+1:3*parini%nat+9),tg(3*parini%nat+1:3*parini%nat+9))
   if(any(parini%fixat))  call elim_fixed_at(parini,parini%nat,tg(1:3*parini%nat))
! call get_fmax(fcart_in,strten_in,fmax,fmax_at,fmax_lat)
!MHM: Write output to file in every step***********************************
       write(*,*) "Pressure, Energy",pressure,etot_in
       ent_pos_0=fp
       en0000=fp-ent_pos_0
       write(fn4,'(i4.4)') int(counter) ! (its)*2-1 
       filename=trim(folder)//"posgeopt."//fn4//".ascii"
       units=units
       write(*,*) "# Writing the positions in :",filename
       call write_atomic_file_ascii(parini,filename,parini%nat,units,xred_in,latvec_in,fcart_in,strten_in,&
            &parini%char_type(1:parini%ntypat_global),parini%ntypat_global,parini%typat_global,parini%fixat,parini%fixlat,etot_in,pressure,fp,en0000)
       write(*,'(a,i4,4(1x,es17.8),1x,es9.2,1x,i4)') " # GEOPT BFGS LS ",int(counter),tfp,fmax,fmax_lat,fmax_at,lambda,iprec
!*********************************************************************
 tg=-tg
! if(fmax.lt.parini%paropt_geopt%fmaxtol) then
!   iexit=1
 if(iexit==1) then
   fp=tfp
   goto 1002
 endif
 
 choice=1
 lambda_1=lambda   ; lambda_2=0.0d0
 etotal_1=tfp      ; etotal_2=fp
! dvin(:)=vin(:)-vin_prev(:)
 dvin(:)=tp(:)-p(:)
 vout=tg
 vout_prev=g
! vout_prev=xi
 dedv_1=dot_product(vout,dvin)
 dedv_2=dot_product(vout_prev,dvin)
 call findmin(choice,dedv_1,dedv_2,dedv_predict,&
& d2edv2_1,d2edv2_2,d2edv2_predict,&
& etotal_1,etotal_2,etotal_predict,&
& lambda_1,lambda_2,lambda_predict,status)


!lambda_predict_prev=min(lambda_predict_prev,lambda)

if(lambda_predict.lt.0.5d0*lambda.or.lambda_predict.gt.1.5d0*lambda) then
   if(lambda_predict.lt.0.d0) lambda_predict=lambda*0.5d0
   lambda_predict=max(1.d-6,min(min(lambda_predict,5.d0),0.5d0*(lambda_predict+lambda_predict_prev)))
!   pnew(:)=p(:)+lambda_predict*xi(:)
   dlatvec=lambda_predict*xi(3*parini%nat+1:3*parini%nat+9)
   dxred=lambda_predict*xi(1:3*parini%nat)
   call propagate(parini,parini%nat,p(1:3*parini%nat),p(3*parini%nat+1:3*parini%nat+9),dxred,dlatvec,pnew(1:3*parini%nat),pnew(3*parini%nat+1:3*parini%nat+9))
   if(parini%usewf_geopt) then
       getwfk=.true.
   else
       getwfk=.false.
   endif
   counter=counter+1.d0
   dg=g       !Save the old gradient,
   call  get_BFGS_forces_strainlatt(parini,parres,pnew,g,fp,getwfk,iprec,latvec0,lattdeg,latvec_in,xred_in,etot_in,fcart_in,strten_in)
   call convcheck(parini,parini%nat,latvec_in,fcart_in,strten_in,parini%target_pressure_habohr,parini%paropt_geopt%strfact,fmax,fmax_at,fmax_lat,parini%paropt_geopt%fmaxtol,iexit)
      !Eliminate components not to be changed
      if(any(parini%fixlat)) call elim_fixed_lat(parini,pnew(3*parini%nat+1:3*parini%nat+9),g(3*parini%nat+1:3*parini%nat+9))
      if(any(parini%fixat))  call elim_fixed_at(parini,parini%nat,g(1:3*parini%nat))
!   call get_fmax(fcart_in,strten_in,fmax,fmax_at,fmax_lat)
!MHM: Write output to file in every step***********************************
       write(*,*) "Pressure, Energy",pressure,etot_in
       ent_pos_0=fp
       en0000=fp-ent_pos_0
       write(fn4,'(i4.4)') int(counter)! its*2
       filename=trim(folder)//"posgeopt."//fn4//".ascii"
       units=units
       write(*,*) "# Writing the positions in :",filename
       call write_atomic_file_ascii(parini,filename,parini%nat,units,xred_in,latvec_in,fcart_in,strten_in,&
            &parini%char_type(1:parini%ntypat_global),parini%ntypat_global,parini%typat_global,parini%fixat,parini%fixlat,etot_in,pressure,fp,en0000)
       write(*,'(a,i4,4(1x,es17.8),1x,es9.2,1x,i4)') " # GEOPT BFGS AC ",int(counter),fp,fmax,fmax_lat,fmax_at,lambda_predict,iprec
!*********************************************************************
       g=-g       !New gradient 
   else
      lambda_predict=lambda
   !   pnew=p+lambda_predict*xi
      dlatvec=lambda_predict*xi(3*parini%nat+1:3*parini%nat+9)
      dxred=lambda_predict*xi(1:3*parini%nat)
      call propagate(parini,parini%nat,p(1:3*parini%nat),p(3*parini%nat+1:3*parini%nat+9),dxred,dlatvec,pnew(1:3*parini%nat),pnew(3*parini%nat+1:3*parini%nat+9))
      if(tp(5).ne.pnew(5)) stop "Womething srong!!!"
      fp=tfp
      dg=g       !Save the old gradient,
      g=tg      !New gradient 
   endif
!Update the line direction, and the current point.
   xi=pnew-p
   p=pnew
   lambda_predict_prev=lambda_predict

1002 continue
   write(16,'(a,1x,I5,1x,1pe21.14,1x,1pe12.5,1x,1pe12.5,1x,1pe12.5)') "  &
   &  BFGS_all",int(counter),fp,fmax,fmax_at,fmax_lat
!   if(fmax.lt.parini%paropt_geopt%fmaxtol) iexit=1
!   if (fnrm < fnrmtol) then  !Test for convergence on zero gradient.
!   latvec=p(3*nat+1:3*nat+9)
!   call backtocell(nat,latvec,rxyz)
   if(iexit==1) then
   write(*,'(a,i4,2(1x,es25.15))') " # BFGS converged", int(counter),fp,fmax
!   call wtpos_inter(nat,rxyz,latvec,555)
   RETURN 
   endif
   if(int(counter).gt.parini%paropt_geopt%nit) then
   write(*,'(a,i5)') " # BFGS did not converge in steps: ", int(counter)
   RETURN
   endif

if(hessupdate==1) then
!Update scheme from wikipedia
! dg=g-dg                !Compute difference of gradients,
! hdg=matmul(hessin,dg)  !and difference times current matrix.
! fac=dot_product(dg,xi) !TERM1
! fae=dot_product(dg,hdg)!TERM2
! sumdg=dot_product(dg,dg)
! sumxi=dot_product(xi,xi)
!  if (fac > sqrt(EPS*sumdg*sumxi)) then !Skip update if fac not sufficiently positive.
!      hessin=hessin-1.d0/fac*(outerprod(hdg,xi)+matmul(outerprod(xi,dg),hessin)) !LAST TERM FIRST
!      hessin=hessin+1.d0/fac**2*(fac+fae)*outerprod(xi,xi)    !FROM WIKI
!!      hessin=hessin+(1.d0+fae/fac)*outerprod(xi,xi)/fac      !FROM EDWIN CHONG
!  else
!  write(*,'(a,4(1x,es15.7))') "WARNING!!!",eps,fac,sumdg,sumxi
!  call init_hessinv(hessin,latvec_in,parmin_bfgs%betax,parmin_bfgs%betax_lat,lattdeg) 
!  endif


!!NOCEDAL FORM************
!dg=g-dg                !Compute difference of gradients,
!call unit_matrix(hess_tmp)
!fac=dot_product(dg,xi)
!hess_tmp=hess_tmp-outerprod(xi,dg)/fac
!hessin=matmul(hess_tmp,hessin)
!hessin=matmul(hessin,hess_tmp)
!hessin=hessin+outerprod(xi,xi)/fac
!!NOCEDAL FORM************


!Update scheme from numerical recipes
   dg=g-dg                !Compute difference of gradients,
   hdg=matmul(hessin,dg)  !and difference times current matrix.
   fac=dot_product(dg,xi) ! Calculate dot products for the denominators.
   fae=dot_product(dg,hdg)
   sumdg=dot_product(dg,dg)
   sumxi=dot_product(xi,xi)
    if (fac > sqrt(EPS*sumdg*sumxi)) then !Skip update if fac not sufficiently positive.
        fac=1.0d0/fac
        fad=1.0d0/fae
        dg=fac*xi-fad*hdg                 !Vector that makes BFGS different from DFP.
        hessin=hessin+fac*outerprod(xi,xi)
        hessin=hessin-fad*outerprod(hdg,hdg)
        hessin=hessin+fae*outerprod(dg,dg)
    else
    write(*,'(a,4(1x,es15.7))') "WARNING!!!",eps,fac,sumdg,sumxi
    call init_hessinv(parini,hessin,latvec_in,parmin_bfgs%betax,parmin_bfgs%betax_lat,lattdeg) 
    end if

elseif(hessupdate==2) then
!Update the hessian, then invert it to get hessin: from WIKI
   dg=g-dg              !Compute difference of gradients,
   hdg=matmul(hess,xi)  !and difference times current matrix.
   hess_tmp=outerprod(xi,xi)
   hess_tmp=matmul(hess,hess_tmp)
   hess_tmp=matmul(hess_tmp,hess)
   fae=dot_product(dg,hdg)
   fac=dot_product(dg,xi) ! Calculate dot products for the denominators.
   hess=hess+outerprod(dg,dg)/fac-hess_tmp/fae
   sumdg=dot_product(dg,dg)
   sumxi=dot_product(xi,xi)
   if (fac > sqrt(EPS*sumdg*sumxi)) then !Skip update if fac not sufficiently positive.
   call invertmat(hess,hessin,3*parini%nat+9)
   else
   write(*,'(a,4(1x,es15.7))') "WARNING!!!",eps,fac,sumdg,sumxi
   call init_hessinv(parini,hessin,latvec_in,parmin_bfgs%betax,parmin_bfgs%betax_lat,lattdeg) 
   endif
!*************************
endif
!if(modulo(int(counter),int((nat*3+9)*0.5d0))==0) then
!Assume that that betax is the average optical phonon frequency in THz, and betax_lat is the expected Bulk modulus in GPa
!call init_hessinv(hessin,latvec_in,parmin_bfgs%betax,parmin_bfgs%betax_lat,lattdeg) 
!call invertmat(hessin,hess,3*nat+9)
!endif
call correct_hessin(hess,hessin,latvec_in,3*parini%nat+9,hessupdate,lattdeg)
!Now calculate the next direction to go
   xi=-matmul(hessin,g)
   write(*,'(a,9(1x,es10.3))') "strain ",p(3*parini%nat+1:3*parini%nat+9)
!and go back for another iteration.
end do
write(*,*) "Too many iterations"

!stop "Too many iterations"
END SUBROUTINE

!************************************************************************************
!************************************************************************************
subroutine GEOPT_MBFGS_MHM_OLD(parini,parres,latvec_in,xred_in,fcart_in,strten_in,etot_in,iprec,counter,folder)
!subroutine bfgs_driver_atoms(latvec_in,xred_in,fcart_in,strten_in,etot_in,iprec,counter,fmax_tol)
 use global, only: units
 use defs_basis
 use minpar

!SUBROUTINE dfpmin_pos(nat,latvec,rxyz,fxyz,stress,pressure,etot,fnrmtol,iter,count)
use mbfgs_interface
use mod_parini, only: typ_parini
implicit none
type(typ_parini), intent(in):: parini
type(typ_parini), intent(inout):: parres
!REAL(8), INTENT(IN) :: fnrmtol!gtol 
REAL(8) :: fret, counter
REAL(8), INTENT(INOUT) :: xred_in(3*parini%nat),latvec_in(9),fcart_in(3*parini%nat),strten_in(6),etot_in
INTEGER, PARAMETER :: ITMAX=4000
REAL(8), PARAMETER :: STPMX=1.0d0,EPS=epsilon(xred_in),TOLX=4.0d0*EPS
!   Given a starting point p that is a vector of length N , the Broyden-Fletcher-Goldfarb-Shanno
!   variant of Davidon-Fletcher-Powell minimization is performed on a function func, using its
!   gradient as calculated by a routine dfunc. The convergence requirement on zeroing the
!   gradient is input as gtol. Returned quantities are p (the location of the minimum), iter
!   (the number of iterations that were performed), and fret (the minimum value of the
!   function). The routine lnsrch is called to perform approximate line minimizations.
!   Parameters: ITMAX is the maximum allowed number of iterations; STPMX is the scaled
!   maximum step length allowed in line searches; EPS is the machine precision; TOLX is the
!   convergence criterion on x values.
INTEGER :: its,i
LOGICAL :: check
REAL(8) :: den,fac,fad,fae,fp,stpmax,sumdg,sumxi
REAL(8):: dg(3*parini%nat+9),g(3*parini%nat+9),hdg(3*parini%nat+9),pnew(3*parini%nat+9),xi(3*parini%nat+9),p(3*parini%nat+9)
REAL(8):: tp(3*parini%nat+9),tg(3*parini%nat+9),dvin(3*parini%nat+9),vout(3*parini%nat+9),vout_prev(3*parini%nat+9)
REAL(8):: vin_min(3*parini%nat+9),vin(3*parini%nat+9)
REAL(8):: vout_min(3*parini%nat+9),dedv_min(3*parini%nat+9)
REAL(8), DIMENSION(3*parini%nat+9,3*parini%nat+9) :: hessin,hessin0
REAL(8) :: alpha_pl
REAL(8) :: gamma0,gammax,lambda_1,lambda_2,tfp,dedv_1,dedv_2,etotal_1,etotal_2,dedv_predict
REAL(8) :: d2edv2_1,d2edv2_2,d2edv2_predict,etotal_predict,lambda_predict,lambda_predict_prev
INTEGER :: choice,status,sumstatus,iprec,iexit
REAL(8) :: latvec0(9),rxyz0(3*parini%nat),eval(3*parini%nat+9),fmax,fmax_at,fmax_lat,pressure
LOGICAL :: getwfk
REAL(8) :: ent_pos_0,enthalpy,en0000
character(40)::filename,folder
character(4) ::fn4
logical:: multiprec
real(8):: tolmxf_switch 
!real(8):: vabs
!real(8):: outerprod
!multiprec is hardcoded and, if true, starts a geopt with iprec==2, and then switches 
!to iprec==1 when the fmax==tolmxf_switch. The switch only occurs once
 multiprec=.true.
 tolmxf_switch=10.d0*parini%paropt_geopt%fmaxtol

 counter=0.d0
write(*,'(a,es15.7,es15.7)') " # BFGS BETAX, BETAX_LAT: ", parmin_bfgs%betax, parmin_bfgs%betax_lat
 
pressure=parini%target_pressure_habohr

open(unit=16,file="geopt.mon")
alpha_pl=1.d-0
gammax=1.d0
fret=0.d0
check=.true.
!call rxyz_cart2int(latvec,p(1:3*nat),rxyz,nat)
p(1:3*parini%nat)=xred_in(:)
p(3*parini%nat+1:3*parini%nat+9)=latvec_in(:)
!Calculate starting function value and gradient.
getwfk=.false.
iprec=1
call get_BFGS_forces_max(parini,parres,p,g,fp,getwfk,iprec,latvec_in,xred_in,etot_in,fcart_in,strten_in)
call get_fmax(parini,fcart_in,strten_in,fmax,fmax_at,fmax_lat)
!MHM: Write output to file in every step***********************************
!INITIAL STEP, STILL THE SAME STRUCTURE AS INPUT
       write(*,*) "Pressure, Energy",pressure,etot_in
       ent_pos_0=fp
       en0000=fp-ent_pos_0
       write(fn4,'(i4.4)') 0
       filename=trim(folder)//"posgeopt."//fn4//".ascii"
       units=units
       write(*,*) "# Writing the positions in :",filename
       call write_atomic_file_ascii(parini,filename,parini%nat,units,xred_in,latvec_in,fcart_in,strten_in,&
            &parini%char_type(1:parini%ntypat_global),parini%ntypat_global,parini%typat_global,parini%fixat,parini%fixlat,etot_in,pressure,fp,en0000)
       write(*,'(a,i4,4(1x,es17.8),1x,es9.2,1x,i4)') " # GEOPT BFGS AC ",0,fp,fmax,fmax_lat,fmax_at,0.d0,iprec
!*********************************************************************
   iexit=0
   if(fmax.lt.parini%paropt_geopt%fmaxtol) iexit=1
   if(iexit==1) then
   write(*,'(a)') " # BFGS converged before entering optimization"
!   call wtpos_inter(nat,rxyz,latvec,555)
   RETURN 
   endif

!Initial iprec after running the first force call
 if(multiprec) iprec=2

!call energyandforces(nat,latvec,rxyz,fxyz,stress,pressure,fp,count)
write(16,*) "Initial energy",fp
!call wtpos_inter(nat,rxyz,latvec,500)
!call fxyz_cart2int(nat,fxyz,g(1:3*nat),latvec)
g(3*parini%nat+1:3*parini%nat+9)=g(3*parini%nat+1:3*parini%nat+9)*alpha_pl
g=-g
call unit_matrix(hessin,3*parini%nat+9) !Initialize inverse Hessian to the unit matrix.

!Initialize Hessian diagonal elements
hessin=hessin*parmin_bfgs%betax
!Initialize Hessian diagonal elements for the cell variables
do i=3*parini%nat+1,3*parini%nat+9
hessin(i,i)=parmin_bfgs%betax_lat
enddo

call init_hessinv(parini,hessin,latvec_in,parmin_bfgs%betax,parmin_bfgs%betax_lat,1) 
xi=-matmul(hessin,g)
!Main loop over the iterations.
! call wtpos_inter(nat,rxyz,latvec,500)
 sumstatus=0

hessin0=hessin
do its=1,ITMAX

 gamma0=1.d0*gammax
1001 continue
 do i=1,3*parini%nat+9
    tp(i)=p(i)+gamma0*xi(i)
 enddo
 write(16,*) "Length of movement along xi", vabs(gamma0*xi)
 if (vabs(gamma0*xi).gt.7.d-1) then
 gamma0=gamma0*0.5d0
 goto 1001
 endif

 if(parini%usewf_geopt) then
     getwfk=.true.
 else
     getwfk=.false.
 endif
 if(iprec==2.and.multiprec.and.fmax.lt.tolmxf_switch) then
     getwfk=.false.
     iprec=1
 endif
 counter=counter+1.d0
 call get_BFGS_forces_max(parini,parres,tp,tg,tfp,getwfk,iprec,latvec_in,xred_in,etot_in,fcart_in,strten_in)
 call get_fmax(parini,fcart_in,strten_in,fmax,fmax_at,fmax_lat)
!MHM: Write output to file in every step***********************************
       write(*,*) "Pressure, Energy",pressure,etot_in
       ent_pos_0=fp
       en0000=fp-ent_pos_0
       write(fn4,'(i4.4)') (its)*2-1 
       filename=trim(folder)//"posgeopt."//fn4//".ascii"
       units=units
       write(*,*) "# Writing the positions in :",filename
       call write_atomic_file_ascii(parini,filename,parini%nat,units,xred_in,latvec_in,fcart_in,strten_in,&
            &parini%char_type(1:parini%ntypat_global),parini%ntypat_global,parini%typat_global,parini%fixat,parini%fixlat,etot_in,pressure,fp,en0000)
       write(*,'(a,i4,4(1x,es17.8),1x,es9.2,1x,i4)') " # GEOPT BFGS LS ",(its)*2-1,tfp,fmax,fmax_lat,fmax_at,gamma0,iprec
!*********************************************************************
! call rxyz_int2cart(tp(3*nat+1:3*nat+9),tp(1:3*nat),rxyz,nat)
! call energyandforces(nat,tp(3*nat+1:3*nat+9),rxyz,fxyz,stress,pressure,tfp,count)
! call wtpos(nat,tp(3*nat+1:3*nat+9),rxyz,real(its,8))
! call fxyz_cart2int(nat,fxyz,tg,tp(3*nat+1:3*nat+9))
! tg(3*nat+1:3*nat+9)=stress(:)*alpha_pl


 tg=-tg 
 if (fmax.lt.1.d-3.or.sumstatus.gt.10000)then
 choice=1
 else
 choice=4
 endif 
 1000 continue
 lambda_1=1.0d0       ; lambda_2=0.0d0
 etotal_1=tfp      ; etotal_2=fp
! dvin(:)=vin(:)-vin_prev(:)
 dvin(:)=tp(:)-p(:)
 vout=tg
 vout_prev=g
! vout_prev=xi
 dedv_1=dot_product(vout,dvin)
 dedv_2=dot_product(vout_prev,dvin)
 call findmin(choice,dedv_1,dedv_2,dedv_predict,&
& d2edv2_1,d2edv2_2,d2edv2_predict,&
& etotal_1,etotal_2,etotal_predict,&
& lambda_1,lambda_2,lambda_predict,status)

if(status==2) sumstatus=sumstatus+status
if(status==3) then
choice=1
call findmin(choice,dedv_1,dedv_2,dedv_predict,&
& d2edv2_1,d2edv2_2,d2edv2_predict,&
& etotal_1,etotal_2,etotal_predict,&
& lambda_1,lambda_2,lambda_predict,status)
endif

!NEW: limit the step to go to 1.d0
lambda_predict=min(lambda_predict,3.5d0)
lambda_predict=max(lambda_predict,-1.d0)


!Generates vin at the minimum, and an interpolated vout, modified
!to have the right value of dedv_predict, from findmin.
! vin_min(:)=vin_prev(:)+lambda_predict*dvin(:)
 vin_min(:)=p(:)+lambda_predict*dvin(:)
 vout_min(:)=vout_prev(:)+lambda_predict*(vout(:)-vout_prev(:))
 dedv_min=dedv_2+lambda_predict*(dedv_1-dedv_2)
!Modify vout_min in order to impose dedv_predict
 vout_min(:)=vout_min(:)+dvin(:)*(dedv_predict-dedv_min)/dot_product(dvin,dvin)
 vin(:)=vin_min(:)
 pnew(:)=vin(:)

!   Update the line direction, and the current point.
   xi=pnew-p
   p=pnew
   dg=g       !Save the old gradient,

   if(parini%usewf_geopt) then
       getwfk=.true.
   else
       getwfk=.false.
   endif
   counter=counter+1.d0
   call get_BFGS_forces_max(parini,parres,p,g,fp,getwfk,iprec,latvec_in,xred_in,etot_in,fcart_in,strten_in)
   call get_fmax(parini,fcart_in,strten_in,fmax,fmax_at,fmax_lat)
!MHM: Write output to file in every step***********************************
       write(*,*) "Pressure, Energy",pressure,etot_in
       ent_pos_0=fp
       en0000=fp-ent_pos_0
       write(fn4,'(i4.4)')  its*2
       filename=trim(folder)//"posgeopt."//fn4//".ascii"
       units=units
       write(*,*) "# Writing the positions in :",filename
       call write_atomic_file_ascii(parini,filename,parini%nat,units,xred_in,latvec_in,fcart_in,strten_in,&
            &parini%char_type(1:parini%ntypat_global),parini%ntypat_global,parini%typat_global,parini%fixat,parini%fixlat,etot_in,pressure,fp,en0000)
       write(*,'(a,i4,4(1x,es17.8),1x,es9.2,1x,i4)') " # GEOPT BFGS AC ",its*2,fp,fmax,fmax_lat,fmax_at,lambda_predict,iprec
!*********************************************************************
!   call rxyz_int2cart(p(3*nat+1:3*nat+9),p(1:3*nat),rxyz,nat)
!   call energyandforces(nat,p(3*nat+1:3*nat+9),rxyz,fxyz,stress,pressure,etot,count)
!   call fxyz_cart2int(nat,fxyz,g,p(3*nat+1:3*nat+9))
!   g(3*nat+1:3*nat+9)=stress(:)*alpha_pl

!   fp=etot  
   g=-g
   den=max(fret,1.0d0)
   write(16,'(a,1x,I5,1x,1pe21.14,1x,1pe12.5,1x,1pe12.5,1x,1pe12.5)') "  &
   &  BFGS_all",its,fp,fmax,fmax_at,fmax_lat
   iexit=0
   if(fmax.lt.parini%paropt_geopt%fmaxtol) iexit=1
!   if (fnrm < fnrmtol) then  !Test for convergence on zero gradient.
!   latvec=p(3*nat+1:3*nat+9)
!   call backtocell(nat,latvec,rxyz)
   if(iexit==1) then
   write(*,'(a,i5,a)') " # BFGS converged in ",its*2," iterations"
!   call wtpos_inter(nat,rxyz,latvec,555)
   RETURN 
   endif
   if(int(counter).gt.parini%paropt_geopt%nit) then
   write(*,'(a,i5)') " # BFGS did not converg in steps: ", int(counter)
   RETURN
   endif

   dg=g-dg                !Compute difference of gradients,
   hdg=matmul(hessin,dg)  !and difference times current matrix.
   fac=dot_product(dg,xi) ! Calculate dot products for the denominators.
   fac=dot_product(dg,xi)
   fae=dot_product(dg,hdg)
   sumdg=dot_product(dg,dg)
   sumxi=dot_product(xi,xi)
    if (fac > sqrt(EPS*sumdg*sumxi)) then !Skip update if fac not sufficiently positive.
        fac=1.0d0/fac
        fad=1.0d0/fae
        dg=fac*xi-fad*hdg                 !Vector that makes BFGS different from DFP.
        hessin=hessin+fac*outerprod(xi,xi)
        hessin=hessin-fad*outerprod(hdg,hdg)
        hessin=hessin+fae*outerprod(dg,dg)
    end if
!Now calculate the next direction to go
    xi=-matmul(hessin,g)
!and go back for another iteration.
end do
write(*,*) "Too many iterations"

!stop "Too many iterations"
END SUBROUTINE

!************************************************************************************

FUNCTION vabs(v) result(res)
implicit none
real(8),dimension(:):: v
real(8)::sumv
integer::i
real(8):: res
sumv=0.d0
do i=1,size(v)
sumv=sumv+v(i)*v(i)
enddo
res=sqrt(sumv)
!read(*,*) i
END FUNCTION vabs

!************************************************************************************

function outerprod(a,b)
implicit none
real(8),dimension(:),intent(in)::a,b
real(8),dimension(size(a),size(b))::outerprod
integer:: i,j
!outerprod=spread(a,dim=2,ncopies=size(b))*spread(b,dim=1,ncopies=size(a))
do i=1,size(a)
do j=1,size(b)
outerprod(i,j)=a(i)*b(j)
enddo
enddo
end function


!************************************************************************************

FUNCTION assert_eq(n1,n2,n3,n4,string) result(res)
implicit none
CHARACTER(LEN=*), INTENT(IN) :: string
INTEGER, INTENT(IN) :: n1,n2,n3,n4
INTEGER :: res
!Action:
!Embedding program dies gracefully with an error message if any of the
!integer arguments are not equal to the first. Otherwise, return the value of
!the first argument. Typical use is for enforcing equality on the sizes of arrays
!passed to a subprogram. nrutil implements and overloads forms with 1, 2,
!3, and 4 integer arguments.
if (n1==n2.and.n2==n3.and.n3==n4) then
res=n1
else
write (*,*) "error: an assert_eq failed with this tag:", string
STOP "program terminated by assert_eq"
end if
END FUNCTION assert_eq

!************************************************************************************

subroutine findmin(choice,dedv_1,dedv_2,dedv_predict,&
& d2edv2_1,d2edv2_2,d2edv2_predict,&
& etotal_1,etotal_2,etotal_predict,&
& lambda_1,lambda_2,lambda_predict,status)
!{\src2tex{textfont=tt}}
!!****f* ABINIT/findmin
!!
!! NAME
!! findmin
!!
!! FUNCTION
!! Compute the minimum of a function whose value
!! and derivative are known at two points,
!! using different algorithms.
!! Also deduce different quantities at this predicted
!! point, and at the two other points
!!
!! COPYRIGHT
!! Copyright (C) 1998-2009 ABINIT group (XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!

!! INPUTS
!! choice=1,uses a linear interpolation of the derivatives
!!       =2,uses a quadratic interpolation based on the
!!        values of the function, and the second derivative at mid-point
!!       =3,uses a cubic interpolation
!!       =4,uses a quartic interpolation, with the supplementary
!!          condition that the second derivative vanishes at one and
!!          only one point (See Schlegel, J. Comp. Chem. 3, 214 (1982).
!!          For this option, lambda_1 must be 1 (new point),
!!          and lambda_2 must be 0 (old point).
!!          Also, if the derivative at the new point is more negative
!!          than the derivative at the old point, the predicted
!!          point cannot correspond to a minimum, but will be lambda=2.5d0,
!!          if the energy of the second point is lower than the energy
!!          of the first point.
!! etotal_1=first value of the function
!! etotal_2=second value of the function
!! dedv_1=first value of the derivative
!! dedv_2=second value of the derivative
!! lambda_1=first value of the argument
!! lambda_2=second value of the argument
!!
!! OUTPUT
!! dedv_predict=predicted value of the derivative (usually zero,
!!  except if choice=4, if it happens that a minimum cannot be located,
!!  and a trial step is taken)
!! d2edv2_predict=predicted value of the second derivative (not if choice=4)
!! d2edv2_1=first value of the second derivative (not if choice=4)
!! d2edv2_2=second value of the second derivative (not if choice=4)
!! etotal_predict=predicted value of the function
!! lambda_predict=predicted value of the argument
!! status= 0 if everything went normally ;
!!         1 if negative second derivative
!!         2 if some other problem
!!!!My part:: if choice=4 and ee=0 then status=3
!!
!! PARENTS
!!      brdene,scfcge
!!
!! CHILDREN
!!      leave_new,wrtout
!!
!! SOURCE


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: choice
 integer,intent(out) :: status
 real(8),intent(in) :: dedv_1,dedv_2,etotal_1,etotal_2,lambda_1,lambda_2
 real(8),intent(out) :: d2edv2_1,d2edv2_2,d2edv2_predict,dedv_predict
 real(8),intent(out) :: etotal_predict,lambda_predict

!Local variables-------------------------------
!scalars
 integer :: printvol
 real(8) :: aa,bb,bbp,cc,ccp,d2edv2_mid,d_lambda,dd,dedv_2bis,dedv_mid1
 real(8) :: dedv_mid2,discr,ee,eep,etotal_2bis,lambda_shift,sum1,sum2,sum3,uu
 real(8) :: uu3,vv,vv3
 real(8) :: tol12
 character(len=500) ::message 

! *************************************************************************
 tol12=0.000000000001d0
 open(unit=16,file="geopt.mon")
!DEBUG
!write(6,*)' findmin : enter'
!write(6,*)' choice,lambda_1,lambda_2=',choice,lambda_1,lambda_2
!ENDDEBUG

 status=0
 d_lambda=lambda_1-lambda_2

!DEBUG
!do choice=3,1,-1
!ENDDEBUG

 if(choice==3)then

! Evaluate cubic interpolation
! etotal = aa + bb * lambda + cc * lambda**2 + dd * lambda**3
  dedv_mid1=(dedv_1+dedv_2)/2.0d0
  dedv_mid2=(etotal_1-etotal_2)/d_lambda
  d2edv2_mid=(dedv_1-dedv_2)/d_lambda
  dd=2.0d0 * ( dedv_mid1 - dedv_mid2 ) / (d_lambda**2)
  cc=0.5d0 * ( d2edv2_mid - 3.0d0*dd*(lambda_1+lambda_2) )
  bb=dedv_2 - 2*cc*lambda_2 - 3*dd*lambda_2**2
  aa=etotal_2 - bb*lambda_2 - cc*lambda_2**2 - dd*lambda_2**3

! Find the lambda at the minimum
  discr=cc*cc-3*bb*dd
  if(discr<0.0d0)then
   write(16, '(a,a,/,a)' )"  ",&
&   ' findmin : BUG -',&
&   '  The 2nd degree equation has no root (choice=3).'
!   call wrtout(06,*,'COLL')
!   call leave_new('COLL')
  end if
  discr=sqrt(discr)
! The root that gives a minimum corresponds to  +discr
  lambda_predict=(-cc+discr)/(3.0d0*dd)

! Predict etotal at that lambda
  etotal_predict=aa+lambda_predict*(bb&
&  +lambda_predict*(cc&
&  +lambda_predict* dd  ))
  dedv_predict=bb+2.0d0*cc*lambda_predict&
&  +3.0d0*dd*lambda_predict**2
  d2edv2_1=2*cc+6*dd*lambda_1
  d2edv2_2=2*cc+6*dd*lambda_2
  d2edv2_predict=2*cc+6*dd*lambda_predict

 else if(choice==4)then

  if(abs(lambda_1-1.0d0)>tol12 .or. abs(lambda_2)>tol12) then
   write(16, '(a,a,/,a)' )"  ",&
&   ' findmin : BUG -',&
&   '  For choice=4, lambda_1 must be 1 and lambda_2 must be 0.'
!   call wrtout(06,*,'COLL')
!   call leave_new('COLL')
  end if

! Evaluate quartic interpolation
! etotal = aa + bb * lambda + cc * lambda**2 + dd * lambda**3 + ee * lambda**4
! Impose positive second derivative everywhere, with
! one point where it vanishes :  3*dd**2=8*cc*ee
  aa=etotal_2
  bb=dedv_2
  sum1=etotal_1-aa-bb
  sum2=dedv_1-bb
  sum3=sum2-2.0d0*sum1

! Build the discriminant of the associated 2nd degree equation
  discr=sum2**2-3.0d0*sum3**2
  if(discr<0.0d0 .or. sum2<0.0d0)then

!  Even if there is a problem, try to keep going ...
   write(16, '(a,a,/,a)' )"  ",&
&   ' findmin : WARNING -',&
&   '  The 2nd degree equation has no positive root (choice=4).'
   status=2
!   call wrtout(06,*,'COLL')
   if(etotal_1<etotal_2)then
    write(16, '(a,a,/,a,/,a)' )"  ",&
&    ' findmin : COMMENT -',&
&    '  Will continue, since the new total energy is lower',&
&    '  than the old. Take a larger step in the same direction.'
!    call wrtout(06,*,'COLL')
    lambda_predict=1.5d0
   else
    write(16, '(a,a,/,a,/,a,/,a)' )"  ",&
&    ' findmin : COMMENT -',&
&    '  There is a problem, since the new total energy is larger',&
&    '  than the old (choice=4).',&
&    '  I take a point between the old and new, close to the old .'
!    call wrtout(06,*,'COLL')
    lambda_predict=0.25d0
   end if
!  Mimick a zero-gradient lambda, in order to avoid spurious
!  action of the inverse hessian (the next line would be a realistic estimation)
   dedv_predict=0.0d0
!  dedv_predict=dedv_2+lambda_predict*(dedv_1-dedv_2)
!  Uses the energies, and the gradient at lambda_2
   etotal_predict=etotal_2+dedv_2*lambda_predict&
&   +(etotal_1-etotal_2-dedv_2)*lambda_predict**2

  else

!  Here, there is an acceptable solution to the 2nd degree equation
   discr=sqrt(discr)
!  The root that gives the smallest ee corresponds to  -discr
!  This is the one to be used: one aims at modelling the
!  behaviour of the function as much as possible with the
!  lowest orders of the polynomial, not the quartic term.
   ee=(sum2-discr)*0.5d0
   dd=sum3-2.0d0*ee
   cc=sum1-dd-ee

!  My additional part
   if(ee==0.d0) then
   status=3
   return
   endif
!  END My additional part   



!  DEBUG
!  write(6,*)'aa,bb,cc,dd,ee',aa,bb,cc,dd,ee
!  ENDDEBUG

!  Now, must find the unique root of
!  $0 = bb + 2*cc * lambda + 3*dd * lambda^2 + 4*ee * lambda^3$
!  This root is unique because it was imposed that the second derivative
!  of the quartic polynomial is everywhere positive.
!  First, remove the quadratic term, by a shift of lambda
!  lambdap=lambda-lambda_shift
!  $0 = bbp + ccp * lambdap + eep * lambdap^3$
   eep=4.0d0*ee
   lambda_shift=-dd/(4.0d0*ee)
   ccp=2.0d0*cc-12.0d0*ee*lambda_shift**2
   bbp=bb+ccp*lambda_shift+eep*lambda_shift**3

!  DEBUG
!  write(6,*)'bbp,ccp,eep,lambda_shift',bbp,ccp,eep,lambda_shift
!  ENDDEBUG

!  The solution of a cubic polynomial equation is as follows :
   discr=(bbp/eep)**2+(4.0d0/27.0d0)*(ccp/eep)**3
   write(16,*)"Discr",discr,(4.0d0/27.0d0)*(ccp/eep)**3,(bbp/eep)**2
   if(discr.lt.0.d0) then
   status=3
   return
   endif
!  In the present case, discr will always be positive
   discr=sqrt(discr)
   uu3=0.5d0*(-bbp/eep+discr) ; uu=sign((abs(uu3))**(1.0d0/3.0d0),uu3)
   vv3=0.5d0*(-bbp/eep-discr) ; vv=sign((abs(vv3))**(1.0d0/3.0d0),vv3)
   lambda_predict=uu+vv
   write(16,*) "Shift",lambda_shift,lambda_predict,bbp,eep,discr
!  Restore the shift
   lambda_predict=lambda_predict+lambda_shift
   etotal_predict=aa+bb*lambda_predict+cc*lambda_predict**2+&
&   dd*lambda_predict**3+ee*lambda_predict**4
   dedv_predict=bb+2.0d0*cc*lambda_predict+3.0d0*dd*lambda_predict**2+&
&   4.0d0*ee*lambda_predict**3
   d2edv2_1=2*cc+6*dd*lambda_1+12*ee*lambda_1**2
   d2edv2_2=2*cc+6*dd*lambda_2+12*ee*lambda_2**2
   d2edv2_predict=2*cc+6*dd*lambda_predict+12*ee*lambda_predict**2

  end if

 else if(choice==1) then
 write(16, '(a,i3)' )'   line minimization, algorithm ',choice

! Use the derivative information to predict lambda
  d2edv2_mid=(dedv_1-dedv_2)/d_lambda
  lambda_predict=lambda_2-dedv_2/d2edv2_mid
  dedv_predict=dedv_2+(lambda_predict-lambda_2)*d2edv2_mid
  d2edv2_1=d2edv2_mid
  d2edv2_2=d2edv2_mid
  d2edv2_predict=d2edv2_mid
! also use the first energy to predict new energy
  etotal_predict=etotal_1+dedv_1*(lambda_predict-lambda_1)&
&  +0.5d0*d2edv2_1*(lambda_predict-lambda_1)**2
  etotal_2bis=etotal_1+dedv_1*(lambda_2-lambda_1)&
&  +0.5d0*d2edv2_1*(lambda_2-lambda_1)**2

  if(d2edv2_mid<0.0d0)then
   write(16, '(a,a,/,a,es18.10,a)' ) "  ",&
&   ' findmin : WARNING -',&
&   '  (scfcge) The second derivative is negative, equal to',d2edv2_mid        ,'.'
!   call wrtout(6,*,'COLL')
   status=1
  end if

 else if(choice==2) then

! Use energies and first derivative information
! etotal = aa + bb * lambda + cc * lambda**2
  dedv_mid2=(etotal_1-etotal_2)/d_lambda
  cc=(dedv_1-dedv_mid2)/d_lambda
  lambda_predict=lambda_1-0.5d0*dedv_1/cc
  d2edv2_1=2*cc
  d2edv2_2=d2edv2_1
  d2edv2_predict=d2edv2_1
  if(d2edv2_predict<0.0d0)then
   write(16, '(a,a,/,a,es18.10,a,/,a)' ) "  ",&
&   ' findmin : WARNING -',&
&   '  (scfcge) The second derivative is negative, equal to',d2edv2_predict,'.',&
&   '  (scfcge) => Pivoting                     '
!  call wrtout(6,*,'COLL')
   status=1
   if(etotal_2 < etotal_1)then
    lambda_predict=lambda_2-0.5d0*(lambda_1-lambda_2)
   else
    lambda_predict=lambda_1-0.5d0*(lambda_2-lambda_1)
   end if
  end if
  dedv_predict=dedv_1+(lambda_predict-lambda_1)*d2edv2_1
  dedv_2bis=dedv_1+(lambda_2-lambda_1)*d2edv2_1
  etotal_predict=etotal_1+dedv_1*(lambda_predict-lambda_1)&
&  +0.5d0*d2edv2_1*(lambda_predict-lambda_1)**2

 end if
 printvol=1
 if(choice==4)printvol=2
 if(printvol==2)then
  write(16, '(a,i3)' )'   line minimization, algorithm ',choice
! call wrtout(6,*,'COLL')
  write(16, '(a,a)' )'                        lambda      etotal ',&
&  '           dedv        d2edv2    '
! call wrtout(6,*,'COLL')
  write(16, '(a,es12.4,es18.10,2es12.4)' )&
&  '   old point         :',lambda_2,etotal_2,dedv_2,d2edv2_2
! call wrtout(6,*,'COLL')
  write(16, '(a,es12.4,es18.10,2es12.4)' )&
&  '   new point         :',lambda_1,etotal_1,dedv_1,d2edv2_1
! call wrtout(6,*,'COLL')
  write(16, '(a,es12.4,es18.10,2es12.4)' )&
&  '   predicted point   :',&
&  lambda_predict,etotal_predict,dedv_predict,d2edv2_predict
! call wrtout(6,*,'COLL')
  if(choice==1) then
   write(16, '(a,es10.4)' ) &
&   ' consistency check :    etotal_2 =',etotal_2bis
!  call wrtout(6,*,'COLL')
  end if
  if(choice==2) then
   write(16, '(a,es10.4)' ) &
&   ' consistency check :    dedv_2 =',dedv_2bis
!  call wrtout(6,*,'COLL')
  end if
  write(16, '(a)' ) ' '
!  call wrtout(6,*,'COLL')
 else if(printvol==1)then
  write(16, '(a,es12.4,a,es18.10)' ) &
&  '   findmin : lambda_predict ',lambda_predict,&
&  '   etotal_predict ',etotal_predict
!  call wrtout(6,*,'COLL')
 end if

end subroutine findmin


