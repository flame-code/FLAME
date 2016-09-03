!*****************************************************************************************
!Module for interfacing to the code openmx in order to provide forces and energy
module energyandforces
implicit none
real(8)::tstime,cell(3)
integer::icount
character(5)::sat(1000)
contains
!*****************************************************************************************
subroutine initialize_openmx(iproc,nproc,nat,rat)
    implicit none
    integer::iproc,nproc,nat,iat,ii
    real(8)::rat(3,1000)
    character(5)::ttsat
    !character(5)::sat_t(1000)
    call init_openmx(iproc,nproc,tstime,nat,rat,sat,cell)
    if(nat>1000) stop 'ERROR: nat>100, program stops in initialize_openmx'
    icount=0
    !sat(1:1000)=sat_t(1:1000)
    !do iat=1,nat
    !    ttsat=sat(iat)
    !    ii=len_trim(ttsat)
    !    sat(iat)=ttsat(1:ii)
    !enddo
    return
endsubroutine initialize_openmx
!*****************************************************************************************
subroutine calenergyforces(iproc,n,rat,fat,epot)
    implicit none
    integer::n,nat,iproc
    real*8 ::rat(3,n/3),fat(3,n/3),epot
    !character*10, intent(in):: caller          
    nat=n/3
    call openmx_energyandforces(iproc,nat,rat,epot,fat,icount)
    !write(*  ,'(a,i9,E,E13.3,5x)')  'iter,E[eV],frms[eV/A]:',iter,epot,sqrt(sum(fat(:,:)**2)/(3*nat))
    return
end subroutine calenergyforces
!*****************************************************************************************
subroutine finalize_openmx(iproc,nproc)   
    implicit none
    integer::iproc,nproc
    call final_openmx(iproc,nproc,tstime)
    return
endsubroutine finalize_openmx
endmodule energyandforces 
!*****************************************************************************************
