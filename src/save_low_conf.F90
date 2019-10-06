subroutine save_low_conf(nat,npmin,npminx,ent_wpos,e_wpos,pos,latvec,spg,spgtol,fdos,elocmin,poslocmin,latlocmin)
!save configuration if it is among the lowest ones in energy
!  implicit real*8 (a-h,o-z)
  implicit none
  integer:: iat,nat, npmin, npminx, kmax, k 
  real(8):: e_wpos, ent_wpos, emax,spg,spgtol,fdos
  real(8):: elocmin(npminx,5)
  real(8):: pos(3,nat),latvec(3,3),poslocmin(3,nat,npminx),latlocmin(3,3,npminx)

  if (npmin.le.npminx) then
     kmax=npmin
     elocmin(kmax,1)=ent_wpos
     elocmin(kmax,2)=e_wpos
     elocmin(kmax,3)=spg
     elocmin(kmax,4)=spgtol
     elocmin(kmax,5)=fdos
     do iat=1,nat
        poslocmin(1,iat,kmax)=pos(1,iat)
        poslocmin(2,iat,kmax)=pos(2,iat)
        poslocmin(3,iat,kmax)=pos(3,iat)
     enddo
     latlocmin(:,:,kmax)=latvec(:,:)
  else
     ! find configuration kmax that is highest in energy
     emax=-1.d100
     do k=1,npminx
        if (elocmin(k,1).gt.emax) then
           emax=elocmin(k,1)
           kmax=k
        endif
     enddo
     if (ent_wpos.lt.elocmin(kmax,1)) then
        elocmin(kmax,1)=ent_wpos
        elocmin(kmax,2)=e_wpos
        elocmin(kmax,3)=spg
        elocmin(kmax,4)=spgtol
        elocmin(kmax,5)=fdos
        do iat=1,nat
           poslocmin(1,iat,kmax)=pos(1,iat)
           poslocmin(2,iat,kmax)=pos(2,iat)
           poslocmin(3,iat,kmax)=pos(3,iat)
        enddo
        latlocmin(:,:,kmax)=latvec(:,:)
     endif
  endif
  return
END SUBROUTINE save_low_conf

