subroutine slab_stress(flat,fix_z)
!This routine will eliminate all z-components of the first two cell forces,
!and all x-y-components of the last cell force 
implicit none
integer:: i,j
real(8):: flat(3,3),ekin1,ekin2
logical:: fix_z
ekin1=0.d0
do i=1,3
 do j=1,3
 ekin1=ekin1+flat(i,j)**2
 enddo
enddo

flat(3,1)=0.d0
flat(3,2)=0.d0
flat(1,3)=0.d0
flat(2,3)=0.d0
if(fix_z) flat(3,3)=0.d0

ekin2=0.d0
do i=1,3
 do j=1,3
 ekin2=ekin2+flat(i,j)**2
 enddo
enddo
write(*,'(a,3(1x,es20.10))') "# Befor/After eliminating restricted components: ",ekin1,ekin2,ekin1-ekin2
end subroutine
