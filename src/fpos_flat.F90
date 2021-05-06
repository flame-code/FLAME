subroutine fpos_flat(parini,pressure,fpos,flat,strten,fcart,latvec,md_type) 
!Computes the pure generalized forces on atom and cell (no contributions from velocities)
use mod_parini, only: typ_parini

implicit none
type(typ_parini), intent(in):: parini
integer:: iat,i,j,md_type
real(8),dimension(3,parini%nat):: fcart,fpos
real(8),dimension(3,3)  :: latvec,tmplat,pressure,a,velmat,sigma,flat,str_matrix
real(8),dimension(3,3)  :: term1,term2,term3,term4,term5,term5_1,term5_2,sigmatrans,f0inv
real(8):: amass(parini%nat),latmass,crossp(3),strten(6),vol,vpostmp(3),volvel,trace3,vol_1_3
           flat=0.d0
!Get volume
           a=latvec
           vol= a(1,1)*a(2,2)*a(3,3)-a(1,1)*a(2,3)*a(3,2)-a(1,2)*a(2,1)*a(3,3)+&
                a(1,2)*a(2,3)*a(3,1)+a(1,3)*a(2,1)*a(3,2)-a(1,3)*a(2,2)*a(3,1)
!Get full stress matrix
           str_matrix(1,1)=strten(1)
           str_matrix(2,2)=strten(2)
           str_matrix(3,3)=strten(3)
           str_matrix(1,2)=strten(6)
           str_matrix(2,1)=strten(6)
           str_matrix(1,3)=strten(5)
           str_matrix(3,1)=strten(5)
           str_matrix(2,3)=strten(4)
           str_matrix(3,2)=strten(4)
!Compute sigma
        call cross_product(latvec(:,2),latvec(:,3),crossp); sigma(:,1)=crossp
        call cross_product(latvec(:,3),latvec(:,1),crossp); sigma(:,2)=crossp
        call cross_product(latvec(:,1),latvec(:,2),crossp); sigma(:,3)=crossp
!Compute the acceleration of the cell
        flat=-str_matrix
!Here the pressure is applied
        flat=flat-pressure
if(md_type.ne.4) then
!Multiply with sigma from left
        flat=matmul(flat,sigma)
!Convert cartesian forces to reduced forces on atoms
        call invertmat(latvec,tmplat,3)
        do iat=1,parini%nat
          fpos(:,iat)=matmul(tmplat,fcart(:,iat))
        enddo
else
!Andersen
!Be careful: the positions are not the real positions but scaled with vol**(1/3)
          vol_1_3=vol**(1.d0/3.d0)
          do iat=1,parini%nat
              fpos(:,iat)=fcart(:,iat)/vol_1_3
          enddo
!Andersen MD
!Compute the hydrostatic pressure stuff
        flat(1,1)=(flat(1,1)+flat(2,2)+flat(3,3))/3.d0
endif
end subroutine
