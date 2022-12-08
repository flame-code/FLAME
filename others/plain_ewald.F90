!****************************************************************************************************
subroutine plain_ewald(atoms,en)
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_atoms), intent(inout):: atoms
    real(8), allocatable, intent(inout)::sfactor_norm(:,:,:) 
    integer, intent(inout):: kx, ky, kz
    !local variables
    real(8):: pi, q(1:atoms%nat)
    real(8):: k, k2, alpha, sum_en, en
    integer:: i, k_max
    pi=4.d0*atan(1.d0)
    sum_en=0.d0
    atoms%fat=0.d0
    do kx=1,k_max
        do ky=1,k_max
            do kz=1,k_max
                k2=kx**2+ky**2+kz**2
                k=sqrt(k2)
                call structur_factor_comput(atoms,kx,ky,kz,s1,s2,sfactor_norm)  
                sum_en=sum_en+(exp(-(k**2*a**2)/4.d0)/(k2))*sfactor_norm
                alpha=kx*atoms%rat(1,i)+ky*atoms%rat(2,i)+kz*atoms%rat(3,i)
                atoms%fat(1,i)=atoms%fat(1,i)+q(i)*kx*(exp(-((k**2)*(a**2))/4.d0)/(k2))*((sin(alpha)*s1)-(cos(alpha)*s2))          
                atoms%fat(2,i)=atoms%fat(2,i)+q(i)*ky*(exp(-((k**2)*(a**2))/4.d0)/(k2))*((sin(alpha)*s1)-(cos(alpha)*s2))          
                atoms%fat(3,i)=atoms%fat(3,i)+q(i)*kz*(exp(-((k**2)*(a**2))/4.d0)/(k2))*((sin(alpha)*s1)-(cos(alpha)*s2))          
            enddo
        enddo
    enddo
    en=sum_en*2.d0*pi/vol
    !electrostic force componants on particle i
    atoms%fat(1,i)= atoms%fat(1,i)*4.d0*pi/vol
    atoms%fat(2,i)= atoms%fat(2,i)*4.d0*pi/vol
    atoms%fat(3,i)= atoms%fat(3,i)*4.d0*pi/vol

end subroutine
!****************************************************************************************************
subroutine structur_factor_comput(atoms,kx,ky,kz,s1,s2,sfactor_norm)
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_atoms), intent(inout):: atoms
    real(8), allocatable, intent(inout)::sfactor_norm(:,:,:) 
    integer, intent(inout):: kx, ky, kz
    !local variables
    real(8), allocatable:: s1(:,:,:), s2(:,:,:)
    real(8):: pi, q(1:atoms%nat)
    integer:: i, k_max
    pi=4.d0*atan(1.d0)
    epot_long=0.d0
    s1=0.d0
    s2=0.d0
    do i=1,atoms%nat
        s1(kx,ky,kz)=s1(kx,ky,kz)+q(k)*cos((kx*atoms%rat(1,i)+ky*atoms%rat(2,i)+kz*atoms%rat(3,i))
        s2(kx,ky,kz)=s2(kx,ky,kz)+q(k)*sin(kx*atoms%rat(1,i)+ky*atoms%rat(2,i)+kz*atoms%rat(3,i))
    enddo
    sfactor_norm(kx,ky,kz)=s1(kx,ky,kz)**2+s2(kx,ky,kz)**2
end subroutine
!****************************************************************************************************
