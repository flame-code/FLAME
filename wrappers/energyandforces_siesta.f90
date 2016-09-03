!*****************************************************************************************
module energyandforces
    implicit none
    integer::icount
    real(8)::cell(3)
    character(5)::sat(1000)
    contains
!*****************************************************************************************
subroutine calenergyforces(iproc,n,rat,fat,epot)
    use m_siesta_init
    !use m_siesta_analysis
    use m_siesta_move
    use m_siesta_end
    use m_siesta_forces
    use m_energies, only:Etot
    use m_forces, only:fa
    use siesta_geom, only:xa
    integer::iproc,n
    real(8)::rat(3,n/3),fat(3,n/3),epot
    xa(1:3,1:n/3)=rat(1:3,1:n/3)
    !write(*,*) 'IOnode ',IOnode
    call siesta_forces(icount)
    fat(1:3,1:n/3)=fa(1:3,1:n/3)
    epot=Etot
    icount=icount+1
end subroutine calenergyforces
end module energyandforces
!*****************************************************************************************
subroutine my_init_siesta(iproc,nproc,nat,natmax,rat)
    use m_siesta_init
    !use m_siesta_analysis
    use m_siesta_move
    use m_siesta_end
    use m_siesta_forces
    use siesta_geom, only:na_s,na_u,xa,ucell,scell,isa !,cisa
    !use atm_types, only:species
    use atmfuncs, only:labelfis
    use parallel, only:node,nodes
    !USE m_steps, only: inicoor, fincoor
    use energyandforces, only:icount,cell,sat
    implicit none
    integer::iproc,nproc,nat,natmax,iat
    real(8)::rat(3,natmax)
    call siesta_init()
    !stop
    nat=na_s
    rat(1:3,1:nat)=xa(1:3,1:nat)
    iproc=node
    nproc=nodes
    cell(1)=scell(1,1)
    cell(2)=scell(2,2)
    cell(3)=scell(3,3)
    do iat=1,nat
        sat(iat)=trim(labelfis(isa(iat)))
        !sat(iat)=trim(cisa(iat))
        !sat(iat)=species(labelfis(isa(iat)))%label
        if(iproc==0) write(*,*) 'SAT ',iat,trim(sat(iat))
    enddo
    icount=0
    if(iproc==0) write(*,*) 'na_s,na_u ',na_s,na_u
    if(iproc==0) write(*,*) 'node,nodes ',node,nodes
    if(iproc==0) write(*,*) 'nat ',nat
    if(iproc==0) write(*,*) 'ucell ',ucell(1,1:3)
    if(iproc==0) write(*,*) 'ucell ',ucell(2,1:3)
    if(iproc==0) write(*,*) 'ucell ',ucell(3,1:3)
    if(iproc==0) write(*,*) 'scell ',scell(1,1:3)
    if(iproc==0) write(*,*) 'scell ',scell(2,1:3)
    if(iproc==0) write(*,*) 'scell ',scell(3,1:3)
    write(*,*) 'iproc,nproc ',node,nodes
end subroutine my_init_siesta
!*****************************************************************************************
subroutine my_final_siesta
    use m_siesta_init
    !use m_siesta_analysis
    use m_siesta_move
    use m_siesta_end
    use m_siesta_forces
    !USE m_steps, only: inicoor, fincoor
    implicit none
    !call siesta_analysis( relaxd )
    call siesta_end()
end subroutine my_final_siesta
!*****************************************************************************************
