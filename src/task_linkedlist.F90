subroutine  linkedlist_test(parini)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, type_pairs, typ_file_info, atom_deallocate_old
    use mod_atoms, only: atom_allocate, atom_allocate_old, update_ratp
    use mod_const, only: bohr2ang
    use mod_linked_lists, only: typ_linked_lists, typ_pia_arr
    use mod_acf, only: acf_read
    implicit none
    type(typ_linked_lists):: linked_lists
    type(typ_parini), intent(in):: parini
    type(typ_atoms):: atoms 
    type(typ_file_info):: file_info
    type(type_pairs),allocatable:: posat1st(:) ,posat1st_list(:)
    type(typ_pia_arr):: pia_arr
    integer:: iat, jat, nx, ny,nz, tt, niat,kat,t1,t2,nx2,ny2,nz2
    integer::limnx, limny ,limnz ,iba,ibr
    integer, allocatable::  nim(:),nim_list(:)
    real(8):: r1, rcut, dx, dy, dz ,xiat, yiat, ziat 
    real(8):: r2, dx2, dy2, dz2 
    real(8):: r3, dx3, dy3, dz3 ,rnd
    real(8):: latvec_x(3,3),rec_nec(3),dxyz(3),r1red(3),r2red(3)
    real(8),allocatable:: rat_int(:,:)
    real(8):: dd , diff(3)
    integer ::conf ,maxconf,nbc
    integer :: i,j,k,nec(3),ii
    logical :: yes
    character(2):: numb
    !**************************************************************
!    call random_seed()
    maxconf=20
    if (parini%posinp_misc== 'known') maxconf=1
do conf=1,maxconf
    if (parini%posinp_misc== 'known') then
        call acf_read(parini,'posinp.acf',1,atoms=atoms)
        allocate(posat1st(atoms%nat), nim(atoms%nat))
        allocate(posat1st_list(atoms%nat), nim_list(atoms%nat))
        rcut=12.00d0!bohr2ang
        do iat= 1,atoms%nat
            allocate(posat1st(iat)%posat2nd(3,9*atoms%nat))
            allocate(posat1st_list(iat)%posat2nd(3,9*atoms%nat))
        enddo
        goto 11
    endif

    if (conf<10)then
        numb= char(48)//char(conf+48)
    else if (conf<100)then
        numb=char(conf/10+48)//char(mod(conf,10)+48)
    endif
    atoms%cellvec(1,2)=0.d0
    atoms%cellvec(1,3)=0.d0
    atoms%cellvec(2,3)=0.d0
    atoms%boundcond = parini%boundcond_misc

    rcut=12.00d0!/bohr2ang
    call random_number(rnd)
    atoms%cellvec(1,1)=6.d0+rnd*30
    call random_number(rnd)
    atoms%cellvec(2,2)=6.d0+rnd*30
    call random_number(rnd)
    atoms%cellvec(3,3)=6.d0+rnd*30
    call random_number(rnd)
    atoms%nat=4+40*rnd
    !allocate(atoms%rat(3,atoms%nat),atoms%sat(atoms%nat))
    call atom_allocate_old(atoms,atoms%nat,0,0)
    allocate(posat1st(atoms%nat), nim(atoms%nat))
    allocate(posat1st_list(atoms%nat), nim_list(atoms%nat))
    do iat= 1,atoms%nat
        allocate(posat1st(iat)%posat2nd(3,9*atoms%nat))
        allocate(posat1st_list(iat)%posat2nd(3,9*atoms%nat))
    enddo
    atoms%sat='Si'
    atoms%bemoved=.true.

    call genrandomconf(atoms,numb,conf)

    call acf_read(parini,'posinp_'//char(49)//numb//'.acf',1,atoms=atoms)
11  call n_rep_dim_alborz(atoms%cellvec,2.d0*rcut,nec(1),nec(2),nec(3))
    allocate(rat_int(3,atoms%nat))

    write(*,*) "nec1,nec2,nec3:",nec(1),nec(2),nec(3)
    nbc =3
    if (trim(atoms%boundcond)=='slab') then
        nbc=2
        nec(3)=1
    elseif (trim(atoms%boundcond)=='wire') then
        nbc=1
        nec(3)=1
        nec(2)=1
    elseif (trim(atoms%boundcond)=='free') then
        nbc=0
        nec(3)=1
        nec(2)=1
        nec(1)=1
    endif
    !========================================== 
    !Expand cell
    do i=1,3
        latvec_x(:,i)=real(nec(i),8)*atoms%cellvec(:,i)
        !Adjust reduced coordinates
        rec_nec(i)=1.d0/real(nec(i),8)
    enddo
    !Get the rest

    call update_ratp(atoms)
    call rxyz_cart2int_alborz(atoms%nat,atoms%cellvec,atoms%ratp,rat_int)
call system_clock(t1)
    nim=0
    do iat= 1,atoms%nat
        r1red(:)=rat_int(:,iat)*rec_nec
        do i=0,nec(1)-1
        do j=0,nec(2)-1
        do k=0,nec(3)-1
        do jat=1,atoms%nat
            r2red(:)=rat_int(:,jat)*rec_nec
            r2red(1)=r2red(1)+real(i,8)*rec_nec(1)
            r2red(2)=r2red(2)+real(j,8)*rec_nec(2)
            r2red(3)=r2red(3)+real(k,8)*rec_nec(3)
!            call pbc_distance1(latvec_x,r1red,r2red,dd,dxyz)

            diff(1:3)=r2red(1:3)-r1red(1:3)

            do ii=1,nbc
                if(.not. diff(ii)>-0.5d0) then
                    diff(ii)=diff(ii)+1.d0
                elseif(diff(ii)>0.5d0) then
                    diff(ii)=diff(ii)-1.d0
                endif
            enddo
            dxyz(1)=latvec_x(1,1)*diff(1)+latvec_x(1,2)*diff(2)+latvec_x(1,3)*diff(3)
            dxyz(2)=latvec_x(2,1)*diff(1)+latvec_x(2,2)*diff(2)+latvec_x(2,3)*diff(3)
            dxyz(3)=latvec_x(3,1)*diff(1)+latvec_x(3,2)*diff(2)+latvec_x(3,3)*diff(3)
            dd=dxyz(1)*dxyz(1)+dxyz(2)*dxyz(2)+dxyz(3)*dxyz(3)
                r1=sqrt(dd)
                if (r1>0.d0 .and. r1<=rcut) then
                !write(*,*)iat,jat,r1,rcut
                    write(2000+conf,*)iat,jat
                    !nim(iat)=nim(iat)+1
                    !posat1st(iat)%posat2nd(1,nim(iat))=jat
                    !posat1st(iat)%posat2nd(2,nim(iat))=r1
                endif
        enddo
        enddo
        enddo
        enddo
    enddo

    deallocate(rat_int)
call system_clock(t2)
    write(*,*) "time for ordinary part = " ,t2-t1
   ! do iat= 1,atoms%nat
   !     do niat=1,nim(iat)
   !         write(2000+conf,*)iat,posat1st(iat)%posat2nd(1,niat)
   !     enddo
   ! enddo
    write(*,*)"============================================================================"
call system_clock(t1)
    linked_lists%triplex=.true.
    linked_lists%rcut=rcut
    call call_linkedlist(parini,atoms,.true.,linked_lists,pia_arr)
    do ibr= 1,linked_lists%maxbound_rad
        write(2200+conf,*)linked_lists%bound_rad(1,ibr),linked_lists%bound_rad(2,ibr)
        !if (atoms%sat(linked_lists%bound_rad(1,ibr))=='O' .and. atoms%sat(linked_lists%bound_rad(2,ibr))=='O' )then
        !    write(2200+conf,'(2i6,6es25.15)')linked_lists%bound_rad(1,ibr),linked_lists%bound_rad(2,ibr),(atoms%rat(:,linked_lists%bound_rad(1,ibr))+atoms%rat(:,linked_lists%bound_rad(2,ibr)))/2*bohr2ang
        !endif
    enddo
    
    do iba= 1,linked_lists%maxbound_ang
        iat = linked_lists%bound_rad(1,linked_lists%bound_ang(1,iba))
        jat = linked_lists%bound_rad(2,linked_lists%bound_ang(1,iba))
        kat = linked_lists%bound_rad(2,linked_lists%bound_ang(2,iba))
        write(3200+conf,'(3i6)')iat,min(jat,kat),max(jat,kat)
        !call sort2(iat,jat,kat,conf,3200)
    enddo
    deallocate(linked_lists%bound_rad)
    deallocate(linked_lists%bound_ang)
    deallocate(linked_lists%prime_bound)
    deallocate(pia_arr%pia)
    call system_clock(t2)
    write(*,*) "time for linked list2 part = " ,t2-t1
!******************************************************************************************************
    do iat= 1,atoms%nat
       deallocate(posat1st(iat)%posat2nd)
       deallocate(posat1st_list(iat)%posat2nd)
    enddo
    deallocate(posat1st, nim)
    deallocate(posat1st_list, nim_list)
    call atom_deallocate_old(atoms)
enddo
write(*,*) "*************************************  end   ************************************"    
end subroutine linkedlist_test
!**************************************************************************************************************
subroutine callinkedlist(parini,atoms,rcut,posat1st,nim,conf)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, type_pairs, update_ratp
    use mod_const, only: bohr2ang
    use mod_linked_lists, only: typ_linked_lists
    implicit none
    type(typ_parini):: parini
    type(typ_atoms):: atoms 
    type(typ_linked_lists):: linked_lists
    type(type_pairs):: posat1st(atoms%nat) 
    integer::nimat ,iat_maincell, jat_maincell
    integer:: istat, nim(atoms%nat)
    integer:: iat, jat, niat,kat,kkz,conf
    integer:: ix,iy,iz,jx,jy,jz,kx,ky,kz,dkjz,dkiy ,tt
    integer:: llx,mmx,lly,mmy ,jp,jl,kp,kl,ip,il ,kpt,jpt
    real(8):: sclinv,cell(3) ,rcut
    real(8):: r1, dx, dy, dz ,xiat, yiat, ziat
    real(8):: r2, dx2, dy2, dz2 
    real(8):: r3, dx3, dy3, dz3 ,hxinv,rat_i(3)
    linked_lists%rcut=rcut
    linked_lists%triplex=.true.
    call linkedlists_init(parini,atoms,cell,linked_lists)
    hxinv=real(linked_lists%mx,8)/cell(1)
    !-------------------------------------------------------
    nim=0
    call update_ratp(linked_lists%typ_atoms)
    include 'act1_cell_linkedlist.inc'
    do  iat=ip,il
        xiat=linked_lists%ratp(1,iat)
        yiat=linked_lists%ratp(2,iat)
        ziat=linked_lists%ratp(3,iat)
        jpt=linked_lists%prime(ix+linked_lists%limnbx(1,jy-iy,jz-iz),jy,jz)
        jp=(iat-ip+1)*((isign(1,ip-jpt)+1)/2)+jpt
        jl=linked_lists%last(ix+linked_lists%limnbx(2,jy-iy,jz-iz),jy,jz)
        do  jat=jp,jl
            call rxyz_cart2int_alborz(1,atoms%cellvec,linked_lists%ratp(:,jat),rat_i)
            !jx=floor(linked_lists%ratp(1,jat)*hxinv)+1
            jx=floor(rat_i(1)*linked_lists%mx)+1
            dx=xiat-linked_lists%ratp(1,jat)
            dy=yiat-linked_lists%ratp(2,jat)
            dz=ziat-linked_lists%ratp(3,jat)
            r1=(dx**2+dy**2+dz**2)
            if (r1<linked_lists%rcut**2) then
                if (linked_lists%maincell(iat)+linked_lists%maincell(jat) >-1 ) then
                    iat_maincell=linked_lists%perm(iat)
                    jat_maincell=linked_lists%perm(jat)
                    if (iat_maincell /= jat_maincell) then
                        nim(iat_maincell)=nim(iat_maincell)+1
                        nim(jat_maincell)=nim(jat_maincell)+1
                    else
                        nim(iat_maincell)=nim(iat_maincell)+1
                    endif
                    posat1st(iat_maincell)%posat2nd(1,nim(iat_maincell))=jat_maincell
                    posat1st(jat_maincell)%posat2nd(1,nim(jat_maincell))=iat_maincell
                !    if (linked_lists%maincell(iat)==-1) then
                !        posat1st(iat_maincell)%posat2nd(3,nim(iat_maincell))=-floor(real(iy-1)/linked_lists%my)
                !        posat1st(jat_maincell)%posat2nd(3,nim(jat_maincell))=floor(real(iy-1)/linked_lists%my)
                !        posat1st(iat_maincell)%posat2nd(2,nim(iat_maincell))=-(1-isign(1,linked_lists%last(linked_lists%mx,iy,iz)-iat))/2
                !        posat1st(jat_maincell)%posat2nd(2,nim(jat_maincell))=(1-isign(1,linked_lists%last(linked_lists%mx,iy,iz)-iat))/2
                !    elseif (linked_lists%maincell(jat)==-1) then
                !        posat1st(iat_maincell)%posat2nd(3,nim(iat_maincell))=floor(real(jy-1)/linked_lists%my)
                !        posat1st(jat_maincell)%posat2nd(3,nim(jat_maincell))=-floor(real(jy-1)/linked_lists%my)
                !        posat1st(iat_maincell)%posat2nd(2,nim(iat_maincell))=(1-isign(1,linked_lists%last(linked_lists%mx,jy,jz)-jat))/2
                !        posat1st(jat_maincell)%posat2nd(2,nim(jat_maincell))=-(1-isign(1,linked_lists%last(linked_lists%mx,jy,jz)-jat))/2
                !    else 
                !        posat1st(iat_maincell)%posat2nd(3,nim(iat_maincell))=0
                !        posat1st(jat_maincell)%posat2nd(3,nim(jat_maincell))=0
                !        posat1st(iat_maincell)%posat2nd(2,nim(iat_maincell))=0
                !        posat1st(jat_maincell)%posat2nd(2,nim(jat_maincell))=0
                !    endif
                endif
    !++++++++++++++++++++++++++  Calculate triplex for angular in half of sphere +++++++++++++++++++++++++++++++++++++++            
                do kz=iz,min(linked_lists%mz+linked_lists%mlimnb3,iz+linked_lists%mlimnb)
                do ky=iy+linked_lists%limnby(1,kz-iz),iy+linked_lists%limnby(2,kz-iz)
                    kpt=linked_lists%prime(ix+linked_lists%limnbx(1,ky-iy,kz-iz),ky,kz)
                    kp=(jat-jp+1)*((isign(1,jp-kpt)+1)/2)+kpt
                    kl=linked_lists%last(ix+linked_lists%limnbx(2,ky-iy,kz-iz),ky,kz)
                    do kat=kp,kl
                        if (kat<=iat .or. kat<=jat) cycle
                        if ((linked_lists%maincell(iat)+linked_lists%maincell(jat)+linked_lists%maincell(kat)) <-2 )cycle
                        dx3=linked_lists%ratp(1,kat)-linked_lists%ratp(1,jat)
                        dy3=linked_lists%ratp(2,kat)-linked_lists%ratp(2,jat)
                        dz3=linked_lists%ratp(3,kat)-linked_lists%ratp(3,jat)
                        r3=(dx3**2+dy3**2+dz3**2)
                        if (r3<linked_lists%rcut**2) cycle
                        dx2=xiat-linked_lists%ratp(1,kat)
                        dy2=yiat-linked_lists%ratp(2,kat)
                        dz2=ziat-linked_lists%ratp(3,kat)
                        r2=(dx2**2+dy2**2+dz2**2)
                        if (r2<linked_lists%rcut**2 ) then
                            call sort_alborz(linked_lists%perm(iat),linked_lists%perm(jat),linked_lists%perm(kat),conf)
                        ! call sort(iat,jat,kat,conf)
                         !write(3100+conf,*)linked_lists%maincell(iat),linked_lists%maincell(jat),linked_lists%maincell(kat)
                        endif
                    enddo
                enddo
                enddo
            endif
    !++++++++++++++++++++++++++  Calculate triplex for angular for full sphere +++++++++++++++++++++++++++++++++++++++            
            do kz=max(1,jz-linked_lists%mlimnb),min(linked_lists%mz+linked_lists%mlimnb3,jz+linked_lists%mlimnb)
            dkjz=abs(kz-jz)      
            do ky=jy-linked_lists%limnby(2,dkjz),jy+linked_lists%limnby(2,dkjz)
                 kpt=linked_lists%prime(jx-linked_lists%limnbx(2,abs(ky-jy),dkjz),ky,kz)
                 kp=kpt
                 kl=linked_lists%last(jx+linked_lists%limnbx(2,abs(ky-jy),dkjz),ky,kz)
                 do kat=kp,kl
                 if (iat>=kat .or. jat==kat) cycle
                 if ((linked_lists%maincell(iat)+linked_lists%maincell(jat)+linked_lists%maincell(kat)) <-2 )cycle
                 dx2=xiat-linked_lists%ratp(1,kat)
                 dy2=yiat-linked_lists%ratp(2,kat)
                 dz2=ziat-linked_lists%ratp(3,kat)
                 r2=(dx2**2+dy2**2+dz2**2)
                 dx3=linked_lists%ratp(1,kat)-linked_lists%ratp(1,jat)
                 dy3=linked_lists%ratp(2,kat)-linked_lists%ratp(2,jat)
                 dz3=linked_lists%ratp(3,kat)-linked_lists%ratp(3,jat)
                 r3=(dx3**2+dy3**2+dz3**2)
                 if (r1<linked_lists%rcut**2 .and. r2<linked_lists%rcut**2 .and. r3<linked_lists%rcut**2)  then
                     if (jat<kat) then
                     call sort_alborz(linked_lists%perm(iat),linked_lists%perm(jat),linked_lists%perm(kat),conf)
                     !    call sort(iat,jat,kat,conf)
                        ! write(3100+conf,*)linked_lists%maincell(iat),linked_lists%maincell(jat),linked_lists%maincell(kat)
                     endif
                 elseif (r1<linked_lists%rcut**2.and. r3<linked_lists%rcut**2) then
                     call sort_alborz(linked_lists%perm(iat),linked_lists%perm(jat),linked_lists%perm(kat),conf)
                     !    call sort(iat,jat,kat,conf)
                        ! write(3100+conf,*)linked_lists%maincell(iat),linked_lists%maincell(jat),linked_lists%maincell(kat)
                 endif
                 enddo
            enddo
            enddo
            !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        enddo
    enddo
    include 'act2_cell_linkedlist.inc'
    call linkedlists_final(linked_lists)
end subroutine callinkedlist 
!**************************************************************************************************************
subroutine sort_alborz(i ,j ,k,conf)
    implicit none
    integer :: i, j, k
    integer ::conf
    if (i<=j) then
        if (j<=k) then
            write(3100+conf,'(3i6)')i,j,k
        else
            if(i<k) then
                write(3100+conf,'(3i6)')i,k,j
            else
                write(3100+conf,'(3i6)')k,i,j
            endif
        endif
    else
        if (j>=k) then
            write(3100+conf,'(3i6)')k,j,i
        elseif (i<k) then
            write(3100+conf,'(3i6)')j,i,k
        else
            write(3100+conf,'(3i6)')j,k,i
        endif
    endif
end subroutine sort_alborz
!**************************************************************************************************************
subroutine sort2_alborz(i ,j ,k,conf,num)
    implicit none
    integer :: i, j, k
    integer ::conf,num
    if (i<=j) then
        if (j<=k) then
            write(num+conf,'(3i6)')i,j,k
        else
            if(i<k) then
                write(num+conf,'(3i6)')i,k,j
            else
                write(num+conf,'(3i6)')k,i,j
            endif
        endif
    else
        if (j>=k) then
            write(num+conf,'(3i6)')k,j,i
        elseif (i<k) then
            write(num+conf,'(3i6)')j,i,k
        else
            write(num+conf,'(3i6)')j,k,i
        endif
    endif
end subroutine sort2_alborz
!*******************************************************************************************************
subroutine genrandomconf(atoms,numb,conf)
    use mod_atoms, only: typ_atoms
    use mod_atoms, only: typ_atoms, typ_file_info, update_rat
    use mod_acf, only: acf_write
    implicit none
    real(8)::amargin ,amargin_xrel ,amargin_yrel ,amargin_zrel
    real(8):: ranxyz(3)
    integer ::mat,conf, itry
    character(2):: numb
    type(typ_atoms):: atoms 
    type(typ_file_info):: file_info
    amargin=0
    amargin_xrel=amargin/atoms%cellvec(1,1)
    amargin_yrel=amargin/atoms%cellvec(2,2)
    amargin_zrel=amargin/atoms%cellvec(3,3)
    mat=0
    itry=0
    do
        mat=mat+1
        if(mat>atoms%nat) exit
        call random_number(ranxyz)
        ranxyz(1)=ranxyz(1)*(1.d0-2.d0*amargin_xrel)+amargin_xrel
        ranxyz(2)=ranxyz(2)*(1.d0-2.d0*amargin_yrel)+amargin_yrel
        ranxyz(3)=ranxyz(3)*(1.d0-2.d0*amargin_zrel)+amargin_zrel
        atoms%ratp(1,mat)=ranxyz(1)*atoms%cellvec(1,1)
        atoms%ratp(2,mat)=ranxyz(2)*atoms%cellvec(2,2)
        atoms%ratp(3,mat)=ranxyz(3)*atoms%cellvec(3,3)
    enddo
    call update_rat(atoms,upall=.true.)
    

    file_info%filename_positions='posinp_'//char(49)//numb//'.acf'
    file_info%file_position='new'
    call acf_write(file_info,atoms=atoms,strkey='posinp')
end subroutine genrandomconf
!*****************************************************************************************
