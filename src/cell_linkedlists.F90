!| natim1        |         natim5                      | natim1   :   natim2 |                                    |               |  
!|=====================================================|========  : =========|====================================|===============|
!| NATIM3/NATIM4 |                                     | natim3   :   natim8 |               NATIM9               |     NATIM8    |
!|---------------+----------------------------+--------|--------  : ---------|--------+---------------------------+---------------|
!|               |                            |        |          :          |                                    |               |        
!|     NATIM4    |                            |        | natim4   :          |                                    |     NATIM6    |  
!|               |                            |        |          :   natim6 |                NATIM7              |               |        
!|---------------+----------------------------+--------|--------  :          |                                    +---------------|
!|     NATIM1    |         NATIM5                      | natim1   :          |                                    | NATIM2/NATIM6 | 
!|=====================================================|========  : =========|====================================|===============|
!|               |                            |        | natim3   :   natim8 |                natim9              |     natim8
!where names with upper case letters indicate number of atoms of
!each segments in the main cell and those with lower case
!letters indicate the number of atoms in periodic images.
!*****************************************************************************************
subroutine linkedlists_init(parini,atoms,cell,linked_lists)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, atom_allocate_old
    use mod_electrostatics, only: typ_linked_lists
    use mod_const, only: bohr2ang
    use yaml_output
    implicit none
    type(typ_parini), intent(in):: parini 
    type(typ_atoms), intent(in):: atoms
    real(8), intent(out):: cell(3)
    type(typ_linked_lists), intent(inout):: linked_lists
    !local variables
    real(8):: sclinv, vol,tmp(3), nrm_tmp
    integer:: istat 
    linked_lists%avgnndis=2.8d0/bohr2ang
    associate(mlimnb=>linked_lists%mlimnb)
    !linked_lists%avgnndis=2.82d0   
    !linked_lists%avgnndis=(cell(1)*cell(2)*cell(3)/atoms%nat)**(1.d0/3.d0)
    mlimnb=ceiling(linked_lists%rcut/linked_lists%avgnndis)
    mlimnb=max(2,2*mlimnb/3)
    linked_lists%scl=linked_lists%rcut/real(mlimnb,8) !length of cells in each direction.
    call cell_vol(atoms%nat,atoms%cellvec,vol) 
    vol=abs(vol)*atoms%nat
    !write(*,*)atoms%cellvec(:,1)
    call cross_product_alborz(atoms%cellvec(1:3,1),atoms%cellvec(1:3,2),tmp)
    nrm_tmp=sqrt(dot_product(tmp,tmp))
    cell(3)=vol/nrm_tmp
    call cross_product_alborz(atoms%cellvec(1:3,2),atoms%cellvec(1:3,3),tmp)
    nrm_tmp=sqrt(dot_product(tmp,tmp))
    cell(1)=vol/nrm_tmp
    call cross_product_alborz(atoms%cellvec(1:3,1),atoms%cellvec(1:3,3),tmp)
    nrm_tmp=sqrt(dot_product(tmp,tmp))
    cell(2)=vol/nrm_tmp
    if(parini%iverbose>1) then
        call yaml_mapping_open('linked list info-1') !,flow=.true.)
        call yaml_map('nat',atoms%nat)
        call yaml_map('scl',linked_lists%scl)
        call yaml_map('mlimnb',mlimnb)
        call yaml_map('cell',(/cell(1),cell(2),cell(3)/))
        call yaml_mapping_close()
        !write(*,*) 'nat,scl,mlimnb',atoms%nat,linked_lists%scl,mlimnb
        !write(*,*) 'cell  ', cell(1),cell(2),cell(3)
    endif
    sclinv=1.d0/linked_lists%scl
    linked_lists%mx=floor(cell(1)*sclinv) !calculating number of cells in x direction.
    linked_lists%my=floor(cell(2)*sclinv) !calculating number of cells in y direction.
    linked_lists%mz=floor(cell(3)*sclinv) !calculating number of cells in z direction.
    if(linked_lists%mx<1 .or. linked_lists%my<1 .or. linked_lists%mz<1) then
        linked_lists%scl=min(linked_lists%scl,cell(1),cell(2),cell(3))
        mlimnb=ceiling(linked_lists%rcut/linked_lists%scl)
        sclinv=1.d0/linked_lists%scl
        linked_lists%mx=max(1,floor(cell(1)*sclinv)) !calculating number of cells in x direction.
        linked_lists%my=max(1,floor(cell(2)*sclinv)) !calculating number of cells in y direction.
        linked_lists%mz=max(1,floor(cell(3)*sclinv)) !calculating number of cells in z direction.
    endif
    mlimnb=ceiling(linked_lists%rcut/linked_lists%scl)
    allocate(linked_lists%limnbx(2,-mlimnb:mlimnb,0:mlimnb),stat=istat)
    if(istat/=0) stop 'ERROR: failure allocating array limnbx.'
    allocate(linked_lists%limnby(2,0:mlimnb),stat=istat)
    if(istat/=0) stop 'ERROR: failure allocating array limnby.'
    end associate
!    if(linked_lists%triplex==.true.) then
!        linked_lists%next=2
!    else
!        linked_lists%next=1
!    endif
    !write(*,*)linked_lists%next
    if(trim(atoms%boundcond)=='bulk') then
        linked_lists%mlimnb1=linked_lists%mlimnb
        linked_lists%mlimnb2=linked_lists%mlimnb
        linked_lists%mlimnb3=linked_lists%mlimnb
    elseif(trim(atoms%boundcond)=='slab') then
        linked_lists%mlimnb1=linked_lists%mlimnb
        linked_lists%mlimnb2=linked_lists%mlimnb
        linked_lists%mlimnb3=0
    elseif(trim(atoms%boundcond)=='wire') then
        linked_lists%mlimnb1=linked_lists%mlimnb
        linked_lists%mlimnb2=0
        linked_lists%mlimnb3=0
    elseif(trim(atoms%boundcond)=='free') then
        linked_lists%mlimnb1=0
        linked_lists%mlimnb2=0
        linked_lists%mlimnb3=0
    else
        write(*,'(2a)') 'ERROR: unknown BC in linkedlists_init ',trim(atoms%boundcond)
        stop
    endif
    call determine_sclimitsphere(linked_lists)
!    if(linked_lists%mx<linked_lists%mlimnb1/2) stop 'cell is the same size as supercell'
!    if(linked_lists%my<linked_lists%mlimnb2/2) stop 'cell is the same size as supercell'
!    if(linked_lists%mz<linked_lists%mlimnb3/2) stop 'cell is the same size as supercell'
    associate(mx=>linked_lists%mx,my=>linked_lists%my,mz=>linked_lists%mz)
    if(parini%iverbose>1) then
        call yaml_mapping_open('linked list info-2',flow=.true.)
        call yaml_map('mx,my,mz',(/mx,my,mz/))
        call yaml_map('total number of subcells',mx*my*mz)
        call yaml_mapping_close()
        !write(*,*) 'mx,my,mz,total number of cells',mx,my,mz,mx*my*mz
    endif
    allocate(linked_lists%head(mx,my,mz),stat=istat)
    if(istat/=0) stop 'ERROR: failure allocating linked_lists%head'
    allocate(linked_lists%list(atoms%nat),stat=istat)
    if(istat/=0) stop 'ERROR: failure allocating linked_lists%list'
    associate(n=>linked_lists%mlimnb)
    if (trim(atoms%boundcond)=='bulk') then
        allocate(linked_lists%prime(1-(n+linked_lists%mlimnb1):mx+(n+linked_lists%mlimnb1),&
                 1-(n+linked_lists%mlimnb2):my+(n+linked_lists%mlimnb2),1:mz+linked_lists%mlimnb3),stat=istat)
    else
        allocate(linked_lists%prime(1-n:mx+(n+linked_lists%mlimnb1),1-(n+linked_lists%mlimnb2):my+(n+linked_lists%mlimnb2),mz),stat=istat)
    endif
    if (trim(atoms%boundcond)=='bulk') then
        allocate(linked_lists%last(1-(n+linked_lists%mlimnb1):mx+(n+linked_lists%mlimnb1),&
                 1-(n+linked_lists%mlimnb2):my+(n+linked_lists%mlimnb2),1:mz+linked_lists%mlimnb3),stat=istat)
    else
        allocate(linked_lists%last(1-n:mx+(n+linked_lists%mlimnb1),1-(n+linked_lists%mlimnb2):my+(n+linked_lists%mlimnb2),mz),stat=istat)
    endif
    end associate
    end associate
    call make_list_new(parini,atoms,linked_lists,cell)
    call atom_allocate_old(linked_lists%typ_atoms,linked_lists%natim,linked_lists%natim,0)
    allocate(linked_lists%maincell(1:linked_lists%natim))
    allocate(linked_lists%perm(1:linked_lists%natim))
    call prepprimelast(atoms,linked_lists)
end subroutine linkedlists_init
!*****************************************************************************************
subroutine linkedlists_final(linked_lists)
    use mod_atoms, only: typ_atoms, atom_deallocate_old
    use mod_electrostatics, only: typ_linked_lists
    implicit none
    type(typ_linked_lists), intent(inout):: linked_lists
    !local variables
    integer:: istat
    deallocate(linked_lists%head,stat=istat)
    if(istat/=0) stop 'ERROR: failure allocating linked_lists%head'
    deallocate(linked_lists%list,stat=istat)
    if(istat/=0) stop 'ERROR: failure allocating linked_lists%list'
    deallocate(linked_lists%prime,stat=istat)
    if(istat/=0) stop 'ERROR: failure allocating linked_lists%prime'
    deallocate(linked_lists%last,stat=istat)
    if(istat/=0) stop 'ERROR: failure allocating linked_lists%last'
    deallocate(linked_lists%limnbx,stat=istat)
    if(istat/=0) then
        write(*,*) 'ERROR: failure deallocating array limnbx.'
    endif
    deallocate(linked_lists%limnby,stat=istat)
    if(istat/=0) then
        write(*,*) 'ERROR: failure deallocating array limnby.'
    endif
    deallocate(linked_lists%perm,stat=istat)

    call atom_deallocate_old(linked_lists%typ_atoms)
    deallocate(linked_lists%maincell)
end subroutine linkedlists_final
!*****************************************************************************************
subroutine prepprimelast(atoms,linked_lists)
    use mod_atoms, only: typ_atoms, get_rat, update_rat
    use mod_electrostatics, only: typ_linked_lists
    implicit none
    type(typ_atoms), intent(in):: atoms
    type(typ_linked_lists), intent(inout):: linked_lists
    !local variables
    integer:: iat, jat, ix, iy, iz, natpopnatim, ishiftx, ishifty, ishift_l, ishiftz
    integer:: nx, ny, nz
    integer:: mx, my, mz, mlimnb, mlimnb1, mlimnb2,bulkbc, ishift_r ,mlimnb3
    real(8), allocatable:: rat_int(:,:)
    real(8), allocatable:: rat(:,:)
    real(8):: rcart(3), rfrac(3)
    allocate(rat_int(1:3,1:atoms%nat))
    allocate(rat(1:3,1:atoms%nat))
    mx=linked_lists%mx
    my=linked_lists%my
    mz=linked_lists%mz
    mlimnb=linked_lists%mlimnb
    mlimnb1=linked_lists%mlimnb1
    mlimnb2=linked_lists%mlimnb2
    mlimnb3=linked_lists%mlimnb3
    bulkbc=0
    if (trim(atoms%boundcond)=='bulk') then
        bulkbc=1
    endif
    natpopnatim=linked_lists%natim+1
    !write(*,*)natpopnatim
    linked_lists%prime=natpopnatim
    !initialize to -1 means the atom is not in main cell and
    !later corresponding atoms in main cell will be +1.
    linked_lists%maincell=-1
    linked_lists%last=0
    call get_rat(atoms,rat)
    call rxyz_cart2int_alborz(atoms%nat,atoms%cellvec,rat,rat_int)
    !we are doing the following since in mutiple application of h and h^0-1, i.e.,
    !in going from reduced to Cartesian and vice versa, there can be rounding
    !which can put atoms outside cell, e.g. 0 can be -2.E-17
    if(trim(atoms%boundcond)=='bulk') then
        do iat=1,atoms%nat
            rat_int(1,iat)=modulo(modulo(rat_int(1,iat),1.d0),1.d0)
            rat_int(2,iat)=modulo(modulo(rat_int(2,iat),1.d0),1.d0)
            rat_int(3,iat)=modulo(modulo(rat_int(3,iat),1.d0),1.d0)
        enddo
    elseif(trim(atoms%boundcond)=='slab') then
        do iat=1,atoms%nat
            rat_int(1,iat)=modulo(modulo(rat_int(1,iat),1.d0),1.d0)
            rat_int(2,iat)=modulo(modulo(rat_int(2,iat),1.d0),1.d0)
        enddo
    elseif(trim(atoms%boundcond)=='wire') then
        do iat=1,atoms%nat
            rat_int(1,iat)=modulo(modulo(rat_int(1,iat),1.d0),1.d0)
        enddo
    elseif(trim(atoms%boundcond)=='free') then
    else
        write(*,'(2a)') 'ERROR: unknown BC in linkedlists_init ',trim(atoms%boundcond)
        stop
    endif
    jat=1
    do iz=1,mz+mlimnb3
        ishiftz=(isign(1,iz-1)-isign(1,mz-iz))/2
        do iy=1-mlimnb2,my+mlimnb2
            ishifty=(isign(1,iy-1)-isign(1,my-iy))/2
            ishift_l=0
            ishift_r=0
            if(ishifty<0 .and. iz<=mz ) ishift_l=1
            if(iz>mz) ishift_r=-mlimnb1
            do ix=1+ishift_l*mx+ishift_r,mx+mlimnb1
                ishiftx=(isign(1,ix-1)-isign(1,mx-ix))/2
                iat=linked_lists%head(modulo(ix-ishiftx*mx-1,mx)+1,modulo(iy-ishifty*my-1,my)+1,modulo(iz-ishiftz*mz-1,mz)+1)
                if(iat==0) then
                    linked_lists%last(ix,iy,iz)=linked_lists%last(ix-1,iy,iz)
                    cycle
                endif
                linked_lists%prime(ix,iy,iz)=jat
                do
                    if(iat==0) exit
                    nx=floor((ix-1)/real(mx))
                    ny=floor((iy-1)/real(my))
                    nz=floor((iz-1)/real(mz))
                    rfrac(1)=rat_int(1,iat)+real(nx,kind=8)
                    rfrac(2)=rat_int(2,iat)+real(ny,kind=8)
                    rfrac(3)=rat_int(3,iat)+real(nz,kind=8)
                    rcart=matmul(atoms%cellvec,rfrac)
                    linked_lists%ratp(1,jat)=rcart(1)
                    linked_lists%ratp(2,jat)=rcart(2)
                    linked_lists%ratp(3,jat)=rcart(3)
                    linked_lists%itypat(jat)=atoms%itypat(iat)
                    linked_lists%perm(jat)=iat
                    linked_lists%qat(jat)=atoms%qat(iat)
                    iat=linked_lists%list(iat)
                    if(ishiftx==0 .and. ishifty==0 .and. ishiftz==0) linked_lists%maincell(jat)=1
                    jat=jat+1
                enddo
                linked_lists%last(ix,iy,iz)=jat-1
            enddo
            do ix=mx+mlimnb1+1,mx+mlimnb1+mlimnb
                linked_lists%last(ix,iy,iz)=linked_lists%last(ix-1,iy,iz)
            enddo
            do ix=mx+mlimnb1+mlimnb-1,1-mlimnb-bulkbc*mlimnb1,-1   
                if(linked_lists%prime(ix,iy,iz)==natpopnatim) then
                    linked_lists%prime(ix,iy,iz)=linked_lists%prime(ix+1,iy,iz)
                endif
            enddo
        enddo !end of loop over iy
    enddo !end of loop over iz
    call update_rat(linked_lists%typ_atoms,upall=.true.)
    !write(*,*)jat,"main"
    !When cell is small such that the number of cells in a direction is
    !small and natim's (number of atoms in segments explained in top of the file)
    !are incorrect the linked list in the current version does not work and
    !the following error message will be shown.
    if(jat/=natpopnatim) then
        !The following is commented since it is not in fact an error.
        !However, we must find out if we can properly count natpopnatim so
        !the message shows up only there really is a problem.
        !write(*,'(a,2i8)') 'ERROR: in number of image atoms: ',jat,natpopnatim
        linked_lists%natim=jat-1
        !stop
    endif
    deallocate(rat_int)
    deallocate(rat)
    !write(*,*) 'jat-1',jat-1
end subroutine prepprimelast
!*****************************************************************************************
!head: heads of cells
!list: lists of particles in cells
subroutine make_list_new(parini,atoms,linked_lists,cell)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, get_rat
    use mod_electrostatics, only: typ_linked_lists
    implicit none
    type(typ_parini), intent(in):: parini 
    type(typ_atoms), intent(in):: atoms
    type(typ_linked_lists), intent(inout):: linked_lists
    real(8), intent(in):: cell(3)
    !local variables
    real(8):: hxinv, hyinv, hzinv, nat7
    integer:: natim, natimarr(9),tt
    integer:: natim1, natim2
    integer:: iat, ix, iy, iz, ind, ishift1, ishift2, ishift3
    integer:: mx, my, mz, mlimnb1, mlimnb2, mlimnb3
    real(8), allocatable:: rat_int(:,:)
    real(8), allocatable:: rat(:,:)
    allocate(rat_int(1:3,1:atoms%nat))
    allocate(rat(1:3,1:atoms%nat))
    mx=linked_lists%mx
    my=linked_lists%my
    mz=linked_lists%mz
    mlimnb1=linked_lists%mlimnb1
    mlimnb2=linked_lists%mlimnb2
    mlimnb3=linked_lists%mlimnb3
    !compute the subcell dimensions
    hxinv=real(mx,8)/cell(1)
    hyinv=real(my,8)/cell(2)
    hzinv=real(mz,8)/cell(3)
    !initialize the heads
    linked_lists%head=0
    !form the cell list
    call get_rat(atoms,rat)
    call rxyz_cart2int_alborz(atoms%nat,atoms%cellvec,rat,rat_int)
    !we are doing the following since in mutiple application of h and h^0-1, i.e.,
    !in going from reduced to Cartesian and vice versa, there can be rounding
    !which can put atoms outside cell, e.g. 0 can be -2.E-17
    if(trim(atoms%boundcond)=='bulk') then
        do iat=1,atoms%nat
            rat_int(1,iat)=modulo(modulo(rat_int(1,iat),1.d0),1.d0)
            rat_int(2,iat)=modulo(modulo(rat_int(2,iat),1.d0),1.d0)
            rat_int(3,iat)=modulo(modulo(rat_int(3,iat),1.d0),1.d0)
        enddo
    elseif(trim(atoms%boundcond)=='slab') then
        do iat=1,atoms%nat
            rat_int(1,iat)=modulo(modulo(rat_int(1,iat),1.d0),1.d0)
            rat_int(2,iat)=modulo(modulo(rat_int(2,iat),1.d0),1.d0)
        enddo
    elseif(trim(atoms%boundcond)=='wire') then
        do iat=1,atoms%nat
            rat_int(1,iat)=modulo(modulo(rat_int(1,iat),1.d0),1.d0)
        enddo
    elseif(trim(atoms%boundcond)=='free') then
    else
        write(*,'(2a)') 'ERROR: unknown BC in linkedlists_init ',trim(atoms%boundcond)
        stop
    endif
    do iat=1,atoms%nat
        ix=floor(rat_int(1,iat)*mx)+1
        iy=floor(rat_int(2,iat)*my)+1
        iz=floor(rat_int(3,iat)*mz)+1
        linked_lists%list(iat)=linked_lists%head(ix,iy,iz)
        linked_lists%head(ix,iy,iz)=iat
    enddo
    linked_lists%natim=atoms%nat
    if(trim(atoms%boundcond)=='free') then
        natim=0
        deallocate(rat_int)
        return
    endif
    natimarr(1)=0 ; natimarr(2)=0 ; natimarr(3)=0 ; natimarr(4)=0 ; natimarr(5)=0
    natimarr(6)=0 ; natimarr(7)=0 ; natimarr(8)=0 ; natimarr(9)=0 
    natim1=0
    do iz=1,min(mlimnb3,mz)
        do iy=1,my
            ishift1=(isign(1,iy)-isign(1,iy-mlimnb2-1))/2
            ishift2=(isign(1,iy)-isign(1,iy-(my-mlimnb2+1)))/2
            ishift3=(isign(1,iy)-isign(1,my-mlimnb2-iy))/2
            do ix=1,mx
                iat=linked_lists%head(ix,iy,iz)
                if(iat==0) cycle
                if(ix <= mx-mlimnb1) then
                    ind=7+2*ishift3
                else
                    ind=6+2*ishift3
                endif
                tt=0
                do 
                    if(iat==0) exit
                    tt=tt+1
                    iat=linked_lists%list(iat)
                enddo
                natimarr(ind)=natimarr(ind)+tt
                if (ishift1==1 .and. ix >mx-mlimnb1) natimarr(2)=natimarr(2)+tt
            enddo
        enddo
    enddo
    natim2=4*natimarr(8)+2*natimarr(6)+natimarr(7)+2*natimarr(9)+natimarr(2)
    natim1=0
    natimarr(1)=0
    natimarr(3:5)=0
    nat7=natimarr(7)
    natimarr=0
    do iz=1,mz
        do iy=1,my
            ishift1=(isign(1,iy)-isign(1,iy-mlimnb2-1))/2
            ishift2=(isign(1,iy)-isign(1,mlimnb2-iy))/2
            ishift3=(isign(1,iy)-isign(1,my-mlimnb2-iy))/2
            do ix=1,min(mx*ishift1+mlimnb1*ishift2,mx)
                iat=linked_lists%head(ix,iy,iz)
                if(iat==0) cycle
                if(ix<mlimnb1+1) then
                    ind=1+ishift2*3
                else
                    ind=5
                endif
                tt=0
                do 
                    if(iat==0) exit
                    tt=tt+1
                    iat=linked_lists%list(iat)
                enddo
                natimarr(ind)=natimarr(ind)+tt
                if (ishift3==1 .and. ix<mlimnb1+1) natimarr(3)=natimarr(3)+tt
            enddo
        enddo
        if (iz<mlimnb3+1) then
            natim1=natim1+3*natimarr(1)+natimarr(3)+natimarr(4)+natimarr(5)
            natimarr=0
        endif
    enddo
    if(parini%iverbose>2) then
        write(*,'(a,5i6)') 'natim1,natim2,natim3,natim4,natim5 ', &
            natimarr(1),natimarr(2),natimarr(3),natimarr(4),natimarr(5)
    endif
    natim=3*natimarr(1)+natimarr(3)+natimarr(4)+natimarr(5)+2*natim1+natim2
    if (mx< mlimnb1) natim=(mlimnb1/mx+1)*(natim) !-nat7) to be check later
    if (my< mlimnb2) natim=(mlimnb2/my+1)*(natim) !-nat7) to be check later
    !if (mz< mlimnb3) natim2=(mlimnb3/mz+1)*(natim2+natim1)
    if (mz< mlimnb3) natim=(mlimnb3/mz+1)*natim
   ! natim=natim+natim2
    if(parini%iverbose>2) then
        write(*,'(a,i6)') 'natim= ',natim
    endif
    linked_lists%natim=linked_lists%natim+natim
    deallocate(rat_int)
    deallocate(rat)
end subroutine make_list_new
!*****************************************************************************************
!This subroutine determines the bounds of ix and iy to count
!all cells in the sphere given by the cutoff radius.
!limnbx(1,:,:) is the lower bound for the cell neighbors in x directions.
!limnbx(2,:,:) is the upper bound for the cell neighbors in x directions.
!iz=0 is the plane which intersects the central cell and determining
!the bounds of iz=0 requires special care.
!iz=0,mlimnb because we count only cell which are after the central cell
!and all cells below the plane iz=0 are considered before the central cell.
subroutine determine_sclimitsphere(linked_lists)
    use mod_electrostatics, only: typ_linked_lists
    implicit none
    type(typ_linked_lists), intent(inout):: linked_lists
    !local variables
    integer:: iy, ix, iz, iyt, izt, mlimnb
    mlimnb=linked_lists%mlimnb
    !First setting the lower bound a larger value than the upper bound.
    do iz=0,mlimnb
        do iy=-mlimnb,mlimnb
            linked_lists%limnbx(1,iy,iz)=1
            linked_lists%limnbx(2,iy,iz)=0
        enddo
    enddo
    !The following three nested loops goes through a quarter of all cells in
    !the cube to find out the cells in the sphere.
    do iz=0,mlimnb
        izt=iz
        if(iz==0) izt=iz+1
        izt=abs(izt)-1
        do iy=-mlimnb,mlimnb
            iyt=iy
            if(iy==0) iyt=iy+1
            iyt=abs(iyt)-1
            do ix=mlimnb,0,-1
                if((abs(ix)-1)**2+iyt**2+izt**2<mlimnb**2) then
                    linked_lists%limnbx(1,iy,iz)=-ix
                    linked_lists%limnbx(2,iy,iz)=ix
                    exit
                endif
            enddo
        enddo
    enddo
    !In the row in plane iz=0 which goes through the central cell, i.e. iy=0,
    !cells after central cell are from 0 to mlimnb. This overwrites the
    !bounds determined in the three nested loops above.
    linked_lists%limnbx(1,0,0)=0
    linked_lists%limnbx(2,0,0)=mlimnb
    !In the plane iz=0, rows with iy<0 are before the central cell and
    !the bounds are set such that to exclude them. This assignment is in fact
    !done in the beginning loops but since the bounds are overwritten in the
    !three nested loops, the assignment must be repeated again.
    do iy=-mlimnb,-1
        linked_lists%limnbx(1,iy,0)=1
        linked_lists%limnbx(2,iy,0)=0
    enddo
    !Now all bounds of ix are determined and the bounds of iy
    !must be determined.
    !In the plane iz=0, iy must be between 0 and mlimnb.
    linked_lists%limnby(1,0)=0
    linked_lists%limnby(2,0)=mlimnb
    !The following loops is conceptually heavily modified and
    !requires to be checked by Samare.
    do iz=1,mlimnb
        linked_lists%limnby(1,iz)=1 !-mlimnb
        linked_lists%limnby(2,iz)=0 !mlimnb
    enddo
    !Determining the bounds of iy for planes above plane iz=0,
    !if the limnbx lower bound is smaller than or equal to the
    !upper bound the bounds of iy must be updated.
    do iz=1,mlimnb
        do iy=1,mlimnb
            if(linked_lists%limnbx(1,iy,iz)<=linked_lists%limnbx(2,iy,iz)) then
                linked_lists%limnby(1,iz)=-iy
                linked_lists%limnby(2,iz)=iy
            endif
        enddo
    enddo
    !do iz=0,mlimnb
    !write(*,'(3i5)') iz,linked_lists%limnby(1,iz),linked_lists%limnby(2,iz)
    !enddo
    !do iz=0,mlimnb
    !do iy=-mlimnb,mlimnb
    !write(*,'(4i5)') iy,iz,linked_lists%limnbx(1,iy,iz),linked_lists%limnbx(2,iy,iz)
    !enddo
    !enddo
end subroutine determine_sclimitsphere
!*****************************************************************************************
!dbl_count if .true., bonds are double counted.
subroutine call_linkedlist(parini,atoms,dbl_count,linked_lists,pia_arr)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, type_pairs, update_ratp
    use mod_const, only: bohr2ang
    use mod_linked_lists, only: typ_linked_lists, typ_pia_arr
    implicit none
    type(typ_parini), intent(in):: parini 
    type(typ_atoms), intent(in):: atoms 
    logical, intent(in):: dbl_count
    type(typ_linked_lists), intent(inout):: linked_lists
    type(typ_pia_arr), intent(inout):: pia_arr
    !local variables
    integer::nimat ,iat_maincell, jat_maincell
    integer:: istat
    integer:: iat, jat, niat,kat,kkz,conf
    integer:: ix,iy,iz,jy,jz,kx,ky,kz,dkjz,dkiy 
    integer:: llx,mmx,lly,mmy ,jp,jl,kp,kl,ip,il ,kpt,jpt
    integer:: nmax, maxnbr, maxnba, maincell_iat
    integer:: ibr, njat, inbr, jnbr, n, ntot
    real(8):: sclinv,cell(3) ,rcutsq
    real(8):: xiat, yiat, ziat
    real(8):: rij, rijsq, drij(3), rik, riksq, drik(3) , rjk, rjksq, drjk(3)
    real(8):: r3, dx3, dy3, dz3 ,hxinv
    integer, allocatable:: bound_rad(:,:)
    real(8), allocatable:: bound_dist(:,:,:)
    integer, allocatable:: neighbor(:)
    logical :: yes
    call linkedlists_init(parini,atoms,cell,linked_lists)
    nmax=1000
    !if (.not. linked_lists%triplex) then
        !allocate(bound_rad(2,min(linked_lists%nat*namx,linked_lists%nat**2)))
        !allocate(bound_dist(4,min(linked_lists%nat*nmax,linked_lists%nat**2)),1)
    !else
        allocate(bound_rad(1:nmax,1:atoms%nat))
        allocate(bound_dist(1:4,1:nmax,1:atoms%nat))
        allocate(linked_lists%prime_bound(1:atoms%nat+1))
        allocate(neighbor(1:atoms%nat))
    !endif

    rcutsq=linked_lists%rcut**2
    maxnbr=0
    !-------------------------------------------------------
    associate(mx=>linked_lists%mx)
    associate(my=>linked_lists%my)
    neighbor(:)=0
    bound_rad=0
    bound_dist=0
    call update_ratp(linked_lists%typ_atoms)
    do iz=1,linked_lists%mz+linked_lists%mlimnb3
    do iy=1-linked_lists%mlimnb2,linked_lists%my+linked_lists%mlimnb2
    do ix=1-((isign(1,iy-1)-1)/2)*linked_lists%mx,linked_lists%mx+linked_lists%mlimnb1
    ip=linked_lists%prime(ix,iy,iz)
    il=linked_lists%last(ix,iy,iz)
    do  iat=ip,il
        xiat=linked_lists%ratp(1,iat)
        yiat=linked_lists%ratp(2,iat)
        ziat=linked_lists%ratp(3,iat)
        maincell_iat=linked_lists%maincell(iat)
        do jz=iz,min(linked_lists%mz+linked_lists%mlimnb3,iz+linked_lists%mlimnb)
        do jy=iy+linked_lists%limnby(1,jz-iz),iy+linked_lists%limnby(2,jz-iz)
        jpt=linked_lists%prime(ix+linked_lists%limnbx(1,jy-iy,jz-iz),jy,jz)
        jp=(iat-ip+1)*((isign(1,ip-jpt)+1)/2)+jpt
        jl=linked_lists%last(ix+linked_lists%limnbx(2,jy-iy,jz-iz),jy,jz)
        do  jat=jp,jl
            drij(1)=linked_lists%ratp(1,jat)-xiat
            drij(2)=linked_lists%ratp(2,jat)-yiat
            drij(3)=linked_lists%ratp(3,jat)-ziat
            rijsq=drij(1)**2+drij(2)**2+drij(3)**2
            if (rijsq< rcutsq .and. (maincell_iat+linked_lists%maincell(jat)) >-1 ) then
                    rij=sqrt(rijsq)
                    maxnbr=maxnbr+1
                    iat_maincell=linked_lists%perm(iat)
                    jat_maincell=linked_lists%perm(jat)

                    neighbor(iat_maincell)=neighbor(iat_maincell)+1
                    bound_rad(neighbor(iat_maincell),iat_maincell)=jat_maincell
                    if (neighbor(iat_maincell) > nmax  ) then
                        write(*,*) " neighbours are more that expected  " 
                        stop
                    endif
                    bound_dist(1,neighbor(iat_maincell),iat_maincell)=rij
                    bound_dist(2,neighbor(iat_maincell),iat_maincell)=drij(1)
                    bound_dist(3,neighbor(iat_maincell),iat_maincell)=drij(2)
                    bound_dist(4,neighbor(iat_maincell),iat_maincell)=drij(3)
                    ! if (iat_maincell==jat_maincell) cycle
                    if(dbl_count) then
                        neighbor(jat_maincell)=neighbor(jat_maincell)+1
                        if (neighbor(jat_maincell) > nmax ) then
                            write(*,*)" neighbours are more that expected  "
                            stop
                        endif
                        bound_rad(neighbor(jat_maincell),jat_maincell)=iat_maincell
                        bound_dist(1,neighbor(jat_maincell),jat_maincell)=rij
                        bound_dist(2,neighbor(jat_maincell),jat_maincell)=-drij(1)
                        bound_dist(3,neighbor(jat_maincell),jat_maincell)=-drij(2)
                        bound_dist(4,neighbor(jat_maincell),jat_maincell)=-drij(3)
                    endif
            endif
        enddo !end of loop over jat
        enddo !end of loop over jy
        enddo !end of loop over jz
        !---------------------------------------------------
    enddo !end of loop over iat
    enddo !end of loop over ix
    enddo !end of loop over iy
    enddo !end of loop over iz
    if(dbl_count) then
        linked_lists%maxbound_rad=maxnbr*2
    else
        linked_lists%maxbound_rad=maxnbr
    endif
    allocate(linked_lists%bound_rad(1:2,1:linked_lists%maxbound_rad))
    allocate(pia_arr%pia(linked_lists%maxbound_rad))
    njat=0
    ibr=0
    iat=0
    do iat=1,atoms%nat
        linked_lists%prime_bound(iat)=ibr+1
        do njat=1,neighbor(iat)
            !The following if added to avoid double counting for bond based linked list
            if(parini%bondbased_ann .and. iat>bound_rad(njat,iat)) then
                cycle
            endif
        !do 
        !    njat=njat+1
        !    if (bound_rad(njat,iat)<1 .or. njat>nmax) then
        !        njat=0
        !        exit
        !    endif
            ibr=ibr+1
            linked_lists%bound_rad(1,ibr)=iat
            linked_lists%bound_rad(2,ibr)=bound_rad(njat,iat)
            pia_arr%pia(ibr)%r=bound_dist(1,njat,iat)
            pia_arr%pia(ibr)%dr(1)=bound_dist(2,njat,iat)
            pia_arr%pia(ibr)%dr(2)=bound_dist(3,njat,iat)
            pia_arr%pia(ibr)%dr(3)=bound_dist(4,njat,iat)
        enddo
    enddo
    !The following if added to avoid double counting for bond based linked list
    !and it must be corrected since it will not work if repulsive atom based
    !linked list is needed. The problem is due to parini%bondbased_ann
    if(parini%bondbased_ann .and. iat>njat) then
        if(mod(linked_lists%maxbound_rad,2)/=0) then
            stop 'ERROR: linked_lists%maxbound_rad must be even.'
        endif
        linked_lists%maxbound_rad=linked_lists%maxbound_rad/2
    endif
    if (ibr/=linked_lists%maxbound_rad) then
        write(*,'(a,2i8)') 'ERROR: in number of bonds ',ibr,linked_lists%maxbound_rad
        stop
    endif
    linked_lists%prime_bound(atoms%nat+1)=linked_lists%maxbound_rad+1
    deallocate(bound_rad)
    deallocate(bound_dist)
    ntot=0
    do iat=1,atoms%nat
        n=linked_lists%prime_bound(iat+1)-linked_lists%prime_bound(iat)
        ntot=ntot+(n*(n-1))/2.d0
    enddo
    linked_lists%maxbound_ang=ntot
    allocate(linked_lists%bound_ang(1:2,1:linked_lists%maxbound_ang))
    maxnba=0
    do iat=1,atoms%nat
        do inbr=linked_lists%prime_bound(iat),linked_lists%prime_bound(iat+1)-1
        do jnbr=inbr+1,linked_lists%prime_bound(iat+1)-1
            maxnba=maxnba+1
            linked_lists%bound_ang(1,maxnba)=inbr
            linked_lists%bound_ang(2,maxnba)=jnbr
        enddo
        enddo
    enddo
    if (maxnba/=linked_lists%maxbound_ang) then
        write(*,'(a,2i8)') 'ERROR: in number of angular bonds ',maxnba,linked_lists%maxbound_ang
        stop
    endif
    end associate
    end associate
    call linkedlists_final(linked_lists)
    deallocate(neighbor)
end subroutine call_linkedlist 
!***************************************************************************************************
