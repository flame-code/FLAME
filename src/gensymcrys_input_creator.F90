!*****************************************************************************************
subroutine gensymcrys_many(parini)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, typ_atoms_arr, atom_allocate_old, atom_deallocate_old
    use mod_atoms, only: atom_copy_old, set_rat, typ_file_info
    use mod_yaml_conf, only: write_yaml_conf
    use mod_const, only: bohr2ang
    implicit none
    type(typ_parini), intent(in):: parini
    !local variables
    integer, allocatable:: cryssys_all(:), brav_all(:), nsymp_all(:), nsym_all(:), ind_rsym_all(:)
    real(8), allocatable:: rsym_all(:,:,:)
    integer:: nsym_tot, ntry, i_spacegroup, NCELLS, NAT_CELL, NKINDS, i, j, ios, nfu
    integer:: iat, iconf, nconf, nconf_allproc, ierr
    integer:: ifu, NAT_CELL_MAX
    real(8):: target_vol_per_atom, ttrand, tt1, tt2
    character(5), allocatable:: stypat(:)
    integer, allocatable:: NAT_KINDS(:)
    real(8), allocatable:: KINDSDIST_MIN(:,:)
    real(8), allocatable:: rxyz_all(:,:,:)
    real(8), allocatable:: cv_all(:,:,:)
    character(5), allocatable:: sat_all(:,:)
    integer, allocatable:: NAT_CELL_ALL(:)
    real(8), allocatable:: rxyz_allproc(:,:,:)
    real(8), allocatable:: cv_allproc(:,:,:)
    character(5), allocatable:: sat_allproc(:,:)
    integer, allocatable:: NAT_CELL_ALLPROC(:)
    logical:: succeeded
    type(typ_atoms):: atoms
    type(typ_file_info):: file_info
#if defined(MPI)
    include 'mpif.h'
#endif
    NKINDS=parini%ntypat
    allocate(stypat(NKINDS))
    allocate(NAT_KINDS(NKINDS))
    allocate(KINDSDIST_MIN(NKINDS,NKINDS))
    write(*,*) "Generating primitive cell"
    KINDSDIST_MIN(1:NKINDS,1:NKINDS)=parini%rmin_pairs(1:NKINDS,1:NKINDS)
    do i=1,NKINDS
        stypat(i)=parini%stypat(i)
    enddo
    nfu=size(parini%list_fu)
    NAT_CELL_MAX=0
    do ifu=1,nfu
    do i=1,NKINDS
        NAT_KINDS(i)=parini%nat_types_fu(i)*parini%list_fu(ifu)
    enddo
    NAT_CELL=sum(NAT_KINDS(1:NKINDS))
    NAT_CELL_MAX=max(NAT_CELL,NAT_CELL_MAX)
    enddo
    ntry=parini%ntry
    ncells=1 !parini%ncells !probably values other than 1 will not work
    !---------------------------------------------------------------------
    nsym_tot=4425
    allocate(cryssys_all(230),brav_all(230),nsymp_all(230),nsym_all(230),ind_rsym_all(230))
    allocate(rsym_all(4,4,nsym_tot))
    if(parini%mpi_env%iproc==0) then
    call read_coeffs_gensymcrys(parini,nsym_tot,rsym_all,cryssys_all,brav_all,nsymp_all,nsym_all,ind_rsym_all)
    endif
    call MPI_BARRIER(parini%mpi_env%mpi_comm,ierr)
    call MPI_BCAST(rsym_all,4*4*nsym_tot,MPI_DOUBLE_PRECISION,0,parini%mpi_env%mpi_comm,ierr)
    call MPI_BCAST(cryssys_all,230,MPI_INTEGER,0,parini%mpi_env%mpi_comm,ierr)
    call MPI_BCAST(brav_all,230,MPI_INTEGER,0,parini%mpi_env%mpi_comm,ierr)
    call MPI_BCAST(nsymp_all,230,MPI_INTEGER,0,parini%mpi_env%mpi_comm,ierr)
    call MPI_BCAST(nsym_all,230,MPI_INTEGER,0,parini%mpi_env%mpi_comm,ierr)
    call MPI_BCAST(ind_rsym_all,230,MPI_INTEGER,0,parini%mpi_env%mpi_comm,ierr)

    !---------------------------------------------------------------------
    nconf=parini%nconf_genconf/parini%mpi_env%nproc
    if(nconf*parini%mpi_env%nproc/=parini%nconf_genconf) nconf=nconf+1

    allocate(rxyz_all(3,NAT_CELL_MAX,nconf),sat_all(NAT_CELL_MAX,nconf))
    allocate(NAT_CELL_ALL(nconf),cv_all(3,3,nconf))

    succeeded=.true.
    iconf=0
    do
    if(succeeded) iconf=iconf+1
    if(parini%ispg<1) then
        call random_number(ttrand)
        i_spacegroup=min(int(ttrand*real(230,kind=8))+1,230)
    else
        i_spacegroup=parini%ispg
    endif
    if(.not. succeeded) i_spacegroup=1
    call random_number(ttrand)
    tt1=parini%volperatom_bounds(1)
    tt2=parini%volperatom_bounds(2)
    target_vol_per_atom=(tt2-tt1)*ttrand+tt1
    call random_number(ttrand)
    ifu=min(int(ttrand*real(nfu,kind=8))+1,nfu)
    do i=1,NKINDS
        NAT_KINDS(i)=parini%nat_types_fu(i)*parini%list_fu(ifu)
    enddo
    NAT_CELL_ALL(iconf)=sum(NAT_KINDS(1:NKINDS))
    call gensymcrys_single(parini,nsym_tot,rsym_all,cryssys_all,brav_all,nsymp_all, &
        nsym_all,ind_rsym_all,target_vol_per_atom,stypat,NCELLS,NAT_CELL_ALL(iconf),ntry, &
        i_spacegroup,NKINDS,NAT_KINDS,KINDSDIST_MIN,cv_all(1,1,iconf),rxyz_all(1,1,iconf),sat_all(1,iconf),succeeded)
    if(succeeded .and. iconf==nconf) exit
    enddo
    call MPI_BARRIER(parini%mpi_env%mpi_comm,ierr)

    nconf_allproc=nconf*parini%mpi_env%nproc
    allocate(rxyz_allproc(3,NAT_CELL_MAX,nconf_allproc),sat_allproc(NAT_CELL_MAX,nconf_allproc))
    allocate(NAT_CELL_ALLPROC(nconf_allproc),cv_allproc(3,3,nconf_allproc))

    call MPI_GATHER(NAT_CELL_ALL,nconf,MPI_INTEGER,NAT_CELL_ALLPROC,nconf,MPI_INTEGER,0,parini%mpi_env%mpi_comm,ierr)
    call MPI_GATHER(rxyz_all,3*NAT_CELL_MAX*nconf,MPI_DOUBLE_PRECISION,rxyz_allproc,3*NAT_CELL_MAX*nconf,MPI_DOUBLE_PRECISION,0,parini%mpi_env%mpi_comm,ierr)
    call MPI_GATHER(cv_all,3*3*nconf,MPI_DOUBLE_PRECISION,cv_allproc,3*3*nconf,MPI_DOUBLE_PRECISION,0,parini%mpi_env%mpi_comm,ierr)
    call MPI_GATHER(sat_all,5*NAT_CELL_MAX*nconf,MPI_CHARACTER,sat_allproc,5*NAT_CELL_MAX*nconf,MPI_CHARACTER,0,parini%mpi_env%mpi_comm,ierr)

    if(parini%mpi_env%iproc==0) then
    file_info%filename_positions='posout.yaml'
    do iconf=1,nconf_allproc
        call atom_allocate_old(atoms,NAT_CELL_ALLPROC(iconf),0,0)
        rxyz_allproc(1:3,1:atoms%nat,iconf)=rxyz_allproc(1:3,1:atoms%nat,iconf)/bohr2ang
        atoms%cellvec(1:3,1:3)=cv_allproc(1:3,1:3,iconf)/bohr2ang
        call set_rat(atoms,rxyz_allproc(1,1,iconf),.true.)
        atoms%boundcond='bulk'
        atoms%units_length_io='angstrom'
        do iat=1,atoms%nat
            atoms%sat(iat)=sat_allproc(iat,iconf)
        enddo
        if(iconf==1) then
            file_info%file_position='new'
        else
            file_info%file_position='append'
        endif
        call write_yaml_conf(file_info,atoms=atoms,strkey='posout')
        call atom_deallocate_old(atoms)
    enddo
    endif

    call MPI_BARRIER(parini%mpi_env%mpi_comm,ierr)
    deallocate(rxyz_all,sat_all,cv_all,NAT_CELL_ALL)
    deallocate(rxyz_allproc,sat_allproc,cv_allproc,NAT_CELL_ALLPROC)
    !---------------------------------------------------------------------
    deallocate(cryssys_all,brav_all,nsymp_all,nsym_all,ind_rsym_all)
    deallocate(rsym_all)
    !---------------------------------------------------------------------
    deallocate(NAT_KINDS)
    deallocate(KINDSDIST_MIN)
end subroutine gensymcrys_many
!*****************************************************************************************
subroutine gensymcrys_single(parini,nsym_tot,rsym_all,cryssys_all,brav_all,nsymp_all, &
        nsym_all,ind_rsym_all,target_vol_per_atom,stypat,NCELLS,NAT_CELL,ntry, &
        i_spacegroup,NKINDS,NAT_KINDS,KINDSDIST_MIN,cv,rxyz,sat,succeeded)
    use mod_parini, only: typ_parini
    implicit none
    type(typ_parini), intent(in):: parini
    integer, intent(in):: nsym_tot
    real(8), intent(in):: rsym_all(4,4,nsym_tot)
    integer, intent(in):: cryssys_all(230), brav_all(230)
    integer, intent(in):: nsymp_all(230), nsym_all(230), ind_rsym_all(230)
    real(8), intent(in):: target_vol_per_atom
    character(5), intent(in):: stypat(3)
    integer, intent(in):: NCELLS, NAT_CELL, ntry, i_spacegroup
    integer, intent(in):: NKINDS, NAT_KINDS(NKINDS)
    real(8), intent(in):: KINDSDIST_MIN(NKINDS,NKINDS)
    real(8), intent(out):: cv(3,3), rxyz(3,NAT_CELL)
    character(5), intent(out):: sat(NAT_CELL)
    logical, intent(out):: succeeded
    !local variables
integer:: LATSGP,NSYM,NSYMP,I,CRYSSYS,BRAV,NGUESS,NSPEC,J,K,L,M
integer, parameter::NSYMMAX=192
real(8):: RSYM(4,4,NSYMMAX),trans(3,3),latvec(3,3),latvec_prim(3,3),target_vol,v(3,3),vol,rotmat(3,3),dproj(6),rand,tmp(3),dout(3)
real(8),allocatable::RED_POS(:,:)
integer:: cells(3),t_cells(3)
integer,allocatable:: KINDS(:)
real(8),allocatable:: ratred(:,:)
real(8),allocatable:: RED_POS_PRIM(:,:)
integer::NPOS_IRRED
integer,allocatable:: NPOS_IRRED_ARR(:)
integer::  isize,idate(8),order(1), iat, itry
real(8):: vol_prim
succeeded=.false.
NGUESS=100
allocate(RED_POS(3,NGUESS))
allocate(NPOS_IRRED_ARR(NSYMMAX))

if(i_spacegroup<0 .or. i_spacegroup>230) then
    stop 'ERROR: invalid value for i_spacegroup'
elseif(i_spacegroup==0) then
    call random_number(rand)
    LATSGP=int(rand*230)+1
else
    LATSGP=i_spacegroup
endif
write(*,*) "Spacegroup", LATSGP

call optcell(ncells,cells)

allocate(KINDS(NAT_CELL))
!allocate(RXYZ(3,NAT_CELL))
allocate(ratred(3,NAT_CELL))
allocate(RED_POS_PRIM(3,NAT_CELL))
!NAT_KINDS(1)=NAT_CELL
!if(NKINDS.gt.1) then
!do i=1,NKINDS
   !write(*,'(a24,i2)') "Number of atoms of kind ",i
!   read(11,*) NAT_KINDS(1:NKINDS)
!enddo
!endif
target_vol=target_vol_per_atom*NAT_CELL
!---------------------------------------------------------------------

itry=0
1050 continue
itry=itry+1
if(itry>ntry) then
    !write(*,'(a,i5,a)') 'ERROR: failed to find a configuration after ',ntry,' tries.'
    return !stop
endif

call sg_ops(nsym_tot,rsym_all,cryssys_all,brav_all,nsymp_all,nsym_all,ind_rsym_all, &
    LATSGP,CRYSSYS,BRAV,NSYMP,NSYM,RSYM,NSYMMAX)
write(*,*) "LATSGP",LATSGP
write(*,*) "NSYMP",NSYMP
write(*,*) "NSYM",NSYM
WRITE(*,*)"BRAV",BRAV
write(*,*) "NGUESS",NGUESS
call random_incell_main(LATSGP,CRYSSYS,NGUESS,RED_POS)

!Generating the random primitive unit cell
call random_lattice(LATSGP,CRYSSYS,BRAV,LATVEC,LATVEC_PRIM,TARGET_VOL)


call find_ir_multiplicity(LATSGP,RSYM,NSYMP,NSYM,BRAV,NSYMMAX,NGUESS,NPOS_IRRED,NPOS_IRRED_ARR)
write(*,*) "Irreducible Multiplicity"
do i=1,NPOS_IRRED
   write(*,*)NPOS_IRRED_ARR(i)
enddo
do j=1,NKINDS
   if(NAT_KINDS(j).lt.minval(NPOS_IRRED_ARR(1:NPOS_IRRED))) then
      write(*,'(a4,i3,a10,i3,a41,i3)') "NAT ",NAT_KINDS(j) ," of KINDS ",j," is smaller than the irred. Multiplicity ",minval(NPOS_IRRED_ARR(1:NPOS_IRRED))
   endif
enddo
do j=1,NKINDS
   do i=1,NPOS_IRRED
         if(NAT_KINDS(j).ge.NPOS_IRRED_ARR(i)) then
           if(modulo(NAT_KINDS(j),NPOS_IRRED_ARR(i))==0) then
             write(*,'(a10,i3,a)') "Atom type ",j," looks fine" 
             exit
           endif
         endif
      write(*,'(a10,i3,a)') "Atom type ",j," seems not to be solvable with this space group" 
     write(*,*) "Going back to start" 
     goto 1050
   enddo
enddo



call random_atom(LATSGP,CRYSSYS,BRAV,NSYMP,NSYM,RSYM,NSYMMAX,NAT_CELL,NAT_KINDS,NKINDS,KINDS,NGUESS,KINDSDIST_MIN,RXYZ,RED_POS_PRIM,LATVEC_PRIM)
if(kinds(1)==0) then
write(*,*) "random atom unsuccessful"
goto 1050
endif

!Multiplicate the primitive cell if necesarry
call dist_dim(latvec_prim,dout)
!Sorting the "cells" array to inverse to the sorting of dout
write(*,*) "dout",dout
write(*,*) "cells",cells
do i=1,3
order=minloc(dout(:))
t_cells(order(1))=cells(i)
dout(order(1))=10000000
enddo
latvec(:,1)=latvec_prim(:,1)*t_cells(1)
latvec(:,2)=latvec_prim(:,2)*t_cells(2)
latvec(:,3)=latvec_prim(:,3)*t_cells(3)
write(*,*) "t_cells",t_cells

succeeded=.true.
write(*,'(a,i5)') 'Succeeded: itry= ',itry

call rxyz_cart2int_gensymcrys(latvec_prim,ratred,RXYZ,NAT_CELL)
do iat=1,NAT_CELL
    ratred(1,iat)=modulo(modulo(ratred(1,iat),1.d0),1.d0)
    ratred(2,iat)=modulo(modulo(ratred(2,iat),1.d0),1.d0)
    ratred(3,iat)=modulo(modulo(ratred(3,iat),1.d0),1.d0)
enddo
call rxyz_int2cart_gensymcrys(latvec_prim,ratred,RXYZ,NAT_CELL)
call latvec2dproj_gensymcrys(dproj,latvec_prim,rotmat,rxyz,nat_cell)

call getvol_gensymcrys(latvec_prim,vol_prim)
call getvol_gensymcrys(latvec,vol)
write(*,'(a,2f10.1)') 'volume of primitive cell and per atom ',vol_prim,vol_prim/NAT_CELL
write(*,'(a,2f10.1)') 'volume of supercell cell and per atom ',vol,vol/(NAT_CELL*t_cells(1)*t_cells(2)*t_cells(3))

cv=0.d0
cv(1,1)=dproj(1)
cv(1,2)=dproj(2)
cv(2,2)=dproj(3)
cv(1,3)=dproj(4)
cv(2,3)=dproj(5)
cv(3,3)=dproj(6)
do i=1,nat_cell
if(KINDS(i)==1) sat(i)=trim(stypat(1))
if(KINDS(i)==2) sat(i)=trim(stypat(2))
if(KINDS(i)==3) sat(i)=trim(stypat(3))
enddo
!open(unit=22,file="prim_out.ascii")
!write(22,*) NAT_CELL,LATSGP
!write(22,*) dproj(1:3)
!write(22,*) dproj(4:6)
!do i=1,nat_cell
!if(KINDS(i)==1) write(22,*)RXYZ(:,i)," ",trim(stypat(1))
!if(KINDS(i)==2) write(22,*)RXYZ(:,i)," ",trim(stypat(2))
!if(KINDS(i)==3) write(22,*)RXYZ(:,i)," ",trim(stypat(3))
!enddo
!write(*,*) "primitive cell written"
!close(22)
!
!
!call latvec2dproj_gensymcrys(dproj,latvec,rotmat,tmp,1)
!open(unit=22,file="lat_out.ascii")
!write(22,*) NAT_CELL*ncells,LATSGP
!write(22,*) dproj(1:3)
!write(22,*) dproj(4:6)
!do k=0,t_cells(1)-1
!do l=0,t_cells(2)-1
!do m=0,t_cells(3)-1
!     do i=1,nat_cell
!!************************************************************************************************************
!!This part will introduce a small noise onto the atomic positions to avoid problems due to perfect symmetry
!       call random_number(tmp)
!       tmp=(tmp-0.5d0)*1.d-6     
!!      tmp=0.d0       
!!************************************************************************************************************
!     if(KINDS(i)==1) write(22,*)RXYZ(:,i)+tmp(:)+latvec_prim(:,1)*k+latvec_prim(:,2)*l+latvec_prim(:,3)*m," ",trim(stypat(1))
!     if(KINDS(i)==2) write(22,*)RXYZ(:,i)+tmp(:)+latvec_prim(:,1)*k+latvec_prim(:,2)*l+latvec_prim(:,3)*m," ",trim(stypat(2))
!     if(KINDS(i)==3) write(22,*)RXYZ(:,i)+tmp(:)+latvec_prim(:,1)*k+latvec_prim(:,2)*l+latvec_prim(:,3)*m," ",trim(stypat(3))
!     enddo
!enddo
!enddo
!enddo
!
!
!write(*,*) "simulation supercell written"
!
!close(22)


end subroutine gensymcrys_single

!********************************************************************************


subroutine get_conv2prim(type_cell,matrix)
implicit none
character(1)::type_cell
real(8):: matrix(3,3)
if(type_cell=="A") then
matrix(1,:)=(/ 1.0d0, 0.0d0, 0.0d0/)
matrix(2,:)=(/ 0.0d0, 0.5d0,-0.5d0/)
matrix(3,:)=(/ 0.0d0, 0.5d0, 0.5d0/)
elseif(type_cell=="B") then
matrix(1,:)=(/ 0.5d0, 0.0d0, 0.5d0/)
matrix(2,:)=(/ 0.0d0, 1.0d0, 0.0d0/)
matrix(3,:)=(/-0.5d0, 0.0d0, 0.5d0/)
elseif(type_cell=="C") then
matrix(1,:)=(/ 0.5d0,-0.5d0, 0.0d0/)
matrix(2,:)=(/ 0.5d0, 0.5d0, 0.0d0/)
matrix(3,:)=(/ 0.0d0, 0.0d0, 1.0d0/)
elseif(type_cell=="F") then
matrix(1,:)=(/ 0.0d0, 0.5d0, 0.5d0/)
matrix(2,:)=(/ 0.5d0, 0.0d0, 0.5d0/)
matrix(3,:)=(/ 0.5d0, 0.5d0, 0.0d0/)
elseif(type_cell=="I") then
matrix(1,:)=(/-0.5d0, 0.5d0, 0.5d0/)
matrix(2,:)=(/ 0.5d0,-0.5d0, 0.5d0/)
matrix(3,:)=(/ 0.5d0, 0.5d0,-0.5d0/)
elseif(type_cell=="R") then
matrix(1,:)=(/ 2.0d0/3.0d0,-1.0d0/3.0d0,-1.0d0/3.0d0/)
matrix(2,:)=(/ 1.0d0/3.0d0, 1.0d0/3.0d0,-2.0d0/3.0d0/)
matrix(3,:)=(/ 1.0d0/3.0d0, 1.0d0/3.0d0, 1.0d0/3.0d0/)
else
stop "Wrong Type"
endif
end subroutine

!********************************************************************************

subroutine optcell(n_cells,cells)
implicit none
integer:: n_cells, cells(3), meancells,m,k,l,t_cells,m_cells,s_cells,p_cells

m_cells=int(n_cells**(1.d0/3.d0))
write(*,*) "Number of cells per dim in average",n_cells**(1.d0/3.d0),m_cells
cells(1)=N_CELLS
cells(2)=1
cells(3)=1
t_cells=1000000000
write(*,*) "Target number of cells:", n_cells
do m=1,n_cells
do k=1,m
do l=1,k
p_cells=m*k*l
s_cells=(m-k)**2+(k-l)**2+(l-m)**2
if(p_cells==n_cells)then
  if(s_cells.lt.t_cells) then
    t_cells=s_cells
    write(*,'(a,4(1x,i4))') "      #cells, x,y,z", s_cells,m,k,l
    cells(1)=m
    cells(2)=k
    cells(3)=l
  endif
endif
enddo
enddo
enddo
write(*,'(a,4(1x,i4))') "Found optimal cell:", m_cells, cells
end subroutine

!********************************************************************************
subroutine dist_dim(latvec,dout)
!This subroutine will return the distances between the planes of the unit cell in latvec 
implicit none
real*8 :: latvec(3,3),nvec(3,3),point(3),point0(3),dist(3),eps,dd
integer:: i
real(8):: dout(3)
! eps=1.d-6
dout=0.d0
call nveclatvec(latvec,nvec)
point0=(/0.d0,0.d0,0.d0/)
do i=1,3
call dist2plane_gensymcrys(latvec(:,mod(i+1,3)+1),nvec(:,i),point0,dist(i))
! write(*,*) "cut",i,cut, dist
enddo
dout(1)=dist(2)
dout(2)=dist(3)
dout(3)=dist(1)
end subroutine

!*********************************************************************************
subroutine read_coeffs_gensymcrys(parini,nsym_tot,rsym_all,cryssys_all,brav_all,nsymp_all,nsym_all,ind_rsym_all)
    use mod_parini, only: typ_parini
    implicit none
    type(typ_parini), intent(in):: parini
    integer, intent(in):: nsym_tot
    real(8), intent(out):: rsym_all(4,4,nsym_tot)
    integer, intent(out):: cryssys_all(230), brav_all(230)
    integer, intent(out):: nsymp_all(230), nsym_all(230), ind_rsym_all(230)
    !local variables
    integer:: ios, ispg, ii
    !real(8):: tt1, tt2, tt3, tt4
    character(256):: datafile
    if(trim(parini%datafilesdir)=="DATAFILESDIR") then
        write(*,*) 'ERROR: perhaps makefile did not succeed to set datafilesdir,'
        write(*,*) '       you can set it in flame_in.yaml in block main.'
        stop
    else
        datafile=trim(parini%datafilesdir)//"/coeffs_gensymcrys.dat"
    endif
    !write(*,*) trim(datafile)
    open(unit=11,file=trim(datafile),status='old',iostat=ios)
    if(ios/=0) then
        write(*,'(a)') 'ERROR: failure openning coeffs_gensymcrys.dat file.'
        stop
    endif
    read(11,*)
    do ispg=1,230
        if(ispg==1) then
            ind_rsym_all(ispg)=1
        else
            ind_rsym_all(ispg)=ind_rsym_all(ispg-1)+nsym_all(ispg-1)
        endif
        read(11,*) nsym_all(ispg),nsymp_all(ispg),cryssys_all(ispg),brav_all(ispg)
    enddo
    if(nsym_tot/=(ind_rsym_all(230)+nsym_all(230)-1)) then
        write(*,'(a)') 'ERROR: nsym_tot not equal to ind_rsym_all(230)+nsym_all(230)-1 !'
        write(*,'(a)') '       perhaps coeffs_gensymcrys.dat is modified.'
        stop
    endif
    read(11,*)
    do ii=1,nsym_tot
        read(11,*) rsym_all(1,1,ii),rsym_all(1,2,ii),rsym_all(1,3,ii),rsym_all(1,4,ii)
        read(11,*) rsym_all(2,1,ii),rsym_all(2,2,ii),rsym_all(2,3,ii),rsym_all(2,4,ii)
        read(11,*) rsym_all(3,1,ii),rsym_all(3,2,ii),rsym_all(3,3,ii),rsym_all(3,4,ii)
        read(11,*) rsym_all(4,1,ii),rsym_all(4,2,ii),rsym_all(4,3,ii),rsym_all(4,4,ii)
    enddo
    close(11)
end subroutine read_coeffs_gensymcrys
