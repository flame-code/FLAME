!*****************************************************************************************
subroutine bader_ongrid(parini)
    use mod_parini, only: typ_parini
    use mod_poisson_ongrid, only: typ_poisson
    implicit none
    type(typ_parini), intent(in):: parini
    !local variables
    type(typ_poisson):: poisson
    logical::nomax=.true.                                !diraction of moment 
    integer::nx,ny,nz                                    !This point is selected
    integer::nx_now,ny_now,nz_now
    integer::nx_now_grid,ny_now_grid,nz_now_grid,nx_now_grad,ny_now_grad,nz_now_grad
    integer::ii,jj,eqg=0,neqg=0,irho_cube(3,3,3),girho
    real(8)::lattice(3,3),lat(3,3),dlat(3),dcar(3),l_dist(-1:1,-1:1,-1:1),i_dist(-1:1,-1:1,-1:1) !MATRIX
    integer::d1,d2,d3                                                  !MATRIX
    integer::aab=0,edge=0                                 !EDAJ,VACUUM
    real(8)::vac=0.d0                                     !EDAJ,VACUUM
    real(8)::hxx,hyy,hzz                                  !hxx=hx , hyy=hy , hzz=hz                  
    integer::nx_cube_up,nx_cube_down                    !bound 
    integer::ny_cube_up,ny_cube_down           
    integer::nz_cube_up,nz_cube_down   
    real(8)::rho_cube(3,3,3)=0                          !cube 3*3*3 for maximaztion .
    integer::listi(3),listj(3),listk(3)                 !rho_cube(1:3,1:3,1:3)=poisson%rho(listi,listj,listk)
    integer::loc_cube(3)                                !showed points is selected in rho_cube
    integer::n_list                                     !Number of list would be assigned
    integer,allocatable::list_nx(:),list_ny(:),list_nz(:)
    integer::ng                                         !Numer of group 
    logical,allocatable::gmax(:,:) 
    real(8),allocatable::des(:)
    real(8),allocatable::charge(:)
    real(8),allocatable::chgat(:)
    real(8)::sum_charge
    integer::nat
    integer:: ix, iy, iz
    logical:: gp_max
    real(8):: ttmax, ttmin, tt1, tt
    real(8):: oat(3), mat(3,2000)=0.d0
    real(8),allocatable::rat(:,:), qat(:), orat(:,:)
    character(len =30)::filename
    filename=trim(parini%filename_bader)
    open(unit=1, file=filename)
      read(1,*)
      read(1,*)
      read(1,*) nat,oat(1),oat(2),oat(3)
      allocate(rat(3,nat), qat(nat),orat(3,nat),des(nat),chgat(nat))
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      do ii=1,3
          read(1,*)tt,lattice(ii,1:3)
      end do
      call transposee(lattice,lat)
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    close(1)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    write(*,*)  
    call cube_read_ongrid(filename,nat,rat,qat,poisson)
    write(*,*)"----------------------------------------------------------------"
    write(*,*)"number of all atom :",sum(poisson%rho)*poisson%hx*poisson%hy*poisson%hz
    hxx=poisson%hx;hyy=poisson%hy;hzz=poisson%hz
    allocate(list_nx(poisson%ngpx),list_ny(poisson%ngpy),list_nz(poisson%ngpz))
    do ii=1,nat
        orat(1:3,ii) =-oat(1:3) +rat(1:3,ii)
    end do    
     write(*,*)"----------------------------------------------------------------"  
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !    partition Initialize ................................       !allocate rho and I don't signed it
    poisson%irho(1:poisson%ngpx,1:poisson%ngpy,1:poisson%ngpz)= -1   !allocate rho and I signed it ,It must will change all by /=0 numbers
    ng=0
!matrix..............................................................................
!*****************************************************************************************
    do d1=-1,1
      dlat(1)=real(d1,8)
      do d2=-1,1
        dlat(2)=real(d2,8)
        do d3=-1,1
          dlat(3)=real(d3,8)
          call mat_vec (lat,dlat,dcar)
          l_dist(d1,d2,d3)=sqrt(sum(dcar*dcar))
         if ((d1 == 0).and.(d2 == 0).and.(d3 == 0)) then
            i_dist(d1,d2,d3)=0.d0
          else
            i_dist(d1,d2,d3)=1.d0/l_dist(d1,d2,d3)
          end if
        end do
      end do
    end do
!*****************************************************************************************
!vacum partition  ...............................................................  
    vac=0;aab=0  
    do nx=1,poisson%ngpx                                     !vacum charge
    do ny=1,poisson%ngpy
    do nz=1,poisson%ngpz
        if(poisson%rho(nx,ny,nz)<=1.d-3)then
            poisson%irho(nx,ny,nz)=-2
            vac=vac+poisson%rho(nx,ny,nz)*poisson%hx*poisson%hy*poisson%hz
            aab=aab+1
        end if
    end do   
    end do 
    end do
    !write(*,*)"NUMBER OF VACUUM VOLUME:",aab*poisson%hx*poisson%hy*poisson%hz
    write(*,*)"NUMBER OF VACUUM CHARGE:",vac
    !................................................................................. 
    !partition  .........................................  
    one:do nx=1,poisson%ngpx                                     !loop in loop in loop for check all rho point
    two:do ny=1,poisson%ngpy
    three:do nz=1,poisson%ngpz
        signed:if (poisson%irho(nx,ny,nz)==-1)then 
            poisson%irho(nx,ny,nz)= 0
            nx_now=nx; ny_now=ny;nz_now=nz
            n_list=0
            ineer:do
                n_list=n_list+1
                list_nx(n_list)=nx_now
                list_ny(n_list)=ny_now
                list_nz(n_list)=nz_now
                ! boundery conditins........................................
                nx_cube_up=nx_now+1; nx_cube_down=nx_now-1
                ny_cube_up=ny_now+1; ny_cube_down=ny_now-1
                nz_cube_up=nz_now+1; nz_cube_down=nz_now-1
                bound:if (nx_now==1 .or. nx_now==poisson%ngpx .or. ny_now==1 .or. ny_now==poisson%ngpy &
                    .or. nz_now==1 .or. nz_now==poisson%ngpz)then       
                    if (nx_now==1) nx_cube_down=poisson%ngpx 
                    if (nx_now==poisson%ngpx) nx_cube_up=1 
                    if (ny_now==1) ny_cube_down=poisson%ngpy 
                    if (ny_now==poisson%ngpy) ny_cube_up=1 
                    if (nz_now==1) nz_cube_down=poisson%ngpz
                    if (nz_now==poisson%ngpz) nz_cube_up=1 
                end if bound
                !..........................................................         
                listi=(/nx_cube_down,nx_now,nx_cube_up/)
                listj=(/ny_cube_down,ny_now,ny_cube_up/)
                listk=(/nz_cube_down,nz_now,nz_cube_up/)
                rho_cube(1:3,1:3,1:3)=poisson%rho(listi,listj,listk)
                loc_cube=maxloc(rho_cube)
                gp_max=.false.
                if(abs(rho_cube(2,2,2))>1.d-4) then
                    gp_max=.true.
                    loop_ix: do ix=1,3
                    do iy=1,3
                    do iz=1,3
                        if(rho_cube(2 ,2 ,2 )<rho_cube(ix,iy,iz)) then
                            gp_max=.false.
                            exit loop_ix
                        endif
                    enddo
                    enddo
                    enddo loop_ix
                endif
                cube:if(gp_max) then
                    nomax=.false.
                else 
                    call ongrid_ongrid(i_dist,rho_cube,nx_now_grid,nx_cube_up,nx_cube_down,ny_now_grid,ny_cube_up&
                    ,ny_cube_down,nz_now_grid,nz_cube_up,nz_cube_down,nx_now,ny_now,nz_now)
                    nx_now=nx_now_grid
                    ny_now=ny_now_grid
                    nz_now=nz_now_grid
                    if(nx_now_grid==nx_now_grad.and.ny_now_grid==ny_now_grad.and.nz_now_grid==nz_now_grad)then
                        eqg=eqg+1
                    else
                        neqg=neqg+1
                    end if    
                end if cube
                !.............................................................         
                iner:if(nomax)then                                                        ! It be valide change
                    iner2:if(poisson%irho(nx_now,ny_now,nz_now)==-1)then                            ! New point not singed 
                        poisson%irho(nx_now,ny_now,nz_now)= 0                                  ! singed it
                    else                                                              ! new point is singed
                        do ii=1,n_list
                            poisson%irho(list_nx(ii),list_ny(ii),list_nz(ii))=&
                            poisson%irho(nx_now,ny_now,nz_now)
                        end do
                        n_list=0
                        exit
                    end if iner2
                else                                                                        ! must be max point 
                    nomax=.true.
                    ng=ng+1 
                    do ii=1,n_list
                        poisson%irho(list_nx(ii),list_ny(ii),list_nz(ii))=ng
                    end do
                    n_list=0
                    mat(1,ng)=(nx_now*1.d0 - 1.d0) *poisson%hx  
                    mat(2,ng)=(ny_now*1.d0 - 1.d0) *poisson%hy  
                    mat(3,ng)=(nz_now*1.d0 - 1.d0) *poisson%hz  
                    !write(*,*) 'MAX_POINT ',poisson%irho(nx_now,ny_now,nz_now)&
                    !,nx_now,ny_now,nz_now,poisson%rho(nx_now,ny_now,nz_now)
                    exit
                end if iner
            end do ineer
        end if signed
    end do three  
    end do two
    end do one
    !Only For Ceak....................................
    allocate(gmax(nat,ng)) 
    gmax(1:nat,1:ng)=.false. 
    do ii=1,ng
        do jj=1,nat
            des(jj)=(mat(1,ii) -orat(1,jj))*(mat(1,ii) -orat(1,jj))+(mat(2,ii) -orat(2,jj))*&
                    (mat(2,ii) -orat(2,jj))+(mat(3,ii) -orat(3,jj))*(mat(3,ii) -orat(3,jj))
        end do 
    gmax(minloc(des),ii)=.true.
    end do
!................................ 
    allocate(charge(0:ng))
    charge(1:ng)=0.d0
    do nx=1,poisson%ngpx
    do ny=1,poisson%ngpy
    do nz=1,poisson%ngpz
        charge(poisson%irho(nx,ny,nz))=charge(poisson%irho(nx,ny,nz))+poisson%rho(nx,ny,nz)
    end do
    end do
    end do
    sum_charge=0.d0
    do ii=1,ng
        sum_charge=sum_charge+poisson%hx*poisson%hy*poisson%hz*charge(ii)
    end do
    !..............................................
    chgat(1:nat)=0.d0
    do jj=1,nat
        do ii=1,ng
            if(gmax(jj,ii)) chgat(jj)=chgat(jj)+charge(ii)
        end do
    end do
     write(*,*)"----------------------------------------------------------------"
    do  ii=1,nat
        write(*,*)"NUMBER OF ATOM ELECTRONS",ii,poisson%hx*poisson%hy*poisson%hz*chgat(ii)
    end do    
    write(*,*)"----------------------------------------------------------------"  
    write(*,*)"NUMBER OF ALL ATOM DETECTED:",sum(chgat)*poisson%hx*poisson%hy*poisson%hz
    write(*,*)"----------------------------------------------------------------"
   
    !..............................................

end subroutine bader_ongrid
!*****************************************************************************************
  subroutine ongrid_ongrid(i_dist,rho_cube,nx_now_grid,nx_cube_up,nx_cube_down,ny_now_grid,ny_cube_up,ny_cube_down,nz_now_grid,nz_cube_up,nz_cube_down,nx_now,ny_now,nz_now)
    real(8),intent(in)::i_dist(-1:1,-1:1,-1:1)
    integer,intent(in)::nx_cube_up,nx_cube_down,ny_cube_up,ny_cube_down,nz_cube_up,nz_cube_down
    real(8),intent(in)::rho_cube(3,3,3)
    INTEGER,INTENT(in) ::nx_now,ny_now,nz_now
    integer,intent(out)::nx_now_grid,ny_now_grid,nz_now_grid
    REAL(8) :: rho_max,rho_t,rho_c
    integer,dimension(3) ::pt, pm
    integer :: d1,d2,d3
    nx_now_grid=nx_now
    ny_now_grid=ny_now
    nz_now_grid=nz_now
    pm = (/2,2,2/)
    rho_c = rho_cube(2,2,2)
    rho_max = rho_c
    do d1 = 1,3
      do d2 = 1,3
        do d3 = 1,3
          pt=(/d1,d2,d3/)
          rho_t = rho_cube(d1,d2,d3)
          rho_t = rho_c+(rho_t-rho_c)*i_dist(d1-2,d2-2,d3-2)
          if (rho_t > rho_max) then
            rho_max = rho_t
            pm = pt
          end if
        end do
      end do
    end do
!    case_1 ......................................... 
   nx:select case(pm(1))
                  case(1);     nx_now_grid=nx_cube_down   
                  case(3);     nx_now_grid=nx_cube_up 
    end select nx
!    case_2 ......................................... 
   ny:select case(pm(2))
                  case(1);     ny_now_grid=ny_cube_down   
                  case(3);     ny_now_grid=ny_cube_up 
    end select ny
!    case_3 ......................................... 
   nz:select case(pm(3))
                  case(1);     nz_now_grid=nz_cube_down  
                  case(3);     nz_now_grid=nz_cube_up
    end select nz
  return
  end subroutine ongrid_ongrid
!*****************************************************************************************
subroutine cube_read_ongrid(filename,nat,rat,qat,poisson)
    use mod_poisson_ongrid, only: typ_poisson
    implicit none
    character(256):: str
    character(*), intent(in):: filename
    real(8), intent(out):: rat(3,nat), qat(nat)
    type(typ_poisson), intent(out):: poisson
    real(8):: tt
    integer:: nat, iat, igpx, igpy, igpz, ind, ios, iline, iatom
    character(2):: separator
    integer:: istat
    separator=achar(32)//achar(9)
    open(unit=1358,file=trim(filename),status='old',iostat=ios)
    if(ios/=0) then;write(*,'(2a)') 'ERROR: failure openning ',trim(filename);stop;endif
    read(1358,*)
    read(1358,*)
    read(1358,*) 
    read(1358,*) poisson%ngpx,poisson%hx
    read(1358,*) poisson%ngpy,tt,poisson%hy
    read(1358,*) poisson%ngpz,tt,tt,poisson%hz
    write(*,'(2a)') 'OPEN ',filename
    do iat=1,nat
        read(1358,*) iatom,qat(iat),rat(1,iat),rat(2,iat),rat(3,iat)
    enddo
    allocate(poisson%rho(poisson%ngpx,poisson%ngpy,poisson%ngpz),stat=istat)
    allocate(poisson%irho(poisson%ngpx,poisson%ngpy,poisson%ngpz),stat=istat)
    if(istat/=0) stop 'ERROR: allocation of rho failed.'
    igpx=1
    igpy=1
    igpz=0
    iline=0
    do
        iline=iline+1
        str=''
        read(1358,'(a)') str
        do
            str=adjustl(str)
            read(str,*,iostat=ios) tt
            ind=scan(str,separator)
            if(ios==0) then
                igpz=igpz+1
                if(igpz>poisson%ngpz) then
                    igpz=1
                    igpy=igpy+1
                    if(igpy>poisson%ngpy) then
                        igpy=1
                        igpx=igpx+1
                    endif
                endif
                poisson%rho(igpx,igpy,igpz)=tt
            else
            endif
            if(ind>1) then
                str(1:ind)=''
            else
                exit
            endif
        enddo
        if(igpx==poisson%ngpx .and. igpy==poisson%ngpy .and. igpz==poisson%ngpz) then
            write(*,*)"CLOSE ",FILENAME
            exit
        endif
    enddo
    close(1358)
end subroutine cube_read_ongrid
!*****************************************************************************************
