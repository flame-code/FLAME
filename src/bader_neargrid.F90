!*****************************************************************************************
subroutine bader_neargrid(parini)
    use mod_parini, only: typ_parini
    use mod_poisson_neargrid, only: typ_poisson
    implicit none
    type(typ_parini), intent(in):: parini
    !local variables
    type(typ_poisson):: poisson
    integer::nx,ny,nz,path_vol,p(3),i                      !this point is selected
    integer::ii,jj,eqg=0,neqg=0,irho_cube(3,3,3),girho
    real(8)::lattice(3,3),dlat(3),dcar(3),l_dist(-1:1,-1:1,-1:1),i_dist(-1:1,-1:1,-1:1) !matrix
    integer::d(3),d1,d2,d3                                                  !matrix
    real(8)::car_lat(3,3),lat_car(3,3),rho_grad_d(3)
    real(8)::dr(3)=0.d0,grad_rho_l(3),grad_rho_c(3), c_grad               !near-grad
    integer::aab=0,edge=0                                 !edaj,vacuum
    real(8)::vac=0.d0                                     !edaj,vacuum
    real(8)::hxx,hyy,hzz                                  !hxx=hx , hyy=hy , hzz=hz                  
    integer::nx_cube_up,nx_cube_down                    !bound 
    integer::ny_cube_up,ny_cube_down           
    integer::nz_cube_up,nz_cube_down   
    real(8)::rho_cube(3,3,3)=0                          !cube 3*3*3 for maximaztion .
    integer,parameter::matl=1e5
    integer::listi(3),listj(3),listk(3)                 !rho_cube(1:3,1:3,1:3)=poisson%rho(listi,listj,listk)
    integer::loc_cube(3)                                !showed points is selected in rho_cube
    integer::n_list                                     !number of list would be assigned
    integer,allocatable::list_nx(:),list_ny(:),list_nz(:)
    integer::ng                                         !numer of group 
    logical::last_iter=.false.                          !edge refinement
    integer::iter=1                                     !edge refinement
    logical,allocatable::gmax(:,:) 
    real(8),allocatable::des(:)
    real(8),allocatable::charge(:)
    real(8),allocatable::chgat(:)
    real(8)::sum_charge,vol,mat_vol
    integer::nat
    integer:: ix, iy, iz
    logical:: gp_max
    logical :: max_point
    real(8):: ttmax, ttmin, tt1, tt
    real(8):: oat(3), mat(3,matl)=0.d0
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
          read(1,*)poisson%ngp(ii),lattice(ii,1:3)
      end do
      call transposee(lattice,lat_car)
      call inverse(lat_car,car_lat)
      do ii=1,3
          lattice(ii,:)=lattice(ii,:)*real(poisson%ngp(ii))
      end do
      vol=mat_vol(lattice)
     close(1)
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    call cube_read_neargrid(filename,nat,rat,qat,poisson,vol)
    write(*,*)"----------------------------------------------------------------"
    write(*,*)"number of all atom :",sum(poisson%rho)*poisson%h(1)*poisson%h(2)*poisson%h(3)/vol
    hxx=poisson%h(1);hyy=poisson%h(2);hzz=poisson%h(3)
    allocate(list_nx(poisson%ngp(1)),list_ny(poisson%ngp(2)),list_nz(poisson%ngp(3)))
    do ii=1,nat
        orat(1:3,ii) =-oat(1:3) +rat(1:3,ii)
    end do    
    write(*,*)"----------------------------------------------------------------"
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !    partition initialize ................................       !allocate rho and i don't signed it
    poisson%irho = 0   !allocate rho and i signed it ,it must will change all by /=0 numbers
    poisson%krho = 0 
!matrix..............................................................................
    do d1=-1,1
      dlat(1)=real(d1,8)
      do d2=-1,1
        dlat(2)=real(d2,8)
        do d3=-1,1
          dlat(3)=real(d3,8)
          call mat_vec (lat_car,dlat,dcar)
          l_dist(d1,d2,d3)=sqrt(sum(dcar*dcar))
         if ((d1 == 0).and.(d2 == 0).and.(d3 == 0)) then
            i_dist(d1,d2,d3)=0.d0
          else
            i_dist(d1,d2,d3)=1.d0/l_dist(d1,d2,d3)
          end if
        end do
      end do
    end do
!vacum partition  ...............................................................  
    vac=0;aab=0  
    if(trim(parini%vacuum_bader)=='yes') then
    do nx=1,poisson%ngp(1)                                     !vacum charge
    do ny=1,poisson%ngp(2)
    do nz=1,poisson%ngp(3)
        if(poisson%rho(nx,ny,nz)/vol <=1.d-3)then
            poisson%irho(nx,ny,nz)=-2
            vac=vac+poisson%rho(nx,ny,nz)
            aab=aab+1
        end if
    end do   
    end do 
    end do
    endif
    !write(*,*)"number of vacuum volume:",aab*poisson%h(1)*poisson%h(2)*poisson%h(3)
    write(*,*)"number of vacuum charge:",vac*poisson%h(1)*poisson%h(2)*poisson%h(3)/vol
    !................................................................................. 
    !partition  .........................................  
    one:do nx=1,poisson%ngp(1)                                     !loop in loop in loop for check all rho point
    two:do ny=1,poisson%ngp(2)
    three:do nz=1,poisson%ngp(3)
        d = (/nx,ny,nz/)
        if (poisson%irho(nx,ny,nz)== 0)then 
            call m_neargrad (poisson,d,i_dist,car_lat)
            path_vol= poisson%irho(d(1),d(2),d(3))  
            if (path_vol == 0) then
                !print*,'new max:',d,poisson%vnum
                mat(:,poisson%vnum+1)=(d*1.d0 - 1.d0) *poisson%h  
                if (poisson%vnum >= poisson%vdim) then
                    call vol_reallocate(poisson,poisson%vdim*2)
                end if
                poisson%vnum = poisson%vnum+1
                path_vol = poisson%vnum
                poisson%volu(poisson%vnum,:) = real(d,8)
            end if
            do i=1,poisson%pnum
                p = (/poisson%path(i,1),poisson%path(i,2),poisson%path(i,3)/)
                if(poisson%irho(p(1),p(2),p(3)) /= -2) then
                    poisson%irho(p(1),p(2),p(3)) = path_vol
                end if
                poisson%krho(p(1),p(2),p(3)) = 0
            end do 
        end if
    end do three  
    end do two
    end do one
!vacum partition  ...............................................................  
    if(trim(parini%vacuum_bader)=='yes') then
    do nx=1,poisson%ngp(1)                                     !vacum charge
    do ny=1,poisson%ngp(2)
    do nz=1,poisson%ngp(3)
        if(poisson%rho(nx,ny,nz)/vol <=1.d-3)then
            poisson%irho(nx,ny,nz)=poisson%vnum+1
        end if
    end do   
    end do 
    end do
    endif
!edge refinement ....................................................................
    last_iter=.false.
    iter=0
    do
    !write(*,*) "number of iteration:    ",iter+1
    call edag_refinement (last_iter,iter,poisson,i_dist,car_lat)
    if (last_iter) exit
    iter=iter+1
    !write(*,*)""
    end do
 !only for ceak bader volums.......................................................................
    allocate(gmax(nat,poisson%vnum)) 
    gmax(1:nat,1:poisson%vnum)=.false. 
    do ii=1,poisson%vnum
        do jj=1,nat
            des(jj)=(mat(1,ii) -orat(1,jj))*(mat(1,ii) -orat(1,jj))+(mat(2,ii) -orat(2,jj))*&
                    (mat(2,ii) -orat(2,jj))+(mat(3,ii) -orat(3,jj))*(mat(3,ii) -orat(3,jj))
        end do 
    gmax(minloc(des),ii)=.true.
    end do
    !..................................................  
    allocate(charge(0:poisson%vnum+1))
    charge(1:poisson%vnum)=0.d0
    do nx=1,poisson%ngp(1)
    do ny=1,poisson%ngp(2)
    do nz=1,poisson%ngp(3)
        charge(poisson%irho(nx,ny,nz))=charge(poisson%irho(nx,ny,nz))+poisson%rho(nx,ny,nz)
    end do
    end do
    end do
    sum_charge=0.d0
    !..............................................
    chgat(1:nat)=0.d0
    do jj=1,nat
        do ii=1,poisson%vnum
            if(gmax(jj,ii)) chgat(jj)=chgat(jj)+charge(ii)
        end do
    end do
    write(*,*)"----------------------------------------------------------------"
    do  ii=1,nat
        write(*,*)"number of atom electrons",ii,poisson%h(1)*poisson%h(2)*poisson%h(3)*chgat(ii) / vol
    end do  
    write(*,*)"----------------------------------------------------------------"  
    write(*,*)"number of all atom detected:",sum(chgat)*poisson%h(1)*poisson%h(2)*poisson%h(3)/vol
    write(*,*)"----------------------------------------------------------------"
    !..............................................
end subroutine bader_neargrid
!*****************************************************************************************
  subroutine ongrid_neargrid(poisson,d,i_dist)
    use mod_poisson_neargrid, only: typ_poisson
    implicit none
    type(typ_poisson), intent(inout):: poisson
    integer, intent(inout):: d(3)
    real(8),intent(in)::i_dist(-1:1,-1:1,-1:1)
    real(8) :: rho_max,rho_t,rho_c
    integer,dimension(3) ::pt, pm
    integer :: p1,p2,p3,ii,jj,kk
    pm = (/d(1),d(2),d(3)/)
    rho_c = poisson%rho(d(1),d(2),d(3))
    rho_max = rho_c
      !....................................................... 
        do ii= -1,1
           p1=d(1)+ii 
           do jj= -1,1
              p2=d(2)+jj
                do kk= -1,1
                   p3=d(3)+kk           
                    call bounds(p1,p2,p3,poisson%ngp(1),poisson%ngp(2),poisson%ngp(3))
                    pt=(/p1,p2,p3/)
                    rho_t = poisson%rho(p1,p2,p3)
                    rho_t = rho_c+(rho_t-rho_c)*i_dist(ii,jj,kk)
                    if (rho_t > rho_max) then
                        rho_max = rho_t
                        pm = pt
                    end if
                end do 
            end do
        end do  
        d=pm               
  return
  end subroutine ongrid_neargrid
!*****************************************************************************************
!max............................................................................................
  function max_point(rho_cube)

    real(8),intent(in) :: rho_cube(3,3,3)
    logical :: max_point

    integer :: d1, d2, d3
  
    max_point=.true. 
    do d1=1,3
      do d2=1,3
        do d3=1,3
          if(d1==0.and.d2==0.and.d3==0)cycle
          if(rho_cube(d1,d2,d3) > rho_cube(2,2,2)) then
            max_point = .false.
          end if
        end do
      end do
    end do

  return 
  end function max_point
!bound...........................................................................................
  subroutine bounds(p1,p2,p3,p_max1,p_max2,p_max3)
  implicit none
    integer,intent(in)::p_max1,p_max2,p_max3
    integer,intent(inout)::p1,p2,p3
    integer::p_max(3),p(3)
    integer :: i
    p(1)=p1;p(2)=p2;p(3)=p3
    p_max(1)=p_max1;p_max(2)=p_max2;p_max(3)=p_max3
    do i=1,3
      do
        if(p(i) > 0) exit
        p(i) = p(i) + p_max(i)
      end do
      do 
        if(p(i) <= p_max(i)) exit
        p(i) = p(i) - p_max(i)
      end do
    end do
   p1=p(1);p2=p(2);p3=p(3)
  return
  end subroutine bounds
!matrix_vector-----------------------------------------------------------------------------------!
  subroutine mat_vec(m,v,vp)
    real(8),intent(in),dimension(3,3) :: m
    real(8),intent(in),dimension(3) :: v
    real(8),intent(out),dimension(3) :: vp
    integer :: i,n
    n=size(v)
    vp=0.d0
    do i=1,n
      vp=vp+v(i)*m(:,i)
    end do
  return
  end subroutine mat_vec
!vector_matrix.......................................................................................
  subroutine vec_mat(v,m,vp)
    real(8),intent(in),dimension(3,3) :: m
    real(8),intent(in),dimension(3) :: v
    real(8),intent(out),dimension(3) :: vp
    integer :: i,n
    n=size(v)
    vp=0.d0
    do i=1,n
      vp=vp+v(i)*m(i,:)
    end do
  return
  end subroutine vec_mat
!.......................................................................................
  subroutine transposee(lattice,lat_car)
    real(8),intent(in),dimension(3,3) :: lattice
    real(8),intent(inout),dimension(3,3) :: lat_car
    integer :: i,j,n,m
    n=size(lattice,1)
    m=size(lat_car,2)
    do i=1,n
      do j=1,m
        lat_car(j,i)=lattice(i,j)
      end do
    end do
  return
  end subroutine transposee
!-----------------------------------------------------------------------------------!
  subroutine inverse(a,b)
    real(8),intent(in),dimension(3,3) :: a
    real(8),intent(out),dimension(3,3) :: b
    real(8) :: det
    integer :: i,j,it,jt
    det=0
    do i=1,3
      it=i-1
      do j=1,3
        jt=j-1
        b(j,i) = & 
        &      a(mod(it+1,3)+1,mod(jt+1,3)+1)*a(mod(it+2,3)+1,mod(jt+2,3)+1)  &
        &     -a(mod(it+1,3)+1,mod(jt+2,3)+1)*a(mod(it+2,3)+1,mod(jt+1,3)+1)
      end do
      det=det+a(i,1)*b(1,i)
    end do

    do i=1,3
      do j=1,3
        b(i,j)=b(i,j)/det
      end do
    end do
  return
  end subroutine inverse
!***************************************************************************************
  function mat_vol(h)

    real(8),intent(in),dimension(3,3) :: h
    real(8) :: mat_vol
    
    mat_vol = h(1,1)*(h(2,2)*h(3,3)-h(2,3)*h(3,2))  &
    &              -h(1,2)*(h(2,1)*h(3,3)-h(3,1)*h(2,3))  &
    &              +h(1,3)*(h(2,1)*h(3,2)-h(3,1)*h(2,2))

    mat_vol = abs(mat_vol)

  return
  end function mat_vol
!***************************************************************************************
subroutine edag_refinement (last_iter,iter,poisson,i_dist,car_lat)
    use mod_poisson_neargrid, only: typ_poisson
    implicit none
    logical, intent(out):: last_iter
    integer, intent(in):: iter
    real(8), intent(in):: i_dist(-1:1,-1:1,-1:1),car_lat(3,3)
    type(typ_poisson), intent(inout):: poisson
    integer::n_edge,n_check,n_ressign
    integer::nx,ny,nz,mx,my,mz,d(3),dp(3)
    integer::iirho,jirho, path_vol,pt(3),i
    logical:: is_edge_neargrid, m_point
    !first stage.............
    if(iter==0)then
        n_edge=0
        do nx=1,poisson%ngp(1)                                     
        do ny=1,poisson%ngp(2)
        do nz=1,poisson%ngp(3)
            d=(/nx,ny,nz/)
            if(poisson%irho(nx,ny,nz)== poisson%vnum+1)cycle
            if(is_edge_neargrid(poisson,d).and.(.not.m_point (poisson,d)))then      
                n_edge=n_edge+1
                poisson%irho(d(1),d(2),d(3))= -poisson%irho(d(1),d(2),d(3))  
                poisson%krho(d(1),d(2),d(3))= 0
            end if
        end do   
        end do 
        end do
    !write(*,*)"number of edge point:     ", n_edge
    end if
    !second stage.............
    if(iter > 0)then
        n_check= 0
        do nx=1,poisson%ngp(1)                                     
        do ny=1,poisson%ngp(2)
        do nz=1,poisson%ngp(3)
            d=(/nx,ny,nz/)
            if(poisson%irho(nx,ny,nz)== poisson%vnum+1)cycle
            if(poisson%irho(nx,ny,nz) < 0 .and. poisson%krho(nx,ny,nz) /= -1 )then
                 do mx=-1,1
                 do my=-1,1
                 do mz=-1,1
                     dp= d+ (/mx,my,mz/)
                     call bounds(dp(1),dp(2),dp(3),poisson%ngp(1),poisson%ngp(2),poisson%ngp(3))  
                     if(poisson%irho(dp(1),dp(2),dp(3))== poisson%vnum+1)cycle
                     if(.not.m_point (poisson,dp))then
                         if(poisson%irho(dp(1),dp(2),dp(3)) > 0 )then
                             poisson%irho(dp(1),dp(2),dp(3))= -poisson%irho(dp(1),dp(2),dp(3))   
                             poisson%krho(dp(1),dp(2),dp(3))= -1
                             n_check= n_check+1
                         else if(poisson%irho(dp(1),dp(2),dp(3)) < 0 .and. poisson%krho(dp(1),dp(2),dp(3))== 0)then
                             poisson%krho(dp(1),dp(2),dp(3))= -2
                             n_check= n_check+1
                         end if
                    end if
                 end do
                 end do
                 end do
                 n_check= n_check-1
            if(poisson%krho(dp(1),dp(2),dp(3)) /= -2) poisson%irho(d(1),d(2),d(3))= abs(poisson%irho(d(1),d(2),d(3))) 
            end if                
        end do   
        end do 
        end do
    !write(*,*)"number of check point:", n_check
    end if  
    !thired stage.............
        n_ressign= 0
        do nx=1,poisson%ngp(1)                                     
        do ny=1,poisson%ngp(2)
        do nz=1,poisson%ngp(3)
            d=(/nx,ny,nz/)
            iirho= poisson%irho(nx,ny,nz)
            if(iirho< 0)then
            call m_neargrad (poisson,d,i_dist,car_lat)
            path_vol=poisson%irho(d(1),d(2),d(3))
            poisson%irho(nx,ny,nz)= path_vol
                if(abs(iirho) /= path_vol) then
                    n_ressign=n_ressign + 1
                    poisson%irho(nx,ny,nz)=-path_vol
                end if
                do i = 1,poisson%pnum
                    pt = (/poisson%path(i,1),poisson%path(i,2),poisson%path(i,3)/)
                       poisson%krho(pt(1),pt(2),pt(3)) = 0
                end do 
            end if                     
        end do   
        end do 
        end do
    !write(*,*)"number of ressign point:", n_ressign
    if(n_ressign == 0) last_iter= .true.        
    return
  end subroutine edag_refinement

!***************************************************************************************
  function is_edge_neargrid (poisson,d) result(is_edge)
    use mod_poisson_neargrid, only: typ_poisson
    type(typ_poisson):: poisson
    logical :: is_edge

    integer,intent(in) :: d(3)
    integer ::pd(3)
    integer :: d1,d2,d3,num,nbr

    num = poisson%irho(d(1),d(2),d(3))
    is_edge = .false.
    nei: do d1 = -1,1
      do d2 = -1,1
        do d3 = -1,1
          pd = d + (/d1,d2,d3/)
          call bounds(pd(1),pd(2),pd(3),poisson%ngp(1),poisson%ngp(2),poisson%ngp(3))
          nbr = poisson%irho(pd(1),pd(2),pd(3))
          if (abs(nbr) /= abs(num)) then
            is_edge = .true.
            exit nei  
          end if
        end do
      end do
    end do nei
  end function is_edge_neargrid
!.................................................................................................
  function m_point (poisson,d)
    use mod_poisson_neargrid, only: typ_poisson
    type(typ_poisson):: poisson
    logical :: m_point
    integer,intent(in) :: d(3)
    integer :: p1, p2 ,p3 ,d1 ,d2, d3,p_max1,p_max2,p_max3	
    m_point=.true. 
    do d1=-1,1
      do d2=-1,1
        do d3=-1,1
          p1=d(1)+d1; p2=d(2)+d2; p3=d(3)+d3
          call bounds(p1,p2,p3,poisson%ngp(1),poisson%ngp(2),poisson%ngp(3))
          if(poisson%rho(p1,p2,p3) > poisson%rho(d(1),d(2),d(3))) then
            m_point = .false.
          end if
        end do
      end do
    end do

  return 
  end function 
!***************************************************************************************
subroutine cube_read_neargrid(filename,nat,rat,qat,poisson,vol)
    use mod_poisson_neargrid, only: typ_poisson
    implicit none
    character(256):: str
    character(*), intent(in):: filename
    real(8), intent(out):: rat(3,nat), qat(nat)
    real(8), intent(in):: vol
    type(typ_poisson), intent(out):: poisson
    real(8):: tt
    integer:: nat, iat, igpx, igpy, igpz, ind, ios, iline, iatom
    character(2):: separator
    integer:: istat,n1,n2,n3
    separator=achar(32)//achar(9)
    open(unit=1358,file=trim(filename),status='old',iostat=ios)
    if(ios/=0) then;write(*,'(2a)') 'error: failure openning ',trim(filename);stop;endif
    read(1358,*)
    read(1358,*)
!    read(1358,*) nat
    read(1358,*) 
    read(1358,*) poisson%ngp(1),poisson%h(1)
    read(1358,*) poisson%ngp(2),tt,poisson%h(2)
    read(1358,*) poisson%ngp(3),tt,tt,poisson%h(3)
    write(*,'(2a)') 'open ',filename
    do iat=1,nat
        read(1358,*) iatom,qat(iat),rat(1,iat),rat(2,iat),rat(3,iat)
    enddo
    allocate(poisson%rho(poisson%ngp(1),poisson%ngp(2),poisson%ngp(3)),stat=istat)
    allocate(poisson%irho(poisson%ngp(1),poisson%ngp(2),poisson%ngp(3)),stat=istat)
    allocate(poisson%krho(poisson%ngp(1),poisson%ngp(2),poisson%ngp(3)),stat=istat)
    poisson%vdim=64; poisson%pdim=64; poisson%vnum=0; poisson%nvol=0; poisson%pnum=0
    allocate(poisson%volu(poisson%vdim,3))
    allocate(poisson%path(poisson%pdim,3))

    if(istat/=0) stop 'error: allocation of rho failed.'
    igpx=1
    igpy=1
    igpz=0
    iline=0
    do
        iline=iline+1
        str=''
        read(1358,'(a)') str
        !str='  1.20645e-04   1.45834e-04   1.76737e-04   2.14873e-04   2.61824e-04   3.19073e-04'
        do
            str=adjustl(str)
            !write(*,'(a)') trim(str)
            read(str,*,iostat=ios) tt
            ind=scan(str,separator)
            if(ios==0) then
                igpz=igpz+1
                if(igpz>poisson%ngp(3)) then
                    igpz=1
                    igpy=igpy+1
                    if(igpy>poisson%ngp(2)) then
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
        if(igpx==poisson%ngp(1) .and. igpy==poisson%ngp(2) .and. igpz==poisson%ngp(3)) then
            write(*,*)"close ",filename
            exit
        endif
    enddo
    !print*, " elec. and ionic charges: " ,  sum(rho)*hx*hy*hz , sum(q)
    close(1358)
      do n1=1,poisson%ngp(1)
      do n2=1,poisson%ngp(2)
        do n3=1,poisson%ngp(3)
          poisson%rho(n1,n2,n3)=vol*poisson%rho(n1,n2,n3)
        end do
      end do
    end do
do n1=1,poisson%ngp(1)
do n2=1,poisson%ngp(2)
do n3=1,poisson%ngp(3)
 !write(31,*)poisson%rho(n1,n2,n3)
end do
end do
end do    
end subroutine cube_read_neargrid
!*****************************************************************************************
subroutine m_neargrad (poisson,d,i_dist,car_lat)
    use mod_poisson_neargrid, only: typ_poisson
    implicit none
    type(typ_poisson), intent(inout):: poisson 
    integer, intent(inout):: d(3)
    real(8), intent(in):: i_dist(-1:1,-1:1,-1:1),car_lat(3,3)
    poisson%pnum=1
    poisson%path(poisson%pnum,:)=d
    do
        call near_grad (poisson,d,i_dist,car_lat)
        if (all(d == poisson%path(poisson%pnum,:))) exit
        if (poisson%pnum >= poisson%pdim) then
            call path_reallocate(poisson, poisson%pdim*2)
        end if
        poisson%pnum = poisson%pnum+1
        poisson%path(poisson%pnum,:)=d                 
   end do
end subroutine m_neargrad
!*****************************************************************************************
subroutine path_reallocate(poisson, allocat2)
    use mod_poisson_neargrid, only: typ_poisson
    implicit none
    type(typ_poisson), intent(inout):: poisson 
    integer :: allocat2
    integer,allocatable :: tpath(:,:)

    allocate(tpath(poisson%pdim,3))
    tpath = poisson%path
    deallocate(poisson%path)
    poisson%pdim = allocat2
    allocate(poisson%path(poisson%pdim,3))
    poisson%path(1:poisson%pnum,:) = tpath(1:poisson%pnum,:)
    deallocate(tpath)

end subroutine path_reallocate
!***************************************************************************************
subroutine vol_reallocate(poisson, allocat2)
    use mod_poisson_neargrid, only: typ_poisson
    implicit none
    type(typ_poisson), intent(inout):: poisson 
    integer :: allocat2
    real(8),allocatable :: tpath(:,:)
    allocate(tpath(poisson%vdim,3))
    tpath = poisson%volu
    deallocate(poisson%volu)
    poisson%vdim = allocat2
    allocate(poisson%volu(poisson%vdim,3))
    poisson%volu(1:poisson%vnum,:) = tpath(1:poisson%vnum,:)
    deallocate(tpath)
end subroutine vol_reallocate
!*****************************************************************************************
subroutine near_grad (poisson,d,i_dist,car_lat)
    use mod_poisson_neargrid, only: typ_poisson
    implicit none
    type(typ_poisson), intent(inout):: poisson
    integer, intent(inout):: d(3)
    real(8), intent(in):: i_dist(-1:1,-1:1,-1:1),car_lat(3,3)
    integer::p1,p2,p3,ii,jj,kk,pm(3)
    real(8):: rho_cube(3,3,3),grad_rho_l(3)
    real(8)::dr(3)=(/0.d0,0.d0,0.d0/), rho_grad_d(3),grad_rho_c(3), c_grad
    logical:: max_point
       save dr
       if (poisson%pnum == 1) then
           dr = (/0.d0,0.d0,0.d0/)
        end if
        !....................................................... 
        do ii= -1,1
           p1=d(1)+ii 
           do jj= -1,1
              p2=d(2)+jj
                do kk= -1,1
                   p3=d(3)+kk           
                    call bounds(p1,p2,p3,poisson%ngp(1),poisson%ngp(2),poisson%ngp(3))
                    rho_cube(ii+2,jj+2,kk+2)=poisson%rho(p1,p2,p3)
                end do 
            end do
        end do                       
        !grad initialize........................................
        grad_rho_l(1) =.5d0 *(rho_cube(3,2,2) -rho_cube(1,2,2))
        grad_rho_l(2) =.5d0 *(rho_cube(2,3,2) -rho_cube(2,1,2))
        grad_rho_l(3) =.5d0 *(rho_cube(2,2,3) -rho_cube(2,2,1))
        !grad initialize........................................
        if(rho_cube(3,2,2)<rho_cube(2,2,2).and.rho_cube(1,2,2)<rho_cube(2,2,2)) grad_rho_l(1) =0.d0
        if(rho_cube(2,3,2)<rho_cube(2,2,2).and.rho_cube(2,1,2)<rho_cube(2,2,2)) grad_rho_l(2) =0.d0
        if(rho_cube(2,2,3)<rho_cube(2,2,2).and.rho_cube(2,2,1)<rho_cube(2,2,2)) grad_rho_l(3) =0.d0
        !.......................................................
        call vec_mat(grad_rho_l,car_lat,grad_rho_c) 
        call mat_vec(car_lat,grad_rho_c,rho_grad_d)
        mic:if(maxval(abs(rho_grad_d)) < 1e-30)then
             if (max_point(rho_cube))then
                 dr = (/0.d0,0.d0,0.d0/)
                 return
             else
                 !on-grid method      
                 pm=d
                 call ongrid_neargrid(poisson,pm,i_dist)   
                 dr = (/0.d0,0.d0,0.d0/)
             end if
         else
             c_grad=1.d0/maxval(abs(rho_grad_d))
             pm =d +anint(c_grad *rho_grad_d)
             dr = dr + (c_grad *rho_grad_d) - anint(c_grad *rho_grad_d)    
             pm =pm +anint(dr)
             dr = dr - anint(dr)
         end if mic
         poisson%krho(d(1),d(2),d(3))=1
         call bounds(pm(1),pm(2),pm(3),poisson%ngp(1),poisson%ngp(2),poisson%ngp(3))  
 
         if(poisson%krho(pm(1),pm(2),pm(3))==1 )then
             pm=d
             call ongrid_neargrid(poisson,pm,i_dist)
             dr = (/0.d0,0.d0,0.d0/) 
          end if
           d=pm
          return      
end subroutine near_grad
!**********************************************************************************************
