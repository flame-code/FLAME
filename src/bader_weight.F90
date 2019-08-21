!*****************************************************************************************
subroutine bader_weight(parini)
    use mod_parini, only: typ_parini
    use mod_poisson_weight, only: typ_poisson
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
    integer::aab=0,edge=0,edaj=0                          !edaj,vacuum
    real(8)::vac=0.d0                                     !edaj,vacuum
    real(8)::hxx,hyy,hzz                                  !hxx=hx , hyy=hy , hzz=hz                  
    integer::nx_cube_up,nx_cube_down                    !bound 
    integer::ny_cube_up,ny_cube_down           
    integer::nz_cube_up,nz_cube_down   
    real(8)::rho_cube(3,3,3)=0                          !cube 3*3*3 for maximaztion .
    real(8)::vlenght(3)                                  !step_size
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
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    close(1)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    write(*,*)  
    ! real(8)::time1,time2,time3
    ! call cpu_time(time1)
    call cube_read_weight(filename,nat,rat,qat,poisson,vol)
    hxx=poisson%h(1);hyy=poisson%h(2);hzz=poisson%h(3)
    allocate(list_nx(poisson%ngp(1)),list_ny(poisson%ngp(2)),list_nz(poisson%ngp(3)))
    do ii=1,nat
        orat(1:3,ii) =-oat(1:3) +rat(1:3,ii)
    end do    
    !call cpu_time(time2)
    !write(*,*)"runtime= ", time2-time1,"second"  
    write(*,*)  
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !    partition initialize ................................       !allocate rho and i don't signed it
    poisson%irho = 0   !allocate rho and i signed it ,it must will change all by /=0 numbers
    poisson%krho = 0 
!matrix..............................................................................
!*****************************************************************************************
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
!*****************************************************************************************
!step_size .......................................................................
do ii=1,3
        vlenght(ii) = sqrt(sum(lat_car(:,ii)*lat_car(:,ii)))
      end do
      poisson%step_size = minval(vlenght(:))
!vacum partition  ...............................................................  
    vac=0;aab=0  
    do nx=1,poisson%ngp(1)                                     !vacum charge
    do ny=1,poisson%ngp(2)
    do nz=1,poisson%ngp(3)
        if(poisson%rho(nx,ny,nz)/vol <=1.d-3)then
            poisson%irho(nx,ny,nz)=-2
            vac=vac+poisson%rho(nx,ny,nz)*poisson%h(1)*poisson%h(2)*poisson%h(3)
            aab=aab+1
        end if
    end do   
    end do 
    end do
    write(*,*)"number of vacuum volume:",aab*poisson%h(1)*poisson%h(2)*poisson%h(3)
    write(*,*)"number of vacuum charge:",vac
    !................................................................................. 
    !partition  .........................................  
    one:do nx=1,poisson%ngp(1)                                     !loop in loop in loop for check all rho point
    two:do ny=1,poisson%ngp(2)
    three:do nz=1,poisson%ngp(3)
        d = (/nx,ny,nz/)
        if (poisson%irho(nx,ny,nz)== 0)then 
            call m_neargrad_weight (poisson,d,i_dist,car_lat)
            path_vol= poisson%irho(d(1),d(2),d(3))  
            if (path_vol == 0) then
                print*,'new max:',d,poisson%vnum
                mat(:,poisson%vnum+1)=(d*1.d0 - 1.d0) *poisson%h  
                if (poisson%vnum >= poisson%vdim) then
                    call vol_reallocate_weight(poisson,poisson%vdim*2)
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
    do nx=1,poisson%ngp(1)                                     !vacum charge
    do ny=1,poisson%ngp(2)
    do nz=1,poisson%ngp(3)
        if(poisson%rho(nx,ny,nz)/vol <=1.d-3)then
            poisson%irho(nx,ny,nz)=poisson%vnum+1
        end if
    end do   
    end do 
    end do
!................................................................................. 
    !edge refinement ....................................................................
    last_iter=.false.
    iter=0
    do
    write(*,*) "number of iteration:    ",iter+1
    stop 'call to edag_refinement is commented, you should fix it before usging it'
    !call edag_refinement (last_iter,iter,poisson,i_dist,car_lat)
    if (last_iter) exit
    iter=iter+1
    write(*,*)""
    end do
!................................................................................. 
 poisson%krho=0
 bone:do nx=1,poisson%ngp(1)                                     !only number of edage charge 
        btwo:do ny=1,poisson%ngp(2)
        bthree:do nz=1,poisson%ngp(3)
            ! boundery conditins........................................ 
            p = (/nx,ny,nz/)
              loop:do d1 = -1,1
                do d2 = -1,1
                  do d3 = -1,1
                    d = p + (/d1,d2,d3/)
                    call bounds(d(1),d(2),d(3), poisson%ngp(1),poisson%ngp(2),poisson%ngp(3))
                    if (poisson%irho(d(1),d(2),d(3)) /= poisson%irho(p(1),p(2),p(3)))then
                    edaj=edaj+1
                    poisson%krho(nx,ny,nz)=-1
                    exit loop
                    end if
                  end do
                end do
              end do loop 
        end do bthree  
        end do btwo
        end do bone   
     sum_charge=0.d0
    do nx=1,poisson%ngp(1)
    do ny=1,poisson%ngp(2)
    do nz=1,poisson%ngp(3)
        if(poisson%krho(nx,ny,nz)==-1)sum_charge=sum_charge+poisson%rho(nx,ny,nz)
    end do
    end do
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
    !..............................................
  do nx=1,poisson%ngp(1)
    do ny=1,poisson%ngp(2)
    do nz=1,poisson%ngp(3)
      ii=poisson%irho(nx,ny,nz)
      do jj=1,nat
            if(gmax(jj,ii))then 		
                poisson%irho(nx,ny,nz)=jj
            end if
        end do
    end do
    end do
    end do
    !check.............................................
  do nx=1,poisson%ngp(1)
    do ny=1,poisson%ngp(2)
    do nz=1,poisson%ngp(3)
      if(poisson%irho(nx,ny,nz)==poisson%vnum+1)cycle
      if(poisson%irho(nx,ny,nz)> nat)write(*,*)"e"
    end do
    end do
    end do
    !weight refinement ....................................................................
    allocate(poisson%weight(poisson%ngp(1),poisson%ngp(2),poisson%ngp(3)))
    do nx=1,poisson%ngp(1)
        do ny=1,poisson%ngp(2)
          do nz=1,poisson%ngp(3)
            allocate (poisson%weight(nx,ny,nz)%w(nat))
              poisson%weight(nx,ny,nz)%w(1:nat)=0.d0
          end do
        end do
      end do
     call calc_weight(poisson,p,nat)
    !..................................................  
    allocate(charge(0:poisson%vnum+1))
    charge(1:nat)=0.d0
    do nx=1,poisson%ngp(1)
    do ny=1,poisson%ngp(2)
    do nz=1,poisson%ngp(3)
        do  ii=1,nat
        charge(ii)=charge(ii)+(poisson%rho(nx,ny,nz)*poisson%weight(nx,ny,nz)%w(ii))
        end do  
    end do
    end do
    end do 
    !weight.............................................
    write(*,*)"weight method.....................................................:"
    do  ii=1,nat
        write(*,*)"number of atom electrons",ii,poisson%h(1)*poisson%h(2)*poisson%h(3)*charge(ii) / vol
    end do    
    write(*,*)"number of all atom detected:",sum(charge)*poisson%h(1)*poisson%h(2)*poisson%h(3)/vol
end subroutine bader_weight
!*****************************************************************************************
  subroutine calc_weight(poisson, p,nat)
    use mod_poisson_weight, only: typ_poisson
    implicit none
    type(typ_poisson), intent(inout):: poisson
    integer,intent(in) :: nat
    integer :: num_edge, n1, n2, n3, d1, d2, d3
    integer :: i, iter, num_change
    integer :: p(3), pn(3)
    real(8) :: sum_top, sum_bottom, length, facet_a, r_
    real(8) :: new_weight, current_weight, wn
    real(8) :: facet_area
    logical:: is_edge_weight
    logical:: is_neighbor
    do i = 1,nat
      num_edge = 0
      do n1 = 1,poisson%ngp(1)
        do n2 = 1,poisson%ngp(2)
          do n3 = 1,poisson%ngp(3)
            p = (/n1,n2,n3/)
            if (poisson%irho(p(1),p(2),p(3)) == i) then
              poisson%weight(p(1),p(2),p(3))%w(i) = 1.d0
            end if
            if (is_edge_weight(poisson,p) .and. &
              &  ((poisson%irho(p(1),p(2),p(3)) == i) .or. is_neighbor(poisson, p, i))) then
              num_edge = num_edge+1
            end if

          end do
        end do
      end do
      write(*,'(2x,a,6x,1i8)') 'volnum = ',i,'edge points:',num_edge
      num_change = 1
      iter = 0
      do while (num_change>0)
        iter = iter + 1
        num_change = 0
        do n1 = 1,poisson%ngp(1)
          do n2 = 1,poisson%ngp(2)
            loop:do n3 = 1,poisson%ngp(3)
              p = (/n1,n2,n3/)
              do d1 = -1,1
                do d2 = -1,1
                  do d3 = -1,1
                    pn = p + (/d1,d2,d3/)
                    call bounds(pn(1),pn(2),pn(3), poisson%ngp(1),poisson%ngp(2),poisson%ngp(3))
                    if (poisson%irho(pn(1),pn(2),pn(3)) == poisson%vnum+1) cycle loop
                 
                  end do
                end do
              end do 
              if (is_edge_weight(poisson,p) .and. &
                &  ((poisson%irho(p(1),p(2),p(3)) == i) .or. is_neighbor(poisson, p, i))) then 
                sum_top = 0
                sum_bottom = 0
                do d1 = -1,1
                  do d2 = -1,1
                    do d3 = -1,1
                      pn = p + (/d1,d2,d3/) !neighbor pt
                      call bounds(pn(1),pn(2),pn(3), poisson%ngp(1),poisson%ngp(2),poisson%ngp(3)) ! just in case pn is out of the boundary
                      length = poisson%step_size 
                      facet_a = facet_area(d1, d2, d3, length)
                      r_ = dim(poisson%rho( pn(1), pn(2), pn(3)), poisson%rho(p(1),p(2),p(3)))
                      wn = poisson%weight(pn(1),pn(2),pn(3))%w(i) ! neighbor weight
                      sum_top = sum_top + facet_a*length*r_*wn
                      sum_bottom = sum_bottom + facet_a*length*r_
                    end do
                  end do
                end do
                new_weight = sum_top/sum_bottom
                current_weight = poisson%weight(p(1),p(2),p(3))%w(i)
                if (abs(new_weight - current_weight) > 0.01) then
                  poisson%weight(p(1),p(2),p(3))%w(i) = new_weight
                  num_change = num_change+1
                end if
              end if
            end do loop 
          end do
        end do
        write(*,'(2x,a,6x,1i8)') 'weight change', num_change
        write(*,'(2x,a,6x,1i8)') 'iteration', iter
      end do 
    end do
 end subroutine	
!*****************************************************************************************
  real(8) function facet_area(d1,d2,d3,length)
      integer d1, d2, d3
      real(8) length
      if (abs(d1)+abs(d2)+abs(d3) == 1) then
         facet_area = length*length
      else
         facet_area = 0
      end if
      return
   end function
!*********************************************************************************
 function is_neighbor(poisson, p, vol)
    use mod_poisson_weight, only: typ_poisson
    type(typ_poisson):: poisson
    logical :: is_neighbor
    integer,dimension(3),intent(in) :: p
    integer,dimension(3) :: pt
    integer :: d1, d2, d3, volneighbor
    integer,intent(in) :: vol
    is_neighbor = .false.
    neighborloop: do d1 = -1,1
      do d2 = -1,1
        do d3 = -1,1
          pt = p + (/d1,d2,d3/)
          if (d1 == 0 .and.d2 == 0 .and.d3 == 0 ) cycle
          call bounds(pt(1),pt(2),pt(3), poisson%ngp(1),poisson%ngp(2),poisson%ngp(3))
          volneighbor = poisson%irho(pt(1),pt(2),pt(3))
          if (volneighbor == vol) then
            is_neighbor = .true.
            exit neighborloop
          end if
        end do
      end do
    end do neighborloop
    return
    end function is_neighbor
!*****************************************************************************************
  subroutine ongrid_weight(poisson,d,i_dist)
    use mod_poisson_weight, only: typ_poisson
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
  end subroutine ongrid_weight
!*****************************************************************************************
subroutine edag_refinement_weight (last_iter,iter,poisson,i_dist,car_lat)
    use mod_poisson_weight, only: typ_poisson
    implicit none
    logical, intent(out):: last_iter
    integer, intent(in):: iter
    real(8), intent(in):: i_dist(-1:1,-1:1,-1:1),car_lat(3,3)
    type(typ_poisson), intent(inout):: poisson
    integer::n_edge,n_check,n_ressign
    integer::nx,ny,nz,mx,my,mz,d(3),dp(3)
    integer::iirho,jirho, path_vol,pt(3),i
    logical:: is_edge_weight, m_point_weight
    !first stage.............
    if(iter==0)then
        n_edge=0
        do nx=1,poisson%ngp(1)                                     
        do ny=1,poisson%ngp(2)
        do nz=1,poisson%ngp(3)
            d=(/nx,ny,nz/)
            if(poisson%irho(nx,ny,nz)== poisson%vnum+1)cycle
            if(is_edge_weight(poisson,d).and.(.not.m_point_weight (poisson,d)))then      
                n_edge=n_edge+1
                poisson%irho(d(1),d(2),d(3))= -poisson%irho(d(1),d(2),d(3))  
                poisson%krho(d(1),d(2),d(3))= 0
            end if
        end do   
        end do 
        end do
    write(*,*)"number of edge point:     ", n_edge
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
                     if(.not.m_point_weight (poisson,dp))then
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
    write(*,*)"number of check point:", n_check
    end if
    !thired stage.............
        n_ressign= 0
        do nx=1,poisson%ngp(1)                                     
        do ny=1,poisson%ngp(2)
        do nz=1,poisson%ngp(3)
            d=(/nx,ny,nz/)
            iirho= poisson%irho(nx,ny,nz)
            if(iirho< 0)then
            call m_neargrad_weight (poisson,d,i_dist,car_lat)
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
    write(*,*)"number of ressign point:", n_ressign
    if(n_ressign == 0) last_iter= .true.        
    return
end subroutine edag_refinement_weight
!***************************************************************************************
  function is_edge_weight (poisson,d) result(is_edge)
    use mod_poisson_weight, only: typ_poisson
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
  end function is_edge_weight
!.................................................................................................
  function m_point_weight (poisson,d) result(m_point)
    use mod_poisson_weight, only: typ_poisson
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
  end function m_point_weight
!***************************************************************************************
subroutine cube_read_weight(filename,nat,rat,qat,poisson,vol)
    use mod_poisson_weight, only: typ_poisson
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
        do
            str=adjustl(str)
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
end do
end do
end do    
end subroutine cube_read_weight
!*****************************************************************************************
subroutine m_neargrad_weight (poisson,d,i_dist,car_lat)
    use mod_poisson_weight, only: typ_poisson
    implicit none
    type(typ_poisson), intent(inout):: poisson 
    integer, intent(inout):: d(3)
    real(8), intent(in):: i_dist(-1:1,-1:1,-1:1),car_lat(3,3)
    poisson%pnum=1
    poisson%path(poisson%pnum,:)=d
    do
        call near_grad_weight (poisson,d,i_dist,car_lat)
        if (all(d == poisson%path(poisson%pnum,:))) exit
        if (poisson%pnum >= poisson%pdim) then
            call path_reallocate_weight(poisson, poisson%pdim*2)
        end if
        poisson%pnum = poisson%pnum+1
        poisson%path(poisson%pnum,:)=d                 
   end do
end subroutine m_neargrad_weight
!*****************************************************************************************
subroutine path_reallocate_weight(poisson, allocat2)
    use mod_poisson_weight, only: typ_poisson
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
end subroutine path_reallocate_weight
!***************************************************************************************
subroutine vol_reallocate_weight(poisson, allocat2)
    use mod_poisson_weight, only: typ_poisson
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
end subroutine vol_reallocate_weight
!*****************************************************************************************
subroutine near_grad_weight (poisson,d,i_dist,car_lat)
    use mod_poisson_weight, only: typ_poisson
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
                 pm=d
                 call ongrid_weight(poisson,pm,i_dist)   
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
             call ongrid_weight(poisson,pm,i_dist)
             dr = (/0.d0,0.d0,0.d0/) 
          end if
           d=pm
          return      
end subroutine near_grad_weight
!**********************************************************************************************
