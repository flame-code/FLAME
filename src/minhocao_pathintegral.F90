subroutine pathintegral(parini,parres,latvec,xred)
 use global, only: units
 use defs_basis
 use interface_code
! Main program to test potential subroutines
!       use parameters 
 use mod_parini, only: typ_parini
       implicit none
       type(typ_parini), intent(in):: parini
       type(typ_parini), intent(inout):: parres
       integer count1,count2,count_rate,count_max,lwork,info,nint,i,iat,idispl,irep,iprec,istr,jstr,ilat,jlat
! nat: number of atoms (e.g. vacancy)
!       parameter(nint=2*32*1+1)
       real(8):: rxyz0(3,parini%nat),fxyz(3,parini%nat),displ(3,parini%nat)
       real(8),allocatable:: work(:),simpson(:)
       real(8):: evals(3),s2(3,3),dmat(3,3),dproj(6),rotmat(3,3),xred(3,parini%nat)
       real(8):: stepsize_at,stepsize_lat,t1,t2,t3,path,xred_in(3,parini%nat),latvec_in(3,3),strten_in(6)
       real(8):: str_matrix(3,3),transformed(3,3),transformed_inv(3,3),fcart_in(3,parini%nat)
       real(8):: ener,vol,sumx,sumy,sumz,ener0,etot_in,tstressall
       real(8):: dlat(6),latvec(3,3),latvecinv(3,3),stress(3,3),displat(3,3),tstress(3,3)
       real(8):: stressmod(3,3),trafo(3,3),latvectrans(3,3),s2t(3,3),metric(3,3)
       real(8):: angmom(3,3),latvecold(3,3),ener_comp1,time,dener,ener_comp2,count
       character(2)::atom
       character(40):: filename
       character(4):: fn4
       logical:: getwfk,atoms,cell
       open(unit=88,file="pathparams.in")
       read(88,*) atoms,cell !Read which degrees of freedom should be tested, usually only one of them
       read(88,*) stepsize_at,stepsize_lat !Stepsize of atomic and cell degrees of freedom
       read(88,*) nint !Number of steps along the path 
       if(.not.atoms) stepsize_at=0.d0
       if(.not.cell) stepsize_lat=0.d0
       if(atoms.and.abs(stepsize_at).le.1.d-12) stop "Increase stepsize for atoms"
       if(cell.and.abs(stepsize_lat).le.1.d-12) stop "Increase stepsize for lattice"
       allocate(simpson(nint)) 

       write(*,*) '# Testing ',trim(parini%potential_potential),' potential'
       if (mod(nint,2).ne.1) stop '# nint has to be odd'
       simpson(1)=1.d0/3.d0
       simpson(2)=4.d0/3.d0
       do i=3,nint-2,2
       simpson(i)=2.d0/3.d0
       simpson(i+1)=4.d0/3.d0
       enddo
       simpson(nint)=1.d0/3.d0

       count=0.d0
       call rxyz_int2cart(latvec,xred,rxyz0,parini%nat)
! create random displacements (use sin() instead of rand())
!       stepsize=1.d-3
       call random_number(displ)
       displ=(displ-5.d-1)*2.d0
       displ=displ*stepsize_at
!       do iat=1,nat
!       displ(1,iat)=stepsize*sin(iat+.2d0)
!       displ(2,iat)=stepsize*sin(iat+.4d0)
!       displ(3,iat)=stepsize*sin(iat+.7d0)
!!       displ(1,iat)=stepsize*abs(sin(iat+.2d0))
!!       displ(2,iat)=stepsize*abs(sin(iat+.4d0))
!!       displ(3,iat)=stepsize*abs(sin(iat+.7d0))
!       enddo
       call random_number(displat)
       displat=(displat-5.d-1)*2.d0
       displat=displat*stepsize_lat   


!      displat=0.d0   
!      displ=0.d0

! calculate energy at equilibrium geometry (diamond structure)
! and at two additional points along the displacements
      call cpu_time(t1)
      call system_clock(count1,count_rate,count_max)
       path=0.d0
       do 1000,idispl=1,nint
       do 10,irep=1,1
       call rxyz_cart2int(latvec,xred_in,rxyz0,parini%nat)
       latvec_in=latvec
       iprec=1; getwfk=.true.; if(idispl==1) getwfk=.false.
       call get_energyandforces_single(parini,parres,latvec_in,xred_in,fcart_in,strten_in,etot_in,iprec,getwfk);count=count+1
       fxyz=fcart_in
       ener=etot_in
!Convert strten into "stress", which is actually the forces on the cell vectors
       str_matrix(1,1)=strten_in(1)
       str_matrix(2,2)=strten_in(2)
       str_matrix(3,3)=strten_in(3)
       str_matrix(1,2)=strten_in(6)
       str_matrix(2,1)=strten_in(6)
       str_matrix(1,3)=strten_in(5)
       str_matrix(3,1)=strten_in(5)
       str_matrix(2,3)=strten_in(4)
       str_matrix(3,2)=strten_in(4)
       call getvol(latvec,vol)
       transformed(:,1)=latvec(1,:) 
       transformed(:,2)=latvec(2,:) 
       transformed(:,3)=latvec(3,:) 
       call invertmat(transformed,transformed_inv,3)
       stress=(-vol*matmul(str_matrix,transformed_inv))
!       call energyandforces(nat,latvec,rxyz0,fxyz,stress,0.d0,ener,count)
!****************************************************************************************************************        
     write(fn4,'(i4.4)') idispl
     filename='pospath_'//fn4//'.ascii'
     call write_atomic_file_ascii(parini,filename,parini%nat,units,xred_in,latvec_in,fcart_in,strten_in,&
          &parini%char_type,parini%ntypat_global,parini%typat_global,parini%fixat,parini%fixlat,etot_in,parini%target_pressure_habohr,etot_in,0.d0)


10      continue
       if (idispl.eq.1) ener0=ener

!       if (idispl.eq.1 .or. idispl.eq.nint) then
! check whether total force vanishes
       sumx=0.d0
       sumy=0.d0
       sumz=0.d0
       do 8483,iat=1,parini%nat
           sumx=sumx+fxyz(1,iat)
           sumy=sumy+fxyz(2,iat)
8483        sumz=sumz+fxyz(3,iat)

       write(*,'(a,i5,2(x,e19.12))') ' # idispl,ener,ener/nat',idispl,ener,ener/parini%nat
       write(*,'(a,3(x,e10.3))') &
         ' # Sum of x, y and z component of forces:',sumx,sumy,sumz
!       endif


! integrate force*displacement
       t1=0.d0
       t2=0.d0
       t3=0.d0
       tstressall=0.d0
       
       do iat=1,parini%nat
       t1=t1-fxyz(1,iat)*displ(1,iat)
       t2=t2-fxyz(2,iat)*displ(2,iat)
       t3=t3-fxyz(3,iat)*displ(3,iat)
       enddo
       
       do istr=1,3
       do jstr=1,3
       tstressall=tstressall-stress(istr,jstr)*displat(istr,jstr)
       enddo
       enddo
      


       path=path+simpson(idispl)*(tstressall+t1+t2+t3)


       trafo=0.d0
      
       latvecold=latvec
       call invertmat(latvec,latvecinv,3) 
       do ilat=1,3
       do jlat=1,3
       latvec(ilat,jlat)=latvec(ilat,jlat)+displat(ilat,jlat)
       enddo
       enddo


       call updaterxyz(latvecold,latvec,rxyz0,parini%nat)


! next positions along path
       do iat=1,parini%nat
       rxyz0(1,iat)=rxyz0(1,iat)+displ(1,iat)
       rxyz0(2,iat)=rxyz0(2,iat)+displ(2,iat)
       rxyz0(3,iat)=rxyz0(3,iat)+displ(3,iat)
       enddo


1000        continue
      call cpu_time(t2)
      call system_clock(count2,count_rate,count_max)

      call latvec2dproj(dproj,latvec,rotmat,rxyz0,parini%nat)

!Calculate first energz with initial cell
      call rxyz_cart2int(latvec,xred_in,rxyz0,parini%nat)
      latvec_in=latvec
      iprec=1; getwfk=.false. 
      call get_energyandforces_single(parini,parres,latvec_in,xred_in,fcart_in,strten_in,etot_in,iprec,getwfk)
      fxyz=fcart_in
      ener_comp1=etot_in
!Convert strten into "stress", which is actually the forces on the cell vectors
       str_matrix(1,1)=strten_in(1)
       str_matrix(2,2)=strten_in(2)
       str_matrix(3,3)=strten_in(3)
       str_matrix(1,2)=strten_in(6)
       str_matrix(2,1)=strten_in(6)
       str_matrix(1,3)=strten_in(5)
       str_matrix(3,1)=strten_in(5)
       str_matrix(2,3)=strten_in(4)
       str_matrix(3,2)=strten_in(4)
       call getvol(latvec,vol)
       transformed(:,1)=latvec(1,:)
       transformed(:,2)=latvec(2,:) 
       transformed(:,3)=latvec(3,:)
       call invertmat(transformed,transformed_inv,3)
       stress=(-vol*matmul(str_matrix,transformed_inv))
!      call energyandforces(nat,latvec,rxyz0,fxyz,stress,0.d0,ener_comp1,count1)
!****************************************************************************************************************   
      

      write(*,*) '# Final forces on the transformed cell'
      do iat=1,3
        write(*,'(a,3(es25.15))') " #  ",stress(:,iat)
      enddo
!Calculate angular momentum
      call cross_product(latvec(:,1),stress(:,1),angmom(:,1))
      call cross_product(latvec(:,2),stress(:,2),angmom(:,2))
      call cross_product(latvec(:,3),stress(:,3),angmom(:,3))
      write(*,*) '# Angular momentum on the cell'
      write(*,'(a,3(es25.15))') " #  ", angmom(:,1)+angmom(:,2)+angmom(:,3) 


      call latvec2dproj(dproj,latvec,rotmat,rxyz0,parini%nat)

! compare energy difference with  force*displacement to check correctness of forces
       dener=ener-ener0
       write(*,*) '# Check correctness of forces'
       write(*,*) '# Difference of total energies ',dener
       write(*,*) '# Integral force*displacement  ',path
       write(*,*) '# Difference ',path-dener
       write(*,*) '# number of force evaluations (count)',int(count)
!       write(*,*) ' #CPU time ', t2-t1
!       time=(count2-count1)/float(count_rate)
       time=(t2-t1)
       write(*,*) '# elapsed time ', time
       write(*,*) '# time/count, time/(count*nat)',time/count, time/(count*parini%nat)
!       write(*,*) 'Energy difference with diagonalized latvec:',ener_comp1-ener_comp2
       end subroutine
