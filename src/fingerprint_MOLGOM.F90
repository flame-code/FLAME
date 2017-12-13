!module parameter_molgom
!   implicit none
!   save
!
!   logical, parameter :: write_files = .false.
!   logical, parameter :: clustering = .false.
!   integer, parameter :: cluster_number = 20
!   integer, parameter :: nat=20
!   integer, parameter :: ntypat=4
!   integer, parameter :: principleev = 6
!   integer, parameter :: lseg=1
!   integer, parameter :: nconf=177
!   integer, parameter :: molecules=4
!   integer, parameter :: fp_18_expaparameter = 4
!   integer, parameter :: fp_18_nex_cutoff = 3
!   integer, parameter :: molecules_sphere = 50
!   real*8,  parameter :: fp_18_width_cutoff = 1.d0
!   real*8,  parameter :: width_overlap = 1.d0
!   real*8,  parameter :: large_vanradius = 1.7d0/0.52917720859d0
!   character(len=100), parameter :: filetype = 'ascii'
!   character(len=100) :: f1, f2,f3, f4, f5, f6, f7, f8, f9, f10, f11, f12
!
!end module parameter_molgom
!
!
!program molecular_fingerprint
!   use parameter_molgom
!   implicit none
!   !implicit real*8 (a-h,o-z)
!
!   integer :: iconf,imol,i,jconf,jmol,j,count,curpos,IOstatus,startconf,n
!   real*8 :: endiff,sum,max_distance
!   real*8 :: alat(3,3),cost(molecules,molecules)
!   integer:: iassign(molecules)
!   real*8 :: tt,dtt,max
!!   dimension alat(3, 3),iassign(molecules),cost(molecules,molecules)
!   real*8 :: rxyz(3, nat*molecules),rvan(nat*molecules),rcov(nat*molecules)
!   real*8 :: energy(nconf), c1_energy(3,nconf)
!   real*8 :: fpsall(lseg*molecules_sphere*principleev,molecules,nconf)
!   character(len=2), dimension(nat*molecules) :: finalchar
!   character(len=10) :: char
!   real*8, dimension(nconf,nconf) :: cluster_arr
!   logical, dimension(nconf) :: surviving
!   logical :: ex
!
!
!
!   if (lseg.eq.1) then
!      write(*,*) 'S orbitals used'
!   else if (lseg.eq.4) then
!      write(*,*) 'S+P orbitals used'
!   else
!      write(*,*) 'wrong lseg'
!   endif
!
!   write(*,*) "It is assumed that only ",ntypat," types of atoms can be found in the system"
!
!   fpsall=0.d0
!   jconf=0
!
!   !read first the fingerprints.txt if existing. this file saves the already acquired fingerprints! make sure the dimesions are
!   !as before (i.e. lseg, molecules_sphere, principleev)
!
!   inquire(file='fingerprints.txt',exist=ex)
!   if (ex .eqv. .true.) then
!      open (17,file='fingerprints.txt')
!      do iconf=1, nconf
!         read(17,*,IOSTAT=IOstatus) char, count
!
!         if (IOstatus/=0) then
!
!            write(*,*) jconf ,"moleculare fingerprints alredy found from the fingerprints.txt"
!            exit
!         endif
!
!         do n=1, lseg*molecules_sphere*principleev
!            read(17,*) (fpsall(n,j,iconf),j=1,molecules)
!         enddo
!
!         read(17,*,IOSTAT=IOstatus)
!
!         if (IOstatus/=0) then
!
!            write(*,*) jconf ,"moleculare fingerprints alredy found from the fingerprints.txt"
!            exit
!         endif
!         jconf=iconf
!      enddo
!      close(17)
!   endif
!
!   open (17,file='fingerprints.txt')
!
!   startconf=jconf+1
!
!
!   do i=1, startconf-1
!      write(17,*) 'conf: ', i
!      do n=1,lseg*molecules_sphere*principleev
!         write(17,*) (fpsall(n,j,i),j=1,molecules)
!      enddo
!      write(17,*)
!   enddo
!
!   write(f1, '(I5.5)') startconf
!   f2 = 'poslocm_'//trim(adjustl(f1))//'.ascii'
!   write(*,*) "fingerprint method starts with:", f2
!
!   if (trim(filetype)=="ascii") then
!      do iconf = startconf, nconf
!
!         write(f1, '(I5.5)') iconf
!         f2 = 'poslocm_'//trim(adjustl(f1))//'.ascii'
!         f3 = 'poslocm_'//trim(adjustl(f1))//'.dat'
!         f4 = 'poslocm_'//trim(adjustl(f1))//'_2.dat'
!         f5 = 'poslocm_'//trim(adjustl(f1))//'_result.txt'
!         f6 = 'poslocm_'//trim(adjustl(f1))//'_sphere.ascii'
!
!         inquire(file=trim(f2),exist=ex)
!         if (ex .eqv. .false.) then
!            write(*,*) 'file', f2 ,'does not exist'
!            CYCLE 
!         endif
!         open (10,file=trim(f2),status='old')
!            read(10,*) c1_energy(1,iconf), energy(iconf)
!         close(10)
!
!!
!! the unit cell is expanded in the subroutine findmolecules in order to find the atoms, which belong to the same molecule.
!! the coordinates of these atoms are safed in rxyz. the type of the atoms is safed in finalchar.
!!
!
!         call findmolecule(rxyz,alat,finalchar)
!
!!
!! Assign the Van-der-Waals radii
!!
!         do i = 1, nat*molecules
!            call sym2rvan(finalchar(i), rvan(i))
!         end do
!
!
!!
!!from this system a fingerprint is taken
!!
!         call periodic_fingerprint(rxyz,alat,finalchar,rvan,fpsall(1,1,iconf))
!
!
!
!         write(17,*) 'conf: ', iconf
!         do n=1,lseg*molecules_sphere*principleev
!            write(17,*) (fpsall(n,j,iconf),j=1,molecules)
!         enddo
!         write(17,*)
!
!      enddo
!   else
!      stop 'only ascii-files can be read'
!   endif
!   close(17)
!
!
!!   open (17,file='fingerprints.txt')
!!      do i=1, nconf
!!         write(17,*) 'conf: ', i
!!         do n=1,lseg*molecules_sphere*principleev
!!            write(17,*) (fpsall(n,j,i),j=1,molecules)
!!         enddo
!!         write(17,*)
!!      enddo
!!   close(17)
!
!
!
!
!
!!
!!short fingerprint
!!
!   open(32,file='results.txt')
!
!   write(*,*) 'short fingerprint'
!   do iconf=1, nconf
!      do jconf=iconf+1, nconf
!!         do imol=1,molecules
!!            do jmol=1,molecules
!!               tt=0.d0
!!               do i=1,lseg*molecules_sphere*principleev
!!                  tt=tt+(fpsall(i,imol,iconf)-fpsall(i,jmol,jconf))**2
!!               enddo
!!               cost(imol,jmol)=tt
!!            enddo
!!         enddo
!!         dtt=0.d0
!!
!!         !hungarian algorithm
!!         call apc(molecules, cost, iassign, dtt)
!         call get_distance_molgom(fpsall(:,:,iconf),fpsall(:,:,jconf),dtt,lseg,molecules,molecules_sphere,principleev)
!         endiff=abs(energy(iconf)-energy(jconf))
!        ! dtt=sqrt(dtt)
!         write(32,*) iconf,jconf,dtt,endiff
!      enddo
!   enddo
!   close(32)
!
!   if (clustering .eqv. .true.) then
!      write(*,*) 'calculate all fingerprint distances for clustering'
!      do iconf=1, nconf
!         do jconf=1, nconf
!            do imol=1,molecules
!               do jmol=1,molecules
!                  tt=0.d0
!                  do i=1,lseg*molecules_sphere*principleev
!                     tt=tt+(fpsall(i,imol,iconf)-fpsall(i,jmol,jconf))**2
!                  enddo
!                  cost(imol,jmol)=tt
!               enddo
!            enddo
!            dtt=0.d0
!
!            !hungarian algorithm
!            call apc(molecules,cost,iassign,dtt)
!            cluster_arr(jconf,iconf)=sqrt(dtt)
!         enddo
!      enddo
!
!
!      write(*,*) 'clustering started'
!      surviving=.true.
!      count=nconf
!      do while (count>cluster_number)
!         max=0.d0
!         do iconf=1, nconf
!            sum=0.d0
!            if (surviving(iconf) .eqv. .true.) then
!               do jconf=1, nconf
!                  if (iconf /= jconf .and. surviving(jconf) .eqv. .true.) then
!                     sum=sum+(1.d0/(cluster_arr(jconf,iconf)+1.d-100))
!                  endif
!               enddo
!               if (sum>max) then
!                  curpos = iconf
!                  max=sum
!               endif
!            endif
!         enddo
!         surviving(curpos)=.false.
!         write(*,*) curpos,'eliminated'
!         count=count-1
!      enddo
!
!      write(*,*)
!      write(*,*) 'surviving structures '
!      do iconf=1, nconf
!         if (surviving(iconf) .eqv. .true.) then
!            write(*,*) iconf
!         endif
!      enddo
!
!   endif
!
!
!end program molecular_fingerprint
!
subroutine get_distance_molgom(fp1,fp2,dist,lseg,molecules,molecules_sphere,principleev)
implicit none
integer:: imol,jmol,i
integer:: lseg,molecules,molecules_sphere,principleev
integer:: iassign(molecules)
real(8):: tt,dtt,dist
real(8):: cost(molecules,molecules)
real(8):: fp1(lseg*molecules_sphere*principleev,molecules),fp2(lseg*molecules_sphere*principleev,molecules)
         do imol=1,molecules
            do jmol=1,molecules
               tt=0.d0
               do i=1,lseg*molecules_sphere*principleev
                  tt=tt+(fp1(i,imol)-fp2(i,jmol))**2
               enddo
               cost(imol,jmol)=tt
            enddo
         enddo
         dtt=0.d0

         !hungarian algorithm
         call apc(molecules, cost, iassign, dtt)
         dist=dtt
end subroutine




subroutine create_contracted_om_1(width_overlap,principleev,nat,molecules,rxyz,rvan,amplitude,fp_t,lseg,write_files)
  implicit none

  logical :: write_files
  integer :: principleev, xyz, alpha, beta, mu, nu, k, i, j, m, n
  integer :: lwork, nat, molecules, lseg, info
  real*8 :: width_overlap
  real*8, dimension(3, molecules*nat):: rxyz
  real*8, dimension(molecules) :: amplitude
  real*8, dimension(3,nat) :: rxyz_temp
  real*8, dimension(molecules*nat) :: rvan
  real*8, dimension(nat) :: rvan_temp
  real*8, dimension(nat,principleev,molecules) :: em
  real*8, dimension(nat,nat) :: om
  real*8, dimension(nat*molecules,nat*molecules) :: om_b
  real*8, dimension(molecules,principleev,molecules,principleev) :: om_t
  real*8, dimension(nat) :: fp
  real*8, dimension(principleev*molecules) :: fp_t


  real*8, dimension(:), allocatable :: work


  if (lseg.eq.1) then
    em=0.d0

    lwork=max(1,3*nat-1)
    k=0
    do i=1, molecules
        do j=1,nat
            k=k+1
            do xyz=1,3
                rxyz_temp(xyz,j)=rxyz(xyz,k)
            enddo
            rvan_temp(j)=rvan(k)
        enddo
        om=0.d0
        call create_molom_1(nat,rxyz_temp,rvan_temp,om,width_overlap)
        allocate(work(lwork))
        fp=0.d0
        call DSYEV('V','L',nat,om,nat,fp,work,lwork,info)
        deallocate(work)

        if (write_files .eqv. .true.) then

            write(14,*) 'eigenvalues of molecule ', i
            do m=1, nat
                write(14,"(e9.4)")fp(m)
            enddo
            write(14,*)
        endif


        do m=1, principleev
            do n=1,nat
                em(n,m,i)=om(n,nat+1-m)
            enddo
        enddo
    enddo

   if(write_files .eqv. .true.) then
      do k=1, molecules
         write(14,*) "principle eigenvektors of molecule: ", k
         do i=1,nat
            write(14,"(6(1x,e9.2),I4)")(em(i,j,k),j=1,principleev), i
           enddo
         write(14,*)
      enddo
   endif

   om_b=0.d0
   call create_molom_1(nat*molecules,rxyz,rvan,om_b,width_overlap)


   om_t=0.d0
   do mu = 1, molecules
      do nu = 1, molecules
         do alpha = 1, principleev
            do beta = 1, principleev
               do i=1, nat
                  do j=1, nat
                     om_t(mu,alpha,nu,beta)=em(i,alpha,mu)*om_b(i+(mu-1)*nat,j+(nu-1)*nat)*em(j,beta,nu)+om_t(mu,alpha,nu,beta)
                  enddo
               enddo
               om_t(mu,alpha,nu,beta)=om_t(mu,alpha,nu,beta)*amplitude(mu)*amplitude(nu)
            enddo
         enddo
      enddo
   enddo



    fp_t=0.d0
    lwork=max(1,3*principleev*molecules-1)
    allocate(work(lwork))
    call DSYEV('N','L',principleev*molecules,om_t,principleev*molecules,fp_t,work,lwork,info)
    if (info/=0) stop 'eigenvalues of om_t not found'
    deallocate(work)




    if(write_files .eqv. .true.) then
      write(14,*)"eigenvalues of om_t"
      do i=1, molecules*principleev
         write(14,*) fp_t(i)
      enddo
      write(14,*)
    endif




  else
    write(*,*) "p- and s-orbitals used, but a method for s-orbitals called!"
    stop
  endif


end subroutine create_contracted_om_1
subroutine create_molom_1(nat,rxyz,rvan,om,width_overlap)
   implicit real*8 (a-h,o-z)

   dimension rxyz(3,nat),rvan(nat),om(nat,nat)
   real*8 :: width_overlap

   rvan=rvan*width_overlap
   ! Gaussian overlap
   !  <sj|si>
   do iat=1,nat
      xi=rxyz(1,iat)
      yi=rxyz(2,iat)
      zi=rxyz(3,iat)

      do jat=1,nat
         d2=(rxyz(1,jat) -xi)**2 +(rxyz(2,jat)-yi)**2+(rxyz(3,jat)-zi)**2
         r=.5d0/(rvan(iat)**2 + rvan(jat)**2)
         om(jat,iat)= sqrt(4.d0*r*(rvan(iat)*rvan(jat)))**3 * exp(-d2*r)
      enddo
   enddo
end subroutine create_molom_1


subroutine periodic_fingerprint(parini,rxyz,alat0,finalchar,rvan,fpsall,nat)
   use mod_parini, only: typ_parini
   use defs_basis, only: bohr_ang
   implicit none
   type(typ_parini), intent(in):: parini
   integer:: nat
   integer, dimension (parini%fp_18_expaparameter+1) :: shifting
   real*8, dimension (3,(parini%fp_18_expaparameter+1)**3) :: possibilites
   real*8, dimension(3,nat*parini%fp_18_molecules):: rxyz
   real*8, dimension(nat*parini%fp_18_molecules) :: rvan
   real*8, dimension(parini%fp_18_lseg*parini%fp_18_molecules_sphere*parini%fp_18_principleev,parini%fp_18_molecules) :: fpsall
   real*8, dimension(:,:,:), allocatable :: rxyz_b
   real*8, dimension(:,:), allocatable :: rxyz_b2
   real*8, dimension(:), allocatable :: rvan_b, amplitude, mol_amplitude, fp_t, fp_b
   character(len=2), dimension(:), allocatable :: finalchar_b
   real*8, dimension (3,3)::alat,alat0
   character(len=2), dimension(nat*parini%fp_18_molecules) :: finalchar
   logical, dimension(parini%fp_18_molecules,(parini%fp_18_expaparameter+1)**3) :: is_copied

   real*8 :: distance, convert, sum
   integer :: i, j, m, n, k, molecules_system, mol
   integer :: iat, curpos
   integer*4 :: position

   !parameters for cutoff function
   real*8 :: radius_cutoff, radius_cutoff2, factor_cutoff, factor2

write(*,*) nat,parini%fp_18_molecules
   !Convert to angstrom
   alat=alat0*bohr_ang





   !termination condition for to big parameters
   if ((parini%fp_18_expaparameter+1)**3>2147483647) then
      write(*,*) "expanding parameter is too large --> change integer type of position"
      stop
   endif

   !termination condition for odd epanding parameter
   if (mod(parini%fp_18_expaparameter,2)/=0) then
      write(*,*) "expanding parameter has to be even --> change parameter"
      stop
   endif

   allocate (rxyz_b(3, nat*parini%fp_18_molecules,(parini%fp_18_expaparameter+1)**3))

   !the unit cell is copied to rxyz_b(:,:,1)
   do i=1, nat*parini%fp_18_molecules
      do j=1,3
         rxyz_b(j,i,1)=rxyz(j,i)
      enddo
   enddo

   !first the sytstem size is expanded
   !the unit cell will be shifted into the the three directions from -fp_18_expaparameter to fp_18_expaparameter
   n=-(parini%fp_18_expaparameter/2)
   do i=1,parini%fp_18_expaparameter+1
      shifting(i)= n
      n=n+1
   enddo

   !write all permuations into possibilites(:,:)
   position=0
   do i=1, parini%fp_18_expaparameter+1
      do j=1, parini%fp_18_expaparameter+1
         do k=1, parini%fp_18_expaparameter+1
            position=position+1
            possibilites(1,position)=shifting(k)
            possibilites(2,position)=shifting(j)
            possibilites(3,position)=shifting(i)
         enddo
      enddo
   enddo

   !the expanded system is stored in rxyz_b(:,:,:)
   k=2
   do m = 1, (parini%fp_18_expaparameter+1)**3
      if (possibilites(1,m)/=0 .or. possibilites(2,m)/=0 .or. possibilites(3,m)/=0) then
         do j=1, nat*parini%fp_18_molecules
            do i=1,3
               if(i==1) then
                  rxyz_b(1,j,k)=possibilites(i,m)*alat(1,1)+rxyz_b(1,j,1)
               else if(i==2) then
                  rxyz_b(1,j,k)=possibilites(i,m)*alat(1,2)+rxyz_b(1,j,k)
                  rxyz_b(2,j,k)=possibilites(i,m)*alat(2,2)+rxyz_b(2,j,1)
               else
                  rxyz_b(1,j,k)=possibilites(i,m)*alat(1,3)+rxyz_b(1,j,k)
                  rxyz_b(2,j,k)=possibilites(i,m)*alat(2,3)+rxyz_b(2,j,k)
                  rxyz_b(3,j,k)=possibilites(i,m)*alat(3,3)+rxyz_b(3,j,1)
               endif
            enddo
         enddo
         k=k+1
      endif
   enddo


   !from here on a ascii file for the expanded system could be generated
!   if(write_files .eqv. .true.) then
!
!      f7 = 'poslocm_'//trim(adjustl(f1))//'_expand.ascii'
!      open (16,file=trim(f7))
!         write(16,*)
!         write(16,*) alat(1,1),alat(1,2),alat(2,2)
!         write(16,*) alat(1,3),alat(2,3),alat(3,3)
!
!         do i=1,(parini%fp_18_expaparameter+1)**3
!            do m=1,nat*fp_18_molecules
!               write(16,"(3E26.15E2,1A4)") (rxyz_b(j,m,i),j=1,3), finalchar(m)
!            enddo
!         enddo
!      close(16)
!
!      f8='poslocm_'//trim(adjustl(f1))//'_expand.dat'
!      open (17,file=trim(f8))
!         do m=1,fp_18_molecules
!            do n=1,nat
!               write(17,"(I2)")m
!            enddo
!         enddo
!
!         do i=2,(parini%fp_18_expaparameter+1)**3
!            do m=1,fp_18_molecules
!               do n=1,nat
!                  write(17,"(I2)")0
!               enddo
!            enddo
!         enddo
!      close(17)
!
!      f9='poslocm_'//trim(adjustl(f1))//'_expand_2.dat'
!      open (18,file=trim(f9))
!         do i=1,(parini%fp_18_expaparameter+1)**3
!            do m=1,fp_18_molecules
!               do n=1,nat
!                  write(18,"(I2)") m
!               enddo
!            enddo
!         enddo
!      close(18)
!   endif



!
!convert the system
!

   convert=1.d0/0.52917720859d0
   do m = 1, (parini%fp_18_expaparameter+1)**3
      do j=1, nat*parini%fp_18_molecules
         do i=1,3
            rxyz_b(i,j,m)=rxyz_b(i,j,m)*convert
         enddo
      enddo
   enddo

   do i=1,3
      do j=1,3
         alat(i,j) = alat(i,j)*convert
      enddo
   enddo


!cutoff radius is defined in order to reduce the system size with a cut off function
   radius_cutoff=sqrt(2.d0*parini%fp_18_nex_cutoff)*parini%fp_18_width_cutoff
   radius_cutoff2=radius_cutoff**2
   factor_cutoff=1.d0/(2.d0*parini%fp_18_nex_cutoff*parini%fp_18_width_cutoff**2)


!
!search for the atoms within radius_cutoff2 and create overlap matrix for the system
!!
!   if(write_files .eqv. .true.) then
!      f11 = 'poslocm_'//trim(adjustl(f1))//'_interaction.txt'
!      open (21,file=trim(f11))
!      open (14,file=trim(f5))
!   endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!cut a sphere of atoms from the geometrical center of the center molecuel
!then create a overlap matrix and applied a cutoff function
!

   do mol=1,parini%fp_18_molecules
      iat=(mol-1)*nat+1

      do i=1,(parini%fp_18_expaparameter+1)**3
         do j=1, parini%fp_18_molecules
            is_copied(j,i)=.false.
         enddo
      enddo

      molecules_system=0
      do i=1,(parini%fp_18_expaparameter+1)**3
         do j=1,parini%fp_18_molecules
            if (is_copied(j,i) .eqv. .false.) then
               do m=iat,iat+nat-1
                  if(is_copied(j,i) .eqv. .true.) exit
                  do n=1,nat
                     distance=0.d0
                     do k=1,3
                        distance=distance+(rxyz_b(k,m,1)-rxyz_b(k,n+(j-1)*nat,i))**2
                     enddo
                     if(distance<=(radius_cutoff2*(rvan(m)+rvan(n+(j-1)*nat))**2)) then
                        is_copied(j,i)=.true.
                        molecules_system=molecules_system+1
                        exit
                     endif
                  enddo
               enddo
            endif
         enddo
      enddo

      write (*,*) 'molecules in sphere:',molecules_system
      if (molecules_system>parini%fp_18_molecules_sphere) stop 'expand molecules_sphere'

!
!write an ascii files with the molecules within radius= radius_cutoff2*radius2
!

!      if(write_files .eqv. .true.) then
!         write(f10,'(I2)') mol
!         f10 = 'poslocm_'//trim(adjustl(f1))//trim(adjustl(f10))//'_spheremolecules.ascii'
!         open (20,file=trim(f10))
!            write(20,*)
!            write(20,*) alat(1,1)/convert,alat(1,2)/convert,alat(2,2)/convert
!            write(20,*) alat(1,3)/convert,alat(2,3)/convert,alat(3,3)/convert
!            do i=1,(parini%fp_18_expaparameter+1)**3
!               do n=1,fp_18_molecules
!                  if(is_copied(n,i) .eqv. .true.) then
!                     do m=1,nat
!                        write(20,"(3E26.15E2,1A4)") (rxyz_b(j,m+(n-1)*nat,i)/convert,j=1,3), finalchar(m+(n-1)*nat)
!                     enddo
!                  endif
!               enddo
!            enddo
!         close(20)
!      endif

      allocate (rxyz_b2(3, nat*molecules_system))
      allocate (rvan_b(nat*molecules_system))
      allocate (amplitude(nat*molecules_system))
      allocate (mol_amplitude(molecules_system))
      allocate (finalchar_b(nat*molecules_system))
      allocate (fp_t(molecules_system*parini%fp_18_principleev*parini%fp_18_lseg))
      allocate (fp_b(molecules_system*parini%fp_18_lseg*nat))


      !array initialization step
      amplitude=0.d0
      rxyz_b2=0.d0
      fp_t=0.d0
      fp_b=0.d0
      mol_amplitude=0.d0
      rvan_b=0.d0

      curpos=0
      do i=1,(parini%fp_18_expaparameter+1)**3
         do j=1,parini%fp_18_molecules
            if (is_copied(j,i) .eqv. .true.) then
               do m=1,nat
                  curpos=curpos+1

                  finalchar_b(curpos)=finalchar(m+(j-1)*nat)
                  rvan_b(curpos)=rvan(m+(j-1)*nat)

                  do k=1,3
                     rxyz_b2(k,curpos)=rxyz_b(k,m+(j-1)*nat,i)
                  enddo

                  do n=iat,iat+nat-1
                     distance=0.d0
                     do k=1,3
                        distance=distance+(rxyz_b(k,n,1)-rxyz_b(k,m+(j-1)*nat,i))**2
                     enddo
                     if (distance<=radius_cutoff2*(rvan(n)+rvan(m+(j-1)*nat))**2) then
                        factor2=1.d0/(rvan(n)+rvan(m+(j-1)*nat))**2
                        amplitude(curpos)=(1.d0-distance*factor_cutoff*factor2)**parini%fp_18_nex_cutoff+amplitude(curpos)
                     endif
                  enddo
               enddo
            endif
         enddo
      enddo


      curpos=0
      do i=1,molecules_system
         do j=1,nat
            curpos=curpos+1
            mol_amplitude(i)=mol_amplitude(i)+amplitude(curpos)
         enddo
      enddo

      sum=0.d0
      do i=1, molecules_system
         if (mol_amplitude(i)>sum) then
            sum=mol_amplitude(i)
         endif
      enddo

      do i=1,molecules_system
         mol_amplitude(i)=mol_amplitude(i)/sum
      enddo

      curpos=0
      do i=1,molecules_system
         do j=1,nat
            curpos=curpos+1
            amplitude(curpos)=mol_amplitude(i)
         enddo
      enddo
!
!
!      if(write_files .eqv. .true.) then
!         write(21,*) 'molecule: ', mol
!         do i=1, molecules_system
!            write(21,*) mol_amplitude(i)
!         enddo
!         write(21,*)
!
!         write(f12,'(I2)') mol
!         f12 = 'poslocm_'//trim(adjustl(f1))//trim(adjustl(f12))//'_spheremolecules.dat'
!         open (22,file=trim(f12))
!            do i=1,molecules_system*nat
!               write(22,*) amplitude(i)
!            enddo
!         close(22)
!      endif

      call create_contracted_om_1(parini%fp_18_width_overlap,parini%fp_18_principleev,nat,molecules_system,rxyz_b2,rvan_b,mol_amplitude,fp_t,&
      parini%fp_18_lseg,.false.)

      do i=1,molecules_system*parini%fp_18_lseg*parini%fp_18_principleev
         fpsall(i,mol)=fp_t(molecules_system*parini%fp_18_lseg*parini%fp_18_principleev+1-i)
      enddo

      do i=molecules_system*parini%fp_18_lseg*parini%fp_18_principleev,parini%fp_18_molecules_sphere*parini%fp_18_lseg*parini%fp_18_principleev
         fpsall(i,mol)=0.d0
      enddo

      deallocate(rxyz_b2,rvan_b,amplitude,mol_amplitude,finalchar_b,fp_t,fp_b)
   enddo


!   if(write_files .eqv. .true.) then
!      close(21)
!      close(14)
!   endif

   deallocate(rxyz_b)
end subroutine periodic_fingerprint
!Input is the structure, output whatever quantity was given originally
!Here the total number of atom is nat*molecules. We should call it nattot or something
subroutine findmolecule(parini,rxyz,alat0,finalchar,xred,char_type,typat,ntypat,nat)
    use mod_parini, only: typ_parini
   use defs_basis, only: bohr_ang
   implicit none
    type(typ_parini), intent(in):: parini
   integer, dimension(nat*parini%fp_18_molecules), intent(in) :: typat
   character(2), dimension(ntypat), intent(in) :: char_type
   integer, intent(in):: ntypat,nat
   real*8, dimension(3,3), intent(in) :: alat0
   real*8 :: alat(3,3)
   real*8, dimension(3,nat*parini%fp_18_molecules), intent(out) :: rxyz
   real*8, dimension(3,nat*parini%fp_18_molecules), intent(in)  :: xred
   character(len=2), dimension(nat*parini%fp_18_molecules), intent(out) :: finalchar

   integer :: molecule_nr,true_nr1,true_nr2,i,j,k,position,m,n,l
   real*8  :: radius
   real*8 :: distance
   real*8, dimension (3,2) :: copy
   real*8, dimension(3, nat*parini%fp_18_molecules,27) :: rxyz_b
   real*8, dimension (3,27) :: possibilites
   real*8, dimension(nat*parini%fp_18_molecules) :: rcov
   logical, dimension (2,nat*parini%fp_18_molecules,27) :: is_true
   logical, dimension (2,nat*parini%fp_18_molecules) :: first_search_is_true
   integer, dimension (nat*parini%fp_18_molecules) :: first_search_molecule_number
   integer, dimension (nat*parini%fp_18_molecules,27) :: molecule_number
   integer, dimension (nat*parini%fp_18_molecules,2) :: atom_number
   integer, dimension(parini%fp_18_molecules,2) :: finalmolecules
   integer, dimension (nat*parini%fp_18_molecules,27) :: finalmolecule_number
   character(len=2), dimension(nat*parini%fp_18_molecules,27) :: is_char
   logical :: file_exists

   integer, dimension(:,:), allocatable :: molecule_id

!Convert to angstrom
  alat=alat0*bohr_ang

!   write(*,*)'findmolecule: poslow'//trim(adjustl(f1))//'.ascii'
!
!   if (trim(filetype)=="ascii") then
!
!      !first the energies from the the ascii files are copied in energy(iconf)
!      open (11,file=trim(f2),status='old')
!         read(11,*)
!
!         !now the lattice cell is safed in alat
!         do i = 1, 2
!            read(11,*) (copy(j, i), j = 1, 3)
!         end do
!
!         alat(1,1)=copy(1,1)
!         alat(2, 1)=0
!         alat(3, 1)=0
!
!         do j=2,3
!            alat(j-1,2)=copy(j,1)
!         end do
!         alat(3,2)=0
!
!         do j=1,3
!            alat(j,3)=copy(j,2)
!         end do


!read the cartesian coordinates of the atoms to rxyz_b and their char value into is_char(i)
!         do i = 1, nat*molecules
!            read(11, *) (rxyz_b(j, i, 1), j = 1, 3), is_char(i,1)
!         end do
!Using MHM data
         write(*,*) nat,parini%fp_18_molecules
         call rxyz_int2cart(alat,xred,rxyz_b,nat*parini%fp_18_molecules)
         do i = 1, nat*parini%fp_18_molecules
            is_char(i,1)=char_type(typat(i))
         end do



         do i = 1, nat*parini%fp_18_molecules
            call sym2rcov(is_char(i,1), rcov(i))
         end do

!
!serach for atoms in the unit cell, that are connected and enumerate the compartments
!this is done in order to weight the molecules. copies of the same moluecules get eliminated according to their weight.
!

!initialise first_search_is_true with .false.
         do i=1,2
            do j=1,nat*parini%fp_18_molecules
               first_search_is_true(i,j)=.false.
            enddo
         enddo

         do i=1, nat*parini%fp_18_molecules
            first_search_molecule_number(i)=0
         enddo

         molecule_nr=0
         do i = 1, nat*parini%fp_18_molecules
            if(first_search_is_true(1,i) .eqv. .false. ) then
               molecule_nr=molecule_nr+1
               first_search_molecule_number(i)=molecule_nr
               first_search_is_true(1,i) = .true.
               first_search_is_true(2,i) = .true.

!calculate distance between the atom (i) and all other atoms

               do j = 1, nat*parini%fp_18_molecules
                  distance=0.d0
                  do l=1,3
                     distance= (rxyz_b(l,i,1)-rxyz_b(l,j,1))**2+distance
                  enddo
                  distance=sqrt(distance)
                  radius=(rcov(i)+rcov(j))*0.52917720859d0
                  radius=radius*1.1
                  if (distance<=radius) then
                     first_search_is_true(1,j)=.true.
                     first_search_molecule_number(j)=molecule_nr
                  endif
               enddo
            endif

            j=1
            do while (j<2)
               do k=1,nat*parini%fp_18_molecules
                  if (first_search_is_true(1,k) .eqv. .true. .and. first_search_is_true(2,k) .eqv. .false.) then
                     first_search_is_true(2,k)=.true.
                     do m=1,nat*parini%fp_18_molecules
                        distance=0.d0
                        do n=1,3
                           distance= (rxyz_b(n,k,1)-rxyz_b(n,m,1))**2+distance
                        enddo
                        distance=sqrt(distance)
                        radius=(rcov(k)+rcov(m))*0.52917720859d0
                        radius=radius*1.1
                        if (distance<=radius) then
                           first_search_is_true(1,m)=.true.
                           first_search_molecule_number(m)=molecule_nr
                        endif
                     enddo
                  endif
               enddo
! check if number of trues in is_true(1,:,:) is equal to number of trues in is_true(1,:,:), else set j = 1
               true_nr1=0
               true_nr2=0
               do k = 1, nat*parini%fp_18_molecules
                  if(first_search_is_true(1,k) .eqv. .true.) then
                     true_nr1=true_nr1+1
                     if(first_search_is_true(2,k) .eqv. .true.) then
                        true_nr2=true_nr2+1
                     endif
                  endif
               enddo
               if(true_nr1/=true_nr2) then
                  j=0
               endif
               j=j+1
            enddo
         enddo

!initialise the rest of the rxyz_b array with zeros
         do j=2,27
            do i=1,nat*parini%fp_18_molecules
               do n=1,3
                  rxyz_b(n,i,j)=0.d0
               enddo
            enddo
         enddo
         
!permutation of the different possibilities are calculated
         position=0
         do i=1,3
            do j=1,3
               do k=1,3
                  position=position+1
                  select case(k)
                     case(1)
                        possibilites(1,position)= -1
                     case(2)
                        possibilites(1,position)=  0
                     case(3)
                        possibilites(1,position)=  1
                  end select

                  select case(j)
                     case(1)
                        possibilites(2,position)= -1
                     case(2)
                        possibilites(2,position)=  0
                     case(3)
                        possibilites(2,position)=  1
                     end select

                  select case(i)
                     case(1)
                        possibilites(3,position)= -1
                     case(2)
                        possibilites(3,position)=  0
                     case(3)
                        possibilites(3,position)=  1
                  end select
               enddo
            enddo
         enddo


         k=2
         do m = 1, 27
            if (possibilites(1,m)/=0 .or. possibilites(2,m)/=0 .or. possibilites(3,m)/=0) then
               do j=1, nat*parini%fp_18_molecules
                  do i=1,3
                     if(i==1) then
                        rxyz_b(1,j,k)=possibilites(i,m)*alat(1,1)+rxyz_b(1,j,1)
                     else if(i==2) then
                        rxyz_b(1,j,k)=possibilites(i,m)*alat(1,2)+rxyz_b(1,j,k)
                        rxyz_b(2,j,k)=possibilites(i,m)*alat(2,2)+rxyz_b(2,j,1)
                     else
                        rxyz_b(1,j,k)=possibilites(i,m)*alat(1,3)+rxyz_b(1,j,k)
                        rxyz_b(2,j,k)=possibilites(i,m)*alat(2,3)+rxyz_b(2,j,k)
                        rxyz_b(3,j,k)=possibilites(i,m)*alat(3,3)+rxyz_b(3,j,1)
                     endif
                  enddo
               enddo
               k=k+1
            endif
         enddo

         do position=2,27
            do n=1, nat*parini%fp_18_molecules
!
!            if (write_files .eqv. .true.) then
!               write(11,"(3E26.15E2,1A4)") (rxyz_b(j,n,position),j=1,3), is_char(n,1)
!            endif

            is_char(n,position)=is_char(n,1)

            enddo
         enddo

      close(11)


!search for connected atoms within radius (2*largest cov_radius of sytem)+10%
!initialise is_true with .false.
      do i=1,2
         do j=1,nat*parini%fp_18_molecules
            do k=1,27
               is_true(i,j,k)=.false.
            enddo
         enddo
      enddo



      do i=1,27
         do j=1,nat*parini%fp_18_molecules
            molecule_number(j,i)=0
         enddo
      enddo


      molecule_nr=0
      do i = 1, nat *parini%fp_18_molecules
         if(is_true(1,i,1) .eqv. .false. ) then
            molecule_nr=molecule_nr+1
            molecule_number(i,1)=molecule_nr
            is_true(1,i,1) = .true.
            is_true(2,i,1) = .true.

!calculate distance between the atom (i,1) and all other atoms
            do k=1,27
               do j = 1, nat*parini%fp_18_molecules
                  distance=0.d0
                  do l=1,3
                     distance= (rxyz_b(l,i,1)-rxyz_b(l,j,k))**2+distance
                  enddo
                  distance=sqrt(distance)
                  radius=(rcov(i)+rcov(j))*0.52917720859d0
                  radius=radius*1.1
                  if (distance<=radius) then
                     is_true(1,j,k)=.true.
                     molecule_number(j,k)=molecule_nr
                  endif
               enddo
            enddo
         endif
         j=1
         do while (j<28)
            do k=1,nat*parini%fp_18_molecules
               if (is_true(1,k,j) .eqv. .true. .and. is_true(2,k,j) .eqv. .false.) then
                  is_true(2,k,j)=.true.
                  do l=1,27
                     do m=1,nat*parini%fp_18_molecules
                        distance=0.d0
                        do n=1,3
                           distance= (rxyz_b(n,k,j)-rxyz_b(n,m,l))**2+distance
                        enddo
                        distance=sqrt(distance)
                        radius=(rcov(k)+rcov(m))*0.52917720859d0
                        radius=radius*1.1
                        if (distance<=radius) then
                           is_true(1,m,l)=.true.
                           molecule_number(m,l)=molecule_nr
                        endif
                     enddo
                  enddo
               endif
            enddo
! check if number of trues in is_true(1,:,:) is equal to number of trues in is_true(1,:,:), else set j = 1
            if (j==27) then
               true_nr1=0
               true_nr2=0
               do k = 1, nat*parini%fp_18_molecules
                  do l =1, 27
                     if(is_true(1,k,l) .eqv. .true.) then
                        true_nr1=true_nr1+1
                        if(is_true(2,k,l) .eqv. .true.) then
                           true_nr2=true_nr2+1
                        endif
                     endif
                  enddo
               enddo
               if(true_nr1/=true_nr2) then
                  j=0
               endif
            endif
            j=j+1
         enddo
      enddo



!!the molecule number of each atom is written into a .dat file in order to colorize the v_sim simulation
!      if(write_files .eqv. .true.) then
!         open (12,file=trim(f3))
!            do j=1, 27
!               do i=1, nat*parini%fp_18_molecules
!                  write(12,"(I2)") (molecule_number(i,j))
!               enddo
!            enddo
!         close(12)
!      endif

!initialise atom_number with 0
      do i=1,nat*parini%fp_18_molecules
         do j=1,2
            atom_number(i,j)=0
         enddo
      enddo

!count the atoms of each molecule and save this number in the atom_number(i,:). save additional the number of the molecule into (:,j)
      do j=1, nat*parini%fp_18_molecules
         if(molecule_number(j,1)>=0) then
            k=molecule_number(j,1)
            atom_number(k,2)=k
            atom_number(k,1)=atom_number(k,1)+1
         endif
      enddo

!all molecules, that don't have the number of "nat" atoms are deleted
      do j=1,nat*parini%fp_18_molecules
         if (atom_number(j,2)/=0) then
            m=0
            do n=1,27
               do i=1,nat*parini%fp_18_molecules
                  if (molecule_number(i,n)==atom_number(j,2)) then
                     m=m+1
                  endif
               enddo
            enddo
            if (m/=nat) then
               atom_number(j,1)=0
               atom_number(j,2)=0
            endif
         endif
      enddo


!bubblesort the number of atoms to know of which molecules the number of atoms is the largest in the unit cell
      do j = nat*parini%fp_18_molecules-1, 1, -1
         do i = 1, j
            if (atom_number(i,1) < atom_number(i+1,1)) then
               k = atom_number(i,1)
               m = atom_number(i,2)
               atom_number(i,1) = atom_number(i+1,1)
               atom_number(i+1,1) = k
               atom_number(i,2) = atom_number(i+1,2)
               atom_number(i+1,2) = m
            endif
         enddo
      enddo

!count the number of molecules in atom_number
      molecule_nr=0
      do i = 1, nat*parini%fp_18_molecules
         if (atom_number(i,2)/=0) then
            molecule_nr=molecule_nr+1
         endif
      enddo
!if molecule_nr ==0 then stop the program because no atoms are found in the file
      if (molecule_nr == 0) then
         write(*,*)"no atoms found in the .ascii file"
         stop
      endif

!each molecule gets its own id_key which is composed of the first_search_molecule_number
!this key gets ordered to filter out the "same" molecules
      allocate (molecule_id(nat,molecule_nr))
      do n=1,molecule_nr
         m=0
         do i=1,27
            do j=1,nat*parini%fp_18_molecules
               if (molecule_number(j,i)==atom_number(n,2)) then
                  m=m+1
                  if(m>nat) then
                     write(*,*) "searching aglortihm found too much atoms for 1 molecule "
                     write(*,*) "--> change searching radius in findmolecule subroutine"
                     stop
                  endif
                     molecule_id(m,n)=first_search_molecule_number(j)
               endif
            enddo
         enddo
      enddo

      do l=1,molecule_nr
         do j = nat-1, 1, -1
            do i = 1, j
               if (molecule_id(i,l) < molecule_id(i+1,l)) then
                  n = molecule_id(i,l)
                  molecule_id(i,l) = molecule_id(i+1,l)
                  molecule_id(i+1,l) = n
               endif
            enddo
         enddo
      enddo

      do i=1,2
         do j=1,parini%fp_18_molecules
            finalmolecules(j,i)=0
         enddo
      enddo

      do i=1,molecule_nr
         do j=1, parini%fp_18_molecules
            if (finalmolecules(j,1)==molecule_id(1,i)) then
               exit
            else if (finalmolecules(j,1)==0) then
               finalmolecules(j,1)=molecule_id(1,i)
               finalmolecules(j,2)=i
               exit
            endif
         enddo
      enddo
      deallocate(molecule_id)

      do i=1,27
         do j=1,nat*parini%fp_18_molecules
            finalmolecule_number(j,i)=0
         enddo
      enddo

!
!check if there are enough molecules in the system
!
      m=0
      do i=1,parini%fp_18_molecules
         if (finalmolecules(i,2)/=0) then
            m=m+1
         endif
      enddo
      if (m/=parini%fp_18_molecules) then
!         write(*,*) 'there are too less molecules in poslow'//trim(adjustl(f1))//'.ascii'
         stop
      endif


      m=1
      do i=1,parini%fp_18_molecules
         do j=1,27
            do k=1,nat*parini%fp_18_molecules
               if(atom_number(finalmolecules(i,2),2)==molecule_number(k,j)) then
                  if (m>nat*parini%fp_18_molecules+1) then
                     write(*,*) "searching aglortihm found too much atoms for molecules in the system "
                     write(*,*) "--> change searching radius in findmolecule subroutine"
                     stop
                  endif

                  do l=1,3
                     rxyz(l,m)=rxyz_b(l,k,j)
                  enddo
                  finalchar(m)=is_char(k,j)
                  finalmolecule_number(k,j)=(molecule_number(k,j))
                  m=m+1
               endif
            enddo
         enddo
      enddo

      if (m/=nat*parini%fp_18_molecules+1) then
         write(*,*) "searching aglortihm found too less atoms for existing molecules"
         write(*,*) "--> change searching radius in findmolecule subroutine"
         stop
      endif
!
!      if(write_files .eqv. .true.) then
!         open (13,file=trim(f4))
!            do i=1,27
!               do j=1,nat*fp_18_molecules
!                  write(13,"(I2)") (finalmolecule_number(j,i))
!               enddo
!            enddo
!         close(13)
!      endif

!   endif

end subroutine findmolecule
!SUBROUTINE APC(N,A,F,Z)
!implicit none
!! Modified by Ali Sadeghi to get real*8 matrix A(N,N) and converted to F90 
!!
!! SOLUTION OF THE LINEAR MIN-SUM ASSIGNMENT PROBLEM.
!! HUNGARIAN METHOD. COMPLEXITY O(N**3).
!!
!! MEANING OF THE INPUT PARAMETERS:
!! N      = NUMBER OF ROWS AND COLUMNS OF THE COST MATRIX.
!! A(I,J) = COST OF THE ASSIGNMENT OF ROW  I  TO COLUMN  J .
!! ON RETURN, THE INPUT PARAMETERS ARE UNCHANGED.
!!
!! MEANING OF THE OUTPUT PARAMETERS:
!! F(I) = COLUMN ASSIGNED TO ROW  I .
!! Z    = COST OF THE OPTIMAL ASSIGNMENT =
!!      = A(1,F(1)) + A(2,F(2)) + ... + A(N,F(N)) .
!!
!!
!! THE CODE IS BASED ON THE HUNGARIAN METHOD AS DESCRIBED BY
!! LAWLER (COMBINATORIAL OPTIMIZATION : NETWORKS AND
!! MATROIDS, HOLT, RINEHART AND WINSTON, NEW YORK, 1976).
!! THE ALGORITHMIC PASCAL-LIKE DESCRIPTION OF THE CODE IS
!! GIVEN IN G.CARPANETO, S.MARTELLO AND P.TOTH, ALGORITHMS AND
!! CODES FOR THE ASSIGNMENT PROBLEM, ANNALS OF OPERATIONS
!! RESEARCH 7, 1988.
!!
!! SUBROUTINE APC DETERMINES THE INITIAL DUAL AND PARTIAL
!! PRIMAL SOLUTIONS AND THEN SEARCHES FOR AUGMENTING PATHS
!! UNTIL ALL ROWS AND COLUMNS ARE ASSIGNED.
!!
!! MEANING OF THE MAIN INTERNAL VARIABLES:
!! FB(J) = ROW ASSIGNED TO COLUMN  J .
!! M     = NUMBER OF INITIAL ASSIGNMENTS.
!! U(I)  = DUAL VARIABLE ASSOCIATED WITH ROW  I .
!! V(J)  = DUAL VARIABLE ASSOCIATED WITH COLUMN  J .
!!
!! APC NEEDS THE FOLLOWING SUBROUTINES: INCR
!!                                      INIT
!!                                      PATH
!!
!! THIS WORK WAS SUPPORTED BY  C.N.R. , ITALY.
!      INTEGER n
!      REAL(8)  A(n,n),Z,U(n),V(n)
!      integer F(N),FB(n), RC(n)
!      INTEGER M,I,J
!! SEARCH FOR THE INITIAL DUAL AND PARTIAL PRIMAL SOLUTIONS.
!      CALL INIT(N,A,F,M,U,V,FB,RC)
!      IF ( M .NE. N ) then 
!! SOLUTION OF THE REDUCED PROBLEM.
!      DO  I=1,N
!        IF ( F(I) == 0 ) THEN 
!! DETERMINATION OF AN AUGMENTING PATH STARTING FROM ROW  I .
!        CALL PATH(N,A,I,F,J,U,V,FB,RC)
!! ASSIGNMENT OF ROW  I  AND COLUMN  J .
!        CALL INCR(n,F,J,FB,RC)
!        ENDIF
!      ENDDO    
!      ENDIF
!
!! COMPUTATION OF THE SOLUTION COST  Z .
!      Z = sum(u(1:N)) + sum(V(1:N))
!      END
!
!
!      SUBROUTINE INCR(n,F,J,FB,RC)
!!
!! ASSIGNMENT OF COLUMN  J .
!!
!      INTEGER n,I,J,JJ,  F(n),FB(n),RC(n)
!   10 I = RC(J)
!      FB(J) = I
!      JJ = F(I)
!      F(I) = J
!      J = JJ
!      IF ( J > 0 ) GO TO 10
!      RETURN
!      END
!
!
!      SUBROUTINE INIT(N,A,F,M,U,V,FB,P)
!!
!! SEARCH FOR THE INITIAL DUAL AND PARTIAL PRIMAL SOLUTIONS.
!!
!! P(I) = FIRST UNSCANNED COLUMN OF ROW  I .
!       IMPLICIT NONE
!!
!      INTEGER n,m, F(n),FB(n),P(n)
!      real(8) A(n,n) , U(n),V(n)
!      REAL(8), parameter :: INF = 1.d9
!      real(8) min, IA
!      integer i,j, k,R, JMIN, KK
!! PHASE 1 .
!      M = 0
!      F(1:N)=0
!      FB(1:N)=0
!! SCANNING OF THE COLUMNS ( INITIALIZATION OF  V(J) ).
!      DO 40 J=1,N
!        MIN = INF
!        DO 30 I=1,N
!          IA = A(I,J)
!          IF ( IA .GT. MIN ) GO TO 30
!          IF ( IA .LT. MIN ) GO TO 20
!          IF ( F(I) .NE. 0 ) GO TO 30
!   20     MIN = IA
!          R = I
!   30   CONTINUE
!        V(J) = MIN
!        IF ( F(R) .NE. 0 ) GO TO 40
!! ASSIGNMENT OF COLUMN  J  TO ROW  R .
!        M = M + 1
!        FB(J) = R
!        F(R) = J
!        U(R) = 0.d0
!        P(R) = J + 1
!   40 CONTINUE
!! PHASE 2 .
!! SCANNING OF THE UNASSIGNED ROWS ( UPDATING OF  U(I) ).
!      DO 110 I=1,N
!        IF ( F(I) .NE. 0 ) GO TO 110
!        MIN = INF
!        DO 60 K=1,N
!          IA = A(I,K) - V(K)
!          IF ( IA .GT. MIN )  GO TO 60
!          IF ( IA .LT. MIN )  GO TO 50
!          IF ( FB(K) .NE. 0 ) GO TO 60
!          IF ( FB(J) .EQ. 0 ) GO TO 60
!   50     MIN = IA
!          J = K
!   60   CONTINUE
!        U(I) = MIN
!        JMIN = J
!        IF ( FB(J) .EQ. 0 ) GO TO 100
!        DO 80 J=JMIN,N
!          IF ( A(I,J) - V(J) .GT. MIN ) GO TO 80
!          R = FB(J)
!          KK = P(R)
!          IF ( KK .GT. N ) GO TO 80
!          DO 70 K=KK,N
!            IF ( FB(K) .GT. 0 ) GO TO 70
!            IF ( A(R,K) - U(R) - V(K) .EQ. 0.d0 ) GO TO 90
!   70     CONTINUE
!          P(R) = N + 1
!   80   CONTINUE
!        GO TO 110
!! REASSIGNMENT OF ROW  R  AND COLUMN  K .
!   90   F(R) = K
!        FB(K) = R
!        P(R) = K + 1
!! ASSIGNMENT OF COLUMN  J  TO ROW  I .
!  100   M = M + 1
!        F(I) = J
!        FB(J)= I
!        P(I) = J + 1
!  110 CONTINUE
!      RETURN
!      END
!      SUBROUTINE PATH(N,A,II,F,JJ,U,V,FB,RC)
!!
!! DETERMINATION OF AN AUGMENTING PATH STARTING FROM
!! UNASSIGNED ROW  II  AND TERMINATING AT UNASSIGNED COLUMN
!! JJ , WITH UPDATING OF DUAL VARIABLES  U(I)  AND  V(J) .
!!
!! MEANING OF THE MAIN INTERNAL VARIABLES:
!! LR(L) = L-TH LABELLED ROW ( L=1,NLR ).
!! PI(J) = MIN ( A(I,J) - U(I) - V(J) , SUCH THAT ROW  I  IS
!!         LABELLED AND NOT EQUAL TO  FB(J) ).
!! RC(J) = ROW PRECEDING COLUMN  J  IN THE CURRENT
!!         ALTERNATING PATH.
!! UC(L) = L-TH UNLABELLED COLUMN ( L=1,NUC ).
!!
!      implicit none
!      INTEGER N 
!      real(8)  A(n,n),U(n),V(N),PI(n), IA, MIN
!      INTEGER F(N),LR(n),UC(n)
!      INTEGER FB(n),RC(n)
!      REAL(8), parameter :: INF = 1.d9
!      integer  i,j,k,L, ii,jj, NUC , NLR, R
!! INITIALIZATION.
!      LR(1) = II
!      DO 10 K=1,N
!        PI(K) = A(II,K) - U(II) - V(K)
!        RC(K) = II
!        UC(K) = K
!   10 CONTINUE
!      NUC = N
!      NLR = 1
!      GO TO 40
!! SCANNING OF THE LABELLED ROWS.
!   20 R = LR(NLR)
!      DO 30 L=1,NUC
!        J = UC(L)
!        IA = A(R,J) - U(R) - V(J)
!        IF ( IA .GE. PI(J) ) GO TO 30
!        PI(J) = IA
!        RC(J) = R
!   30 CONTINUE
!! SEARCH FOR A ZERO ELEMENT IN AN UNLABELLED COLUMN.
!   40 DO 50 L=1,NUC
!        J = UC(L)
!        IF ( PI(J) .EQ. 0.d0 ) GO TO 100
!   50 CONTINUE
!! UPDATING OF THE DUAL VARIABLES  U(I)  AND  V(J) .
!      MIN = INF
!      DO 60 L=1,NUC
!        J = UC(L)
!        IF ( MIN .GT. PI(J) ) MIN = PI(J)
!   60 CONTINUE
!      DO 70 L=1,NLR
!        R = LR(L)
!        U(R) = U(R) + MIN
!   70 CONTINUE
!      DO 90 J=1,N
!        IF ( PI(J) .EQ. 0.d0 ) GO TO 80
!        PI(J) = PI(J) - MIN
!        GO TO 90
!   80   V(J) = V(J) - MIN
!   90 CONTINUE
!      GO TO 40
!  100 IF ( FB(J) .EQ. 0 ) GO TO 110
!! LABELLING OF ROW  FB(J)  AND REMOVAL OF THE LABEL  OF COLUMN  J .
!      NLR = NLR + 1
!      LR(NLR) = FB(J)
!      UC(L) = UC(NUC)
!      NUC = NUC - 1
!      GO TO 20
!! DETERMINATION OF THE UNASSIGNED COLUMN  J .
!  110 JJ = J
!      RETURN
!      END
 

! returns the covalent radius of atom with chemical symbol sym 
subroutine sym2rcov(sym,rcov)
  real(8)  :: rcov
  character(len=2) :: sym  ! chemical symbol 
  select case (trim(adjustl(sym)))
     ! covalet radius in Angstrom taken from WebElements: http://www.webelements.com/periodicity/covalent_radius/                      
  case('H')
     rcov= 0.37d0 
  case('He')
     rcov= 0.32d0 
  case('Li')
     rcov= 1.34d0 
  case('Be')
     rcov= 0.90d0 
  case('B')
     rcov= 0.82d0 
  case('C')
     rcov= 0.77d0 
  case('N')
     rcov= 0.75d0 
  case('O')
     rcov= 0.73d0 
  case('F')
     rcov= 0.71d0 
  case('Ne')
     rcov= 0.69d0 
  case('Na')
     rcov= 1.54d0 
  case('Mg')
     rcov= 1.30d0 
  case('Al')
     rcov= 1.18d0 
  case('Si')
     rcov= 1.11d0 
  case('P')
     rcov= 1.06d0 
  case('S')
     rcov= 1.02d0 
  case('Cl')
     rcov= 0.99d0 
  case('Ar')
     rcov= 0.97d0 
  case('K')
     rcov= 1.96d0 
  case('Ca')
     rcov= 1.74d0 
  case('Sc')
     rcov= 1.44d0 
  case('Ti')
     rcov= 1.36d0 
  case('V')
     rcov= 1.25d0 
  case('Cr')
     rcov= 1.27d0 
  case('Mn')
     rcov= 1.39d0 
  case('Fe')
     rcov= 1.25d0 
  case('Co')
     rcov= 1.26d0 
  case('Ni')
     rcov= 1.21d0 
  case('Cu')
     rcov= 1.38d0 
  case('Zn')
     rcov= 1.31d0 
  case('Ga')
     rcov= 1.26d0 
  case('Ge')
     rcov= 1.22d0 
  case('As')
     rcov= 1.19d0 
  case('Se')
     rcov= 1.16d0 
  case('Br')
     rcov= 1.14d0 
  case('Kr')
     rcov= 1.10d0 
  case('Rb')
     rcov= 2.11d0 
  case('Sr')
     rcov= 1.92d0 
  case('Y')
     rcov= 1.62d0 
  case('Zr')
     rcov= 1.48d0 
  case('Nb')
     rcov= 1.37d0 
  case('Mo')
     rcov= 1.45d0 
  case('Tc')
     rcov= 1.56d0 
  case('Ru')
     rcov= 1.26d0 
  case('Rh')
     rcov= 1.35d0 
  case('Pd')
     rcov= 1.31d0 
  case('Ag')
     rcov= 1.53d0 
  case('Cd')
     rcov= 1.48d0 
  case('In')
     rcov= 1.44d0 
  case('Sn')
     rcov= 1.41d0 
  case('Sb')
     rcov= 1.38d0 
  case('Te')
     rcov= 1.35d0 
  case('I')
     rcov= 1.33d0 
  case('Xe')
     rcov= 1.30d0 
  case('Cs')
     rcov= 2.25d0 
  case('Ba')
     rcov= 1.98d0 
  case('La')
     rcov= 1.69d0 
     !     case('Ce')
     !     case('Pr')
     !     case('Nd')
     !     case('Pm')
     !     case('Sm')
     !     case('Eu')
     !     case('Gd')
     !     case('Tb')
     !     case('Dy')
     !     case('Ho')
     !     case('Er')
     !     case('Tm')
     !     case('Yb')
  case('Lu')
     rcov= 1.60d0 
  case('Hf')
     rcov= 1.50d0 
  case('Ta')
     rcov= 1.38d0 
  case('W')
     rcov= 1.46d0 
  case('Re')
     rcov= 1.59d0 
  case('Os')
     rcov= 1.28d0 
  case('Ir')
     rcov= 1.37d0 
  case('Pt')
     rcov= 1.28d0 
  case('Au')
     rcov= 1.44d0 
  case('Hg')
     rcov= 1.49d0 
  case('Tl')
     rcov= 1.48d0 
  case('Pb')
     rcov= 1.47d0 
  case('Bi')
     rcov= 1.46d0 
     !     case('Po')
     !     case('At')
  case('Rn')
     rcov= 1.45d0 
  case('LA') 
     rcov= 1.122462048309373d0
  case('LB')
     rcov= 0.9877666025122482d0
     !     case('Fr')
     !     case('Ra')
     !     case('Ac')
     !     case('Th')
     !     case('Pa')
     !     case('U')
     !     case('Np')
     !     case('Pu')
     !     case('Am')
     !     case('Cm')
  case default                     
     print*, " Not recognized atomic type "//sym ; stop
  endselect

  rcov = rcov /  0.52917720859d0   ! convert to atomic units 

endsubroutine sym2rcov


! returns the Van-der-Waals radius of atom with chemical symbol sym
subroutine sym2rvan(sym,rvan)
  real(8)  :: rvan
  character(len=2) :: sym  ! chemical symbol 
  select case (trim(adjustl(sym)))
     ! Radius in Angstrom taken from WebElements: http://www.webelements.com/periodicity/covalent_radius/
  case('H')
     rvan= 1.20d0
  case('He')
     rvan= 1.40d0
  case('Li')
     rvan= 1.82d0
  case('Be')
     rvan= 0.00d0
  case('B')
     rvan= 0.00d0
  case('C')
     rvan= 1.70d0
  case('N')
     rvan= 1.55d0
  case('O')
     rvan= 1.52d0
  case('F')
     rvan= 1.47d0
  case('Ne')
     rvan= 1.54d0
  case('Na')
     rvan= 2.27d0
  case('Mg')
     rvan= 1.73d0
  case('Al')
     rvan= 0.00d0
  case('Si')
     rvan= 2.10d0
  case('P')
     rvan= 1.80d0
  case('S')
     rvan= 1.80d0
  case('Cl')
     rvan= 1.75d0
  case('Ar')
     rvan= 1.88d0
  case('K')
     rvan= 2.75d0
  case('Ca')
     rvan= 0.00d0
  case('Sc')
     rvan= 0.00d0
  case('Ti')
     rvan= 0.00d0
  case('V')
     rvan= 0.00d0
  case('Cr')
     rvan= 0.00d0
  case('Mn')
     rvan= 0.00d0
  case('Fe')
     rvan= 0.00d0
  case('Co')
     rvan= 0.00d0
  case('Ni')
     rvan= 1.63d0
  case('Cu')
     rvan= 1.40d0
  case('Zn')
     rvan= 1.39d0
  case('Ga')
     rvan= 1.87d0
  case('Ge')
     rvan= 0.00d0
  case('As')
     rvan= 1.85d0
  case('Se')
     rvan= 1.90d0
  case('Br')
     rvan= 1.85d0
  case('Kr')
     rvan= 2.02d0
  case('Rb')
     rvan= 0.00d0
  case('Sr')
     rvan= 0.00d0
  case('Y')
     rvan= 0.00d0
  case('Zr')
     rvan= 0.00d0
  case('Nb')
     rvan= 0.00d0
  case('Mo')
     rvan= 0.00d0
  case('Tc')
     rvan= 0.00d0
  case('Ru')
     rvan= 0.00d0
  case('Rh')
     rvan= 0.00d0
  case('Pd')
     rvan= 1.63d0
  case('Ag')
     rvan= 1.72d0
  case('Cd')
     rvan= 1.58d0
  case('In')
     rvan= 1.93d0
  case('Sn')
     rvan= 2.17d0
  case('Sb')
     rvan= 0.00d0
  case('Te')
     rvan= 2.06d0
  case('I')
     rvan= 1.98d0
  case('Xe')
     rvan= 2.16d0
  case('Cs')
     rvan= 0.00d0
  case('Ba')
     rvan= 0.00d0
  case('La')
     rvan= 0.00d0
     !     case('Ce')
     !     case('Pr')
     !     case('Nd')
     !     case('Pm')
     !     case('Sm')
     !     case('Eu')
     !     case('Gd')
     !     case('Tb')
     !     case('Dy')
     !     case('Ho')
     !     case('Er')
     !     case('Tm')
     !     case('Yb')
  case('Lu')
     rvan= 0.00d0
  case('Hf')
     rvan= 0.00d0
  case('Ta')
     rvan= 0.00d0
  case('W')
     rvan= 0.00d0
  case('Re')
     rvan= 0.00d0
  case('Os')
     rvan= 0.00d0
  case('Ir')
     rvan= 0.00d0
  case('Pt')
     rvan= 0.00d0
  case('Au')
     rvan= 0.00d0
  case('Hg')
     rvan= 0.00d0
  case('Tl')
     rvan= 0.00d0
  case('Pb')
     rvan= 0.00d0
  case('Bi')
     rvan= 0.00d0
     !     case('Po')
     !     case('At')
  case('Rn')
     rvan= 0.00d0
  case('LA') 
     rvan= 0.00d0
  case('LB')
     rvan= 0.00d0
     !     case('Fr')
     !     case('Ra')
     !     case('Ac')
     !     case('Th')
     !     case('Pa')
     !     case('U')
     !     case('Np')
     !     case('Pu')
     !     case('Am')
     !     case('Cm')
  case default                     
     print*, " Not recognized atomic type "//sym ; stop
  endselect

  if (rvan.EQ.0) then
     print*, " Van-der-Waals radius is zero! " ; stop
  end if

  rvan = rvan /  0.52917720859d0   ! convert to atomic units

endsubroutine sym2rvan

