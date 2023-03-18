program main
implicit none
integer:: i,n,k,nvert,nstrmax,istr,ifile,nfile
integer,allocatable::iwk(:),vertex(:),nat(:),nattypat(:,:),ntypat(:),nstr(:),typat(:),frac(:,:)
real(8),allocatable:: enthalpies(:,:),concentrations(:),confrac(:),znucl(:,:),x(:),y(:)
character(250):: tmp,folder,filename
character(450):: all_line
real(8):: znucl_types(2),e_elemental(2)
real(8):: tmp_r,minconc,maxconc,amu,rcov
integer:: tmp_i
logical:: found,nat_found,ntypat_found,file_exists
character(15):: fdy,fdx,tmp_c
real(8),parameter:: Ha_eV=27.21138386d0
character(2):: sym_elemental(2),fnn1,fnn2

nstrmax=1000
  file_exists=.false.
  filename="folders"
  INQUIRE(FILE=trim(filename), EXIST=file_exists)
     if(.not.file_exists) then
             stop "Provide the folders file"
     endif


!Count the number of files 
open(unit=4,file="folders")
nfile=0
do 
read(4,'(a)',end=99,err=99) folder
nfile=nfile+1
enddo
99 continue
close(4)
if(nfile.lt.2) stop "At least two runs must be provided"

!Allocate the files
allocate(enthalpies(nstrmax,nfile),concentrations(nfile),confrac(nfile),znucl(2,nfile),&
         &nat(nfile),ntypat(nfile),nattypat(2,nfile),nstr(nfile),frac(2,nfile))
!Initiallize
enthalpies=1.d10
concentrations=0.d0
confrac=0.d0
znucl=0.d0
nat=0
ntypat=0
nattypat=0

open(unit=4,file="folders")
do ifile=1,nfile
   ntypat_found=.false.
   read(4,'(a)',end=11,err=11) folder
!Get the energy of a compound
!Use the Enthalpy file if refined
   file_exists=.false.
   filename=trim(folder)//"/Enthalpies_poslow"
   INQUIRE(FILE=trim(filename), EXIST=file_exists)
   if(file_exists) then
      open(unit=66,file=trim(folder)//"/Enthalpies_poslow")
      nstr(ifile)=0
      read(66,*,end=13,err=13) tmp_c
      do while(.true.)
         nstr(ifile)=nstr(ifile)+1
         read(66,*,end=13,err=13) tmp_c,tmp_r,enthalpies(nstr(ifile),ifile)
         write(*,*) trim(folder), enthalpies(nstr(ifile),ifile)
      enddo
      13 continue
      close(66)
      enthalpies(:,ifile)=enthalpies(:,ifile)/Ha_eV
      nstr(ifile)=nstr(ifile)-1
   else 
      file_exists=.false.
      filename=trim(folder)//"/earr.dat"
      INQUIRE(FILE=trim(filename), EXIST=file_exists)
      if(.not.file_exists) stop "Enthalpies_poslow or earr.dat must be provided!"
      open(unit=5,file=trim(folder)//"/earr.dat")
      read(5,*) nstr(ifile)
      read(5,*) tmp
      nstr(ifile)=min(nstr(ifile),nstrmax)
      do istr=1,nstr(ifile)
        read(5,*)n,enthalpies(istr,ifile)
      enddo
      close(5)
   endif
!Get the composition
   open(unit=12,file=trim(folder)//"/params_new.in")
     do while(.true.)
     read(12,'(a450)',end=22)all_line
     n = len_trim(all_line)
  !NAT
     call parsescalar_int("NAT",3,all_line(1:n),n,nat(ifile),found)
     if(found) nat_found=.true.
     if(found) cycle
  !NTYPE
     call parsescalar_int("NTYPE",5,all_line(1:n),n,ntypat(ifile),found)
     if(found) ntypat_found=.true.
     if(found) cycle
   enddo
   22 continue
   close(12)
if(.not.nat_found.or..not.ntypat_found) then
   write(*,*) "NAT and NTYPE must be provided in params_new.in"
   stop
endif
!Allocate the arrays
   allocate(typat(nat(ifile)))
   open(unit=12,file=trim(folder)//"/params_new.in")
   found=.false.
   do while(.true.)
      read(12,'(a450)',end=33)all_line
      n = len_trim(all_line)
   !Znucl
      call parsearray_real("ZNUCL",5,all_line(1:n),n,znucl(1:ntypat(ifile),ifile),ntypat(ifile),found)
      if(found) cycle
   !TYPAT
      call parsearray_int("TYPAT",5,all_line(1:n),n,typat(1:nat(ifile)),nat(ifile),found)
      if(found) cycle
   enddo
   33 continue
   close(12)
   !Count the types
   do i=1,nat(ifile)
      nattypat(typat(i),ifile)=nattypat(typat(i),ifile)+1
   enddo
   deallocate(typat)
enddo

11 continue

!Get the two znucl that are around
   znucl_types(1) = MINVAL(znucl, MASK = znucl .GT.0.d0) 
   znucl_types(2) = MAXVAL(znucl, MASK = znucl .GT.0.d0) 

!Sort the arrays according to the znucl
do ifile=1,nfile
   if(znucl(1,ifile).ne.znucl_types(1)) then
           tmp_r=znucl(1,ifile)
           znucl(1,ifile)=znucl(2,ifile)
           znucl(2,ifile)=tmp_r
           tmp_i=nattypat(1,ifile)
           nattypat(1,ifile)=nattypat(2,ifile)
           nattypat(2,ifile)=tmp_i
    endif
enddo

!Get elemental symbol
call atmdata(amu,rcov,sym_elemental(1),maxval(znucl(1,:)))
call atmdata(amu,rcov,sym_elemental(2),maxval(znucl(2,:)))


!Get the concentration in terms of znucl_types(1)
do ifile=1,nfile
   concentrations(ifile)=real(nattypat(1,ifile),8)/real(nat(ifile),8)
   confrac(ifile)=real(nattypat(1,ifile),8)/real(nattypat(2,ifile),8)
   call findfrac(confrac(ifile),frac(1,ifile),frac(2,ifile)) 
   enthalpies(:,ifile)=enthalpies(:,ifile)/real(nat(ifile),8)
enddo
minconc=minval(concentrations)
maxconc=maxval(concentrations)
!Get the elemental (or nearest to elemental) energies

e_elemental=1.d10
do ifile=1,nfile
   if(concentrations(ifile)==minconc) then
           e_elemental(2)=min(minval(enthalpies(:,ifile)),e_elemental(2))
   endif
   if(concentrations(ifile)==maxconc) then
           e_elemental(1)=min(minval(enthalpies(:,ifile)),e_elemental(1))
   endif
enddo


!Here we generate the scaled enthalpies
do ifile=1,nfile
   enthalpies(:,ifile)=enthalpies(:,ifile)-concentrations(ifile)*e_elemental(1)-(1.d0-concentrations(ifile))*e_elemental(2)
   enthalpies(:,ifile)=enthalpies(:,ifile)*Ha_eV
enddo




!To generate the plot
n=nfile
allocate(x(n),y(n),iwk(n),vertex(n))
do ifile=1,nfile
   x(ifile)=concentrations(ifile)
   y(ifile)=min(minval(enthalpies(:,ifile)),0.d0)   
enddo


open(unit=2,file="hull_list")
do ifile=1,nfile
   do i=1,nstr(ifile)
      write(2,*) x(ifile),enthalpies(i,ifile)
   enddo
enddo
call envelope(x, y, n, vertex, nvert, iwk)

open(unit=3,file="hull")
do i=1,nvert
k=modulo(i,nvert)+1
write(3,*) x(vertex(k)),y(vertex(k))
enddo

!Write a plotting script
open(unit=88,file="PlotHull.gnu")
write(88,*) "set style line  1 lt 1 lw 1.5 pt 6 ps 2 lc rgb '#000000'"
write(88,*) "# start value for H"
write(88,*) "h1 = 117/360.0"
write(88,*) "# end value for H"
write(88,*) "h2 = 227/360.0"
write(88,*) "# creating the palette by specifying H,S,V"
write(88,*) "set palette model HSV functions (1-gray)*(h2-h1)+h1,1,0.8"
write(88,*) "unset colorbox"
write(88,*) "set xrange[-0.05:1.05]"
write(fdx,'(es15.7)') FLOOR(10.d0*(minval(enthalpies)-0.3d0))/10.d0
write(fdy,'(es15.7)') CEILING(10.d0*(maxval(enthalpies,MASK = enthalpies .LT.1.d5)))/10.d0
write(88,*) "set yrange["//fdx//":"//fdy//"]"
do ifile=1,nfile
  write(fdx,'(es15.7)') x(ifile)
  write(fdy,'(es15.7)') y(ifile)-0.2d0
  if(concentrations(ifile)==0.d0) then
     write(fnn1,'(a)')  ""
     write(fnn2,'(a)')  ""
  elseif(concentrations(ifile)==1.d0) then
     write(fnn1,'(a)')  ""
     write(fnn2,'(a)')  ""
  else
     write(fnn1,'(i2)')  frac(1,ifile)
     write(fnn2,'(i2)')  frac(2,ifile)
  endif
  if(.not.(ifile.gt.1.and.any(concentrations(1:ifile-1)==concentrations(ifile)))) then
    if(concentrations(ifile)==0.d0) then
  write(88,'(a)') "set label '"//adjustl(trim(sym_elemental(2)))//"' at first "//fdx//","//fdy//" right rotate by 90"
    elseif(concentrations(ifile)==1.d0) then
  write(88,'(a)') "set label '"//adjustl(trim(sym_elemental(1)))//"' at first "//fdx//","//fdy//" right rotate by 90"
    else
  write(88,'(a)') "set label '"//adjustl(trim(sym_elemental(2)))//adjustl(trim(fnn2))//&
  &adjustl(trim(sym_elemental(1)))//adjustl(trim(fnn1))//"' at first "//fdx//","//fdy//" right rotate by 90"
    endif
  endif
enddo
write(88,*) "set ylabel 'Energy (eV)'"
write(88,*) "set xlabel 'Fraction "//sym_elemental(1)//"'"
write(88,*) "plot 'hull_list'  u 1:2:1  lc palette notitle"
write(88,*) "replot 'hull' u 1:2 w lp ls 1 notitle"

end program

subroutine findfrac(dec,num,denom)
implicit none
real(8):: dec,frac,error
integer:: inum,idenom,num,denom

error=1.d10
do inum=0,255
  do idenom=1,255
    frac=real(inum,8)/real(idenom,8)
    if(abs(frac-dec).lt.error) then
      error=abs(frac-dec)
      num=inum
      denom=idenom
    endif
  enddo
enddo
end

