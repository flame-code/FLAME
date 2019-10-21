program driver_espresso
implicit none
character(250):: all_line,arg
integer:: niggli
character(40):: filename,filename_out_struct,filename_out_kpt
niggli=0
if(iargc().lt.1) stop "No arguments passed"
   CALL getarg(1, arg)
   WRITE (*,*) "Input file ",trim(arg)
   filename=trim(arg)
if(iargc().ge.2) then
   CALL getarg(2, arg)
   if(trim(arg)=="niggli") then
      niggli=1
   else
      stop "The second argument must be empty or 'niggli' for the Niggli reduction"
   endif
endif
filename_out_struct="espresso.STRUCT"
filename_out_kpt="espresso.KPOINTS"
call rewrite_struct_espresso(filename,filename_out_struct,filename_out_kpt,niggli)
end program



  subroutine rewrite_struct_espresso(filename_in,filename_out_struct,filename_out_kpt,niggli)
  !use defs_basis
  !Since its a single call, we only have forces and stresses from one configuration!
  implicit none
  integer:: io,i,iat,n,k,l,m,i_tmp,nat,ntypat,itype,niggli
  real(8):: energy,strten(6),value,latvec(3,3),str_matrix(3,3),vol,a(3,3),scaling,alat,r_tmp
  real(8),allocatable:: fcart(:,:),xred(:,:),xcart(:,:),amu(:)
  integer,allocatable:: fix(:,:)
  character(2),allocatable:: char_type(:),typat_char(:)
  character(11):: ch_tmp
  character(150)::all_line
  character(40) ::filename_in,filename_out_struct,filename_out_kpt
  real(8):: dkpt1,dkpt2,dkpt_12(2)
  integer:: ka,kb,kc,kpt_abc(3)
  logical:: tmpl,auto_kpt,found
  !if vasp is version 5.x, use vasp_5=.true.
  
  
  ch_tmp="old"
  open(unit=32,file=trim(filename_in))
  !First loop to get the number of atoms, number of types, etc  
  do while(.true.)
   read(32,'(a150)',end=29)all_line
  !!write(*,*) all_line
   n = len_trim(all_line)
   k = index(all_line(1:n),"number of atoms/cell      =")
    if(k.ne.0) then
      !write(*,*) "NAT"
      l = scan(all_line(1:n),"=",.true.)
      read(all_line(l+1:n),*) nat
      cycle
    endif
    k = index(all_line(1:n),"number of atomic types    =")
    if(k.ne.0) then
      !write(*,*) "NTYPAT"
      l = scan(all_line(1:n),"=",.true.)
      read(all_line(l+1:n),*) ntypat
      cycle
    endif
   k = index(all_line(1:n),"nat =")
    if(k.ne.0) then
      !write(*,*) "NAT"
      l = scan(all_line(1:n),"=",.true.)
      read(all_line(l+1:n),*) nat
      cycle
    endif
   k = index(all_line(1:n),"ntyp =")
    if(k.ne.0) then
      !write(*,*) "NTYP"
      l = scan(all_line(1:n),"=",.true.)
      read(all_line(l+1:n),*) ntypat
      cycle
    endif


   enddo
29 continue
   close(32)
!Allocate arrays
   allocate(fix(3,nat),fcart(3,nat),xred(3,nat),xcart(3,nat),char_type(ntypat),amu(ntypat),typat_char(nat))

  energy=1.d10
  fcart=1.d10
  strten=1.d10

!Do second loop to get the other stuff
  open(unit=32,file=trim(filename_in))
  !First loop to get the number of atoms, number of types, etc  
  do while(.true.)
   read(32,'(a150)',end=99)all_line
  !!write(*,*) all_line
   n = len_trim(all_line)
   k = index(all_line(1:n),"total   stress ")
    if(k.ne.0) then
      !write(*,*) "STRESSES FOUND"
      read(32,*)str_matrix(:,1)
      read(32,*)str_matrix(:,2)
      read(32,*)str_matrix(:,3)
      !write(*,*)str_matrix(:,1)
      !write(*,*)str_matrix(:,2)
      !write(*,*)str_matrix(:,3)
       strten(1)=-str_matrix(1,1)
       strten(2)=-str_matrix(2,2)
       strten(3)=-str_matrix(3,3)
       strten(6)=-str_matrix(1,2)
       strten(5)=-str_matrix(1,3)
       strten(4)=-str_matrix(2,3)
      cycle
    endif
   
   k = index(all_line(1:n),"!    total energy              =")
    if(k.ne.0) then
            !write(*,*) "ENERGY FOUND"
            read(all_line(1:n),*) ch_tmp,ch_tmp,ch_tmp,ch_tmp,energy
            !write(*,*) energy
    cycle
    endif

   k = index(all_line(1:n),"lattice parameter (alat)")
    if(k.ne.0) then
    !write(*,*) "alat found"
    !write(*,*) all_line(1:n)
      l = scan(all_line(1:n),"=",.true.)
      read(all_line(l+1:n),*) alat
    !write(*,*) alat
    cycle
    endif
   
   k = index(all_line(1:n),"crystal axes: (cart. coord. in units of alat)")
    if(k.ne.0) then
      read(32,'(a150)',end=99)all_line
      m = len_trim(all_line)
      l = scan(all_line(1:m),"(",.true.)
      read(all_line(l+1:m),*) latvec(:,1)
      read(32,'(a150)',end=99)all_line
      read(all_line(l+1:m),*) latvec(:,2)
      read(32,'(a150)',end=99)all_line
      read(all_line(l+1:m),*) latvec(:,3)
    !write(*,*) latvec
    latvec=latvec*alat
    cycle
    endif

   k = index(all_line(1:n),"atomic species   valence    mass     pseudopotential")
    if(k.ne.0) then
            do itype=1,ntypat
            read(32,*,end=99) char_type(itype),r_tmp,amu(itype)
            enddo 
    cycle
    endif



   k = index(all_line(1:n),"site n.     atom                  positions (alat units)")
    if(k.ne.0) then
      !write(*,*) "xcart found"
      do iat=1,nat
      read(32,'(a150)',end=99)all_line
      m = len_trim(all_line)
      l = scan(all_line(1:m),"(",.true.)
      read(all_line(1:m),*)i_tmp,typat_char(iat)
      read(all_line(l+1:m),*) xcart(:,iat)
      enddo
      !write(*,*) xcart
      xcart=xcart*alat
      call rxyz_cart2int(latvec,xred,xcart,nat)
    cycle
    endif
   k = index(all_line(1:n),"ATOMIC_SPECIES")
    if(k.ne.0) then
            do itype=1,ntypat
            read(32,*,end=99) char_type(itype),amu(itype)
            enddo 
    !write(*,*) char_type
    cycle
    endif
    k = index(all_line(1:n),"CELL_PARAMETERS")
    if(k.ne.0) then
      read(32,*,end=99) latvec(:,1)
      read(32,*,end=99) latvec(:,2)
      read(32,*,end=99) latvec(:,3)
    !write(*,*) latvec
    cycle
    endif

   k = index(all_line(1:n),"ATOMIC_POSITIONS")
    if(k.ne.0) then
      fix=1
      !write(*,*) "xcart found"
      do iat=1,nat
      read(32,'(a150)',end=99)all_line
      read(all_line,*,iostat=io) typat_char(iat), xred(:,iat),fix(:,iat)
      if(io.lt.0)  fix(:,iat)=1
!        do itype=1,ntypat
!          if(trim(ch_tmp)==trim(char_type(itype))) typat=itype
!        enddo
      enddo
      !write(*,*) xred
      !write(*,*) fix
    cycle
    endif
  enddo
  
  99 continue 
  close(32)
  
  !Transform all to bohr
  energy=energy*0.5d0
  strten=strten*0.5d0
  fcart=fcart*0.5d0


!Get infor from params_new.in
open(unit=12,file="params_new.in")
  do while(.true.)
   read(12,'(a150)',end=97)all_line
   n = len_trim(all_line)
!Block KPT****************
!AUTO_KPT
   call parse_logical("AUTO_KPT",8,all_line(1:n),n,auto_kpt,found)
   if(found) cycle
!KPTMESH
   call parsearray_int("KPTMESH",7,all_line(1:n),n,kpt_abc(1:3),3,found)
   if(found) then
     ka=kpt_abc(1)
     kb=kpt_abc(2)
     kc=kpt_abc(3)
   endif
   if(found) cycle
!KPTDEN
   call parsearray_real("KPTDEN",6,all_line(1:n),n,dkpt_12(1:2),2,found)
   if(found) then
     dkpt1=dkpt_12(1)
     dkpt2=dkpt_12(2)
   endif
   if(found) cycle
!Block KPT****************
enddo
 97 continue
  if(AUTO_KPT) then
    ka=0;kb=0;kc=0
  else
    dkpt1=0.d0
    dkpt2=0.d0
  endif

 close(12)
!Do niggli reduction
if(niggli==1) write(*,*) "Performing Niggli reduction..."
if(niggli==1) call fixcell_niggli(nat,latvec,xred)

!Set up kpt
if(dkpt1.ne.0.d0) then
call find_kpt(ka,kb,kc,latvec,dkpt1)
endif

!Create Block with atoms, structure, etc.
    open(unit=87,file=trim(filename_out_struct))
        write(87,'(a)') "ATOMIC_SPECIES"
        do itype=1,ntypat
           write(87,'(a,2x,f10.5,2x,a)') trim(char_type(itype)),amu(itype),trim(char_type(itype))//".PSP"
        enddo
        write(87,'(a)') "ATOMIC_POSITIONS crystal"
        do iat=1,nat
              write(87,'(a,2x,3(es25.15),2x,3(i5))') trim(typat_char(iat)),xred(:,iat),fix(:,iat)
        enddo
        write(87,'(a)') "CELL_PARAMETERS bohr"
           write(87,'(3(es25.15))') latvec(:,1)
           write(87,'(3(es25.15))') latvec(:,2)
           write(87,'(3(es25.15))') latvec(:,3)
    close(87)
    open(unit=87,file=trim(filename_out_kpt))
        write(87,'(a)') "K_POINTS automatic"
           write(87,'(i5,i5,i5,5x,i5,i5,i5)') ka,kb,kc,0,0,0
    close(87)
  end subroutine
subroutine fixcell_niggli(nat,latvec,xred)
implicit none
!This routine will apply the niggli reduction to the cell and to the reduced coordinates and transform them back into the cell
integer:: nat,iat
real(8):: latvec(3,3),xred(3,nat),latvec_out(3,3),eps,transmat(3,3),invmat(3,3),epsvol,a(3,3),vol

epsvol=1.d-6
    a=latvec
    vol= a(1,1)*a(2,2)*a(3,3)-a(1,1)*a(2,3)*a(3,2)-a(1,2)*a(2,1)*a(3,3)+&
         a(1,2)*a(2,3)*a(3,1)+a(1,3)*a(2,1)*a(3,2)-a(1,3)*a(2,2)*a(3,1)
    eps=epsvol*vol**(1.d0/3.d0)
    call niggli(latvec,latvec_out,transmat,eps)
    call invertmat(transmat,invmat,3)
    do iat=1,nat
       xred(:,iat)=matmul(invmat,xred(:,iat))
       xred(1,iat)=modulo(xred(1,iat),1.d0)
       xred(2,iat)=modulo(xred(2,iat),1.d0)
       xred(3,iat)=modulo(xred(3,iat),1.d0)
    enddo
latvec=latvec_out
end subroutine


subroutine niggli(latvec_in,latvec_out,transmat,eps)
!This soubroutine will produce the niggli reduced cell based on the improved algorithm of R. W. Grosse-Kunstleve,* N. K. Sauter and P. D. Adams
!On exit, latvec_out will contain the reduced cell, and transmat will provide the transformation operator in latvec_in
implicit none
real(8):: latvec_in(3,3),latvec_out(3,3),transmat(3,3),eps,tmpmat(3,3)
real(8):: pi,nigmat(6),dist_ang(6),nigmat_check(6)
integer:: iout,l,m,n
logical:: def_gt_0,feq,flt,fgt,is_niggli_cell,debug

!debug=.true.
debug=.false.

pi=acos(-1.d0)
transmat=0.d0
transmat(1,1)=1.d0
transmat(2,2)=1.d0
transmat(3,3)=1.d0
!0 step
!Initiallize the niggli matrix
call init_nigmat(latvec_in,nigmat,l,m,n,eps,pi)
write(888,'(a,6(1x,es15.7))') "step 0",nigmat
!1 step
1000 continue
    ! A1
    if (fgt(nigmat(1),nigmat(2),eps) .or. (feq(nigmat(1), nigmat(2), eps) .and. fgt(abs(nigmat(4)), abs(nigmat(5)),eps))) then
      call a1_action(nigmat,tmpmat,eps)
      transmat=matmul(transmat,tmpmat)
if(debug) then
      write(888,'(a,6(1x,es15.7))') "step 1",nigmat
!Check the transformation matrix
            latvec_out=matmul(latvec_in,transmat)
            call init_nigmat(latvec_out,nigmat_check,l,m,n,eps,pi)
            if(.not.all(abs(nigmat_check-nigmat).lt.1.d-7)) then
            write(*,*) "Transformation matrix wrong"
            write(888,'(a,6(1x,es15.7))') "step 1",nigmat_check
            stop
            endif
endif
    endif
    ! A2
    if (fgt(nigmat(2), nigmat(3),eps) .or. (feq(nigmat(2),nigmat(3),eps) .and. fgt(abs(nigmat(5)), abs(nigmat(6)),eps))) then
      call a2_action(nigmat,tmpmat,eps)
      transmat=matmul(transmat,tmpmat)
if(debug) then
      write(888,'(a,6(1x,es15.7))') "step 2",nigmat
!Check the transformation matrix
            latvec_out=matmul(latvec_in,transmat)
            call init_nigmat(latvec_out,nigmat_check,l,m,n,eps,pi)
            if(.not.all(abs(nigmat_check-nigmat).lt.1.d-7)) then
            write(*,*) "Transformation matrix wrong"
            write(888,'(a,6(1x,es15.7))') "step 2",nigmat_check
            stop
            endif
endif
      goto 1000
    endif
    ! A3
    if (def_gt_0(nigmat,eps)) then
      call a3_action(nigmat,tmpmat,eps)
      transmat=matmul(transmat,tmpmat)
if(debug) then
      write(888,'(a,6(1x,es15.7))') "step 3",nigmat
!Check the transformation matrix
            latvec_out=matmul(latvec_in,transmat)
            call init_nigmat(latvec_out,nigmat_check,l,m,n,eps,pi)
            if(.not.all(abs(nigmat_check-nigmat).lt.1.d-7)) then
            write(*,*) "Transformation matrix wrong"
            write(888,'(a,6(1x,es15.7))') "step 3",nigmat_check
            stop
            endif
endif
    else
    ! A4
      call a4_action(nigmat,tmpmat,eps)
      transmat=matmul(transmat,tmpmat)
if(debug) then
      write(888,'(a,6(1x,es15.7))') "step 4",nigmat
!Check the transformation matrix
            latvec_out=matmul(latvec_in,transmat)
            call init_nigmat(latvec_out,nigmat_check,l,m,n,eps,pi)
            if(.not.all(abs(nigmat_check-nigmat).lt.1.d-7)) then
            write(*,*) "Transformation matrix wrong"
            write(888,'(a,6(1x,es15.7))') "step 4",nigmat_check
            stop
            endif
endif
    endif 
    ! A5
    if (fgt(abs(nigmat(4)), nigmat(2),eps) &
        &.or. (feq(nigmat(4), nigmat(2),eps) .and. flt(nigmat(5)+nigmat(5), nigmat(6),eps))&
        &.or. (feq(nigmat(4),-nigmat(2),eps) .and. flt(nigmat(6), 0.d0,eps))) then
      call a5_action(nigmat,tmpmat,eps)
      transmat=matmul(transmat,tmpmat)
if(debug) then
      write(888,'(a,6(1x,es15.7))') "step 5",nigmat
!Check the transformation matrix
            latvec_out=matmul(latvec_in,transmat)
            call init_nigmat(latvec_out,nigmat_check,l,m,n,eps,pi)
            if(.not.all(abs(nigmat_check-nigmat).lt.1.d-7)) then
            write(*,*) "Transformation matrix wrong"
            write(888,'(a,6(1x,es15.7))') "step 5",nigmat_check
            stop
            endif
endif
      goto 1000
    endif
    ! A6
    if (fgt(abs(nigmat(5)), nigmat(1),eps)&
        &.or. (feq(nigmat(5), nigmat(1),eps) .and. flt(nigmat(4)+nigmat(4), nigmat(6),eps))&
        &.or. (feq(nigmat(5),-nigmat(1),eps) .and. flt(nigmat(6), 0.d0,eps))) then
      call a6_action(nigmat,tmpmat,eps)
      transmat=matmul(transmat,tmpmat)
if(debug) then
      write(888,'(a,6(1x,es15.7))') "step 6",nigmat
!Check the transformation matrix
            latvec_out=matmul(latvec_in,transmat)
            call init_nigmat(latvec_out,nigmat_check,l,m,n,eps,pi)
            if(.not.all(abs(nigmat_check-nigmat).lt.1.d-7)) then
            write(*,*) "Transformation matrix wrong"
            write(888,'(a,6(1x,es15.7))') "step 6",nigmat_check
            stop
            endif
endif
      goto 1000
    endif
    ! A7
    if (fgt(abs(nigmat(6)),nigmat(1),eps)&
        &.or. (feq(nigmat(6), nigmat(1),eps) .and. flt(nigmat(4)+nigmat(4),nigmat(5),eps))&
        &.or. (feq(nigmat(6),-nigmat(1),eps) .and. flt(nigmat(5), 0.d0,eps))) then
      call a7_action(nigmat,tmpmat,eps)
      transmat=matmul(transmat,tmpmat)
if(debug) then
      write(888,'(a,6(1x,es15.7))') "step 7",nigmat
!Check the transformation matrix
            latvec_out=matmul(latvec_in,transmat)
            call init_nigmat(latvec_out,nigmat_check,l,m,n,eps,pi)
            if(.not.all(abs(nigmat_check-nigmat).lt.1.d-7)) then
            write(*,*) "Transformation matrix wrong"
            write(888,'(a,6(1x,es15.7))') "step 7",nigmat_check
            stop
endif
            endif
      goto 1000
    endif
    ! A8
    if (flt(nigmat(4)+nigmat(5)+nigmat(6)+nigmat(1)+nigmat(2),0.d0,eps)&
        &.or. (feq(nigmat(4)+nigmat(5)+nigmat(6)+nigmat(1)+nigmat(2),0.d0,eps)) .and.&
        & fgt(nigmat(1)+nigmat(1)+nigmat(5)+nigmat(5)+nigmat(6),0.d0,eps))then
      call a8_action(nigmat,tmpmat,eps)
      transmat=matmul(transmat,tmpmat)
if(debug) then
      write(888,'(a,6(1x,es15.7))') "step 8",nigmat
!Check the transformation matrix
            latvec_out=matmul(latvec_in,transmat)
            call init_nigmat(latvec_out,nigmat_check,l,m,n,eps,pi)
            if(.not.all(abs(nigmat_check-nigmat).lt.1.d-7)) then
            write(*,*) "Transformation matrix wrong"
            write(888,'(a,6(1x,es15.7))') "step 8",nigmat_check
            stop
            endif
endif
      goto 1000
    endif
latvec_out=matmul(latvec_in,transmat)
if(debug) then
!Check the transformation matrix
    if(.not.is_niggli_cell(nigmat,eps)) stop "We did not generate a niggli cell"
    call init_nigmat(latvec_out,nigmat_check,l,m,n,eps,pi)
    if(.not.is_niggli_cell(nigmat_check,eps))then
    write(*,*) "No niggli cell in check!"
    endif
    if(.not.all(abs(nigmat_check-nigmat).lt.1.d-7)) then
    write(*,*) "Transformation matrix wrong"
    nigmat=nigmat_check
    stop
    endif
endif
end subroutine

subroutine  a1_action(nigmat,tmpmat,eps)
implicit none
real(8):: nigmat(6),tmpmat(3,3),eps
call n1_action(nigmat,tmpmat,eps)
end subroutine
subroutine  a2_action(nigmat,tmpmat,eps)
implicit none
real(8):: nigmat(6),tmpmat(3,3),eps
call n2_action(nigmat,tmpmat,eps)
end subroutine
subroutine  a3_action(nigmat,tmpmat,eps)
implicit none
real(8):: nigmat(6),tmpmat(3,3),eps
call n3_true_action(nigmat,tmpmat,eps)
end subroutine
subroutine  a4_action(nigmat,tmpmat,eps)
implicit none
real(8):: nigmat(6),tmpmat(3,3),eps
call n3_false_action(nigmat,tmpmat,eps)
end subroutine
subroutine  a5_action(nigmat,tmpmat,eps)
implicit none
real(8):: nigmat(6),tmpmat(3,3),eps
    if (nigmat(4) > 0.d0) then
      tmpmat(1,:)=(/1.d0,0.d0, 0.d0/)
      tmpmat(2,:)=(/0.d0,1.d0,-1.d0/)
      tmpmat(3,:)=(/0.d0,0.d0, 1.d0/)
      nigmat(3)=nigmat(3)+ nigmat(2)- nigmat(4)
      nigmat(4)=nigmat(4)-(nigmat(2)+nigmat(2))
      nigmat(5)=nigmat(5)-(nigmat(6))
    else
      tmpmat(1,:)=(/1.d0,0.d0, 0.d0/)
      tmpmat(2,:)=(/0.d0,1.d0, 1.d0/)
      tmpmat(3,:)=(/0.d0,0.d0, 1.d0/)
      nigmat(3)=nigmat(3)+ nigmat(2)+ nigmat(4)
      nigmat(4)=nigmat(4)+(nigmat(2)+nigmat(2))
      nigmat(5)=nigmat(5)+(nigmat(6))
    endif
    if(nigmat(3).lt.0.d0) stop "Error in action a5"
end subroutine

subroutine  a6_action(nigmat,tmpmat,eps)
implicit none
real(8):: nigmat(6),tmpmat(3,3),eps
logical:: flt,fgt,feq
    if (nigmat(5) > 0.d0) then
      tmpmat(1,:)=(/1.d0,0.d0,-1.d0/)
      tmpmat(2,:)=(/0.d0,1.d0, 0.d0/)
      tmpmat(3,:)=(/0.d0,0.d0, 1.d0/)
      nigmat(3)=nigmat(3)+ nigmat(1)- nigmat(5)
      nigmat(4)=nigmat(4)-(nigmat(6))
      nigmat(5)=nigmat(5)-(nigmat(1)+nigmat(1))
    else
      tmpmat(1,:)=(/1.d0,0.d0, 1.d0/)
      tmpmat(2,:)=(/0.d0,1.d0, 0.d0/)
      tmpmat(3,:)=(/0.d0,0.d0, 1.d0/)
      nigmat(3)=nigmat(3)+ nigmat(1)+ nigmat(5)
      nigmat(4)=nigmat(4)+(nigmat(6))
      nigmat(5)=nigmat(5)+(nigmat(1)+nigmat(1))
    endif
    if(nigmat(3).lt.0.d0) stop "Error in action a6" 
end subroutine

subroutine  a7_action(nigmat,tmpmat,eps)
implicit none
real(8):: nigmat(6),tmpmat(3,3),eps
    if (nigmat(6) > 0.d0) then
      tmpmat(1,:)=(/1.d0,-1.d0,0.d0/)
      tmpmat(2,:)=(/0.d0,1.d0, 0.d0/)
      tmpmat(3,:)=(/0.d0,0.d0, 1.d0/)
      nigmat(2)=nigmat(2)+ nigmat(1)- nigmat(6)
      nigmat(4)=nigmat(4)-(nigmat(5))
      nigmat(6)=nigmat(6)-(nigmat(1)+nigmat(1))
    else
      tmpmat(1,:)=(/1.d0, 1.d0,0.d0/)
      tmpmat(2,:)=(/0.d0,1.d0, 0.d0/)
      tmpmat(3,:)=(/0.d0,0.d0, 1.d0/)
      nigmat(2)=nigmat(2)+ nigmat(1)+ nigmat(6)
      nigmat(4)=nigmat(4)+(nigmat(5))
      nigmat(6)=nigmat(6)+(nigmat(1)+nigmat(1))
    endif
    if(nigmat(2).lt.0.d0) stop "Error in action a7"
end subroutine


subroutine  a8_action(nigmat,tmpmat,eps)
implicit none
real(8):: nigmat(6),tmpmat(3,3),eps
      tmpmat(1,:)=(/1.d0,0.d0, 1.d0/)
      tmpmat(2,:)=(/0.d0,1.d0, 1.d0/)
      tmpmat(3,:)=(/0.d0,0.d0, 1.d0/)
      nigmat(3)=nigmat(3)+ nigmat(1)+ nigmat(2)+ nigmat(4)+ nigmat(5)+ nigmat(6)
      nigmat(4)=nigmat(4)+ nigmat(2)+ nigmat(2)+ nigmat(6)
      nigmat(5)=nigmat(5)+ nigmat(1)+ nigmat(1)+ nigmat(6)
    if(nigmat(3).lt.0.d0) stop "Error in action a7"
end subroutine

subroutine  n1_action(nigmat,tmpmat,eps)
implicit none
real(8):: nigmat(6),tmpmat(3,3),eps,tmp
      tmpmat(1,:)=(/ 0.d0,-1.d0, 0.d0/)
      tmpmat(2,:)=(/-1.d0, 0.d0, 0.d0/)
      tmpmat(3,:)=(/ 0.d0, 0.d0,-1.d0/)
    tmp=nigmat(1) 
    nigmat(1)=nigmat(2)
    nigmat(2)=tmp
    tmp=nigmat(4)
    nigmat(4)=nigmat(5)
    nigmat(5)=tmp
end subroutine

subroutine  n2_action(nigmat,tmpmat,eps)
implicit none
real(8):: nigmat(6),tmpmat(3,3),eps,tmp
      tmpmat(1,:)=(/-1.d0, 0.d0, 0.d0/)
      tmpmat(2,:)=(/ 0.d0, 0.d0,-1.d0/)
      tmpmat(3,:)=(/ 0.d0,-1.d0, 0.d0/)
    tmp=nigmat(2) 
    nigmat(2)=nigmat(3)
    nigmat(3)=tmp
    tmp=nigmat(5)
    nigmat(5)=nigmat(6)
    nigmat(6)=tmp
end subroutine

subroutine  n3_true_action(nigmat,tmpmat,eps)
implicit none
real(8):: nigmat(6),tmpmat(3,3),eps,tmp
logical:: flt
real(8):: i,j,k
i=1.d0
j=1.d0
k=1.d0
    if (flt(nigmat(4), 0.d0,eps)) i = -1.d0
    if (flt(nigmat(5), 0.d0,eps)) j = -1.d0
    if (flt(nigmat(6), 0.d0,eps)) k = -1.d0
      tmpmat(1,:)=(/ i,   0.d0, 0.d0/)
      tmpmat(2,:)=(/ 0.d0, j  , 0.d0/)
      tmpmat(3,:)=(/ 0.d0,0.d0,    k/)
    nigmat(4)=abs(nigmat(4))
    nigmat(5)=abs(nigmat(5))
    nigmat(6)=abs(nigmat(6))
end subroutine

subroutine  n3_false_action(nigmat,tmpmat,eps)
implicit none
real(8):: nigmat(6),tmpmat(3,3),eps,tmp
logical:: flt,fgt
real(8):: i,j,k,f(3)
integer:: z
    f = 1.d0
    z = -1
    if (fgt(nigmat(4), 0.d0,eps)) then 
        f(1) = -1.d0
    elseif (.not. flt(nigmat(4), 0.d0,eps))then
        z = 1
    endif
    if (fgt(nigmat(5), 0.d0,eps)) then 
        f(2) = -1.d0
    elseif (.not. flt(nigmat(5), 0.d0,eps))then 
        z = 2
    endif
    if (fgt(nigmat(6), 0.d0,eps)) then 
        f(3) = -1.d0
    elseif (.not. flt(nigmat(6), 0.d0,eps))then 
        z = 3
    endif
    if (f(1)*f(2)*f(3) < 0.d0) then
       if(z==-1) stop "Wrong z in n4"
       f(z) = -1.d0
    endif
    tmpmat(1,:)=(/ f(1), 0.d0, 0.d0/)
    tmpmat(2,:)=(/ 0.d0, f(2), 0.d0/)
    tmpmat(3,:)=(/ 0.d0,0.d0,  f(3)/)
    nigmat(4)=-abs(nigmat(4))
    nigmat(5)=-abs(nigmat(5))
    nigmat(6)=-abs(nigmat(6))
end subroutine

function def_gt_0(nigmat,eps)
implicit none
real(8):: nigmat(6),eps
integer:: nzero,npositive,i
logical:: flt,def_gt_0
def_gt_0=.false.
nzero=0
npositive=0
do i=4,6
!if(0.d0.lt.nigmat(i)) then
if(flt(0.d0,nigmat(i),eps)) then
  npositive=npositive+1
  elseif (.not.flt(nigmat(i),0.d0,eps)) then
  nzero=nzero+1
endif
enddo
if(npositive==3.or.(nzero==0.and.npositive==1)) def_gt_0=.true.
end function

function fpos(nigmat,eps)
implicit none
real(8):: nigmat(6),eps
integer:: nzero,npositive,i
logical:: fpos,flt
fpos=.false.
nzero=0
npositive=0
do i=4,6
!if(0.d0.lt.nigmat(i)) then
if(flt(0.d0,nigmat(i),eps)) then
  npositive=npositive+1
  elseif (.not.flt(nigmat(i),0.d0,eps)) then
  nzero=nzero+1
endif
enddo
if(npositive==3.or.(nzero==0.and.npositive==1)) fpos=.true.
end function

function flt(a,b,eps)
implicit none
real(8)::a,b,eps
logical:: flt
if(a.lt.b-eps) then
  flt=.true.
else
  flt=.false.
endif
end function

function fgt(a,b,eps)
implicit none
real(8)::a,b,eps
logical:: fgt,flt
fgt=flt(b,a,eps)
end function

function fle(a,b,eps)
implicit none
real(8)::a,b,eps
logical:: fle
if(.not.(b.lt.a-eps)) then
  fle=.true.
else
  fle=.false.
endif
end function

function fge(a,b,eps)
implicit none
real(8)::a,b,eps
logical:: fge
if(.not.(a.lt.b-eps)) then
  fge=.true.
else
  fge=.false.
endif
end function

function feq(a,b,eps)
implicit none
real(8)::a,b,eps
logical:: feq,flt
feq=(.not. (flt(a,b,eps).or.flt(b, a,eps)))
end function

!subroutine init_nigmat(dist_ang,nigmat,l,m,n,eps,pi)
subroutine init_nigmat(latvec,nigmat,l,m,n,eps,pi)
!Initiallization step (step 0I of the Krivy Gruber algo 
implicit none
real(8):: dist_ang(6),nigmat(6),pi,convang,eps,latvec(3,3)
integer:: l,m,n
convang=pi/180.d0
nigmat(1)=dot_product(latvec(:,1),latvec(:,1))!dist_ang(1)**2
nigmat(2)=dot_product(latvec(:,2),latvec(:,2))!dist_ang(2)**2
nigmat(3)=dot_product(latvec(:,3),latvec(:,3))!dist_ang(3)**2
nigmat(4)=2.d0*dot_product(latvec(:,2),latvec(:,3))!2.d0*dist_ang(2)*dist_ang(3)*cos(convang*dist_ang(4))
nigmat(5)=2.d0*dot_product(latvec(:,3),latvec(:,1))!2.d0*dist_ang(3)*dist_ang(1)*cos(convang*dist_ang(5))
nigmat(6)=2.d0*dot_product(latvec(:,1),latvec(:,2))!2.d0*dist_ang(1)*dist_ang(2)*cos(convang*dist_ang(6))
l=0;m=0;n=0
end subroutine

subroutine init_lmn(nigmat,l,m,n,eps)
implicit none
real(8):: dist_ang(6),nigmat(6),eps
integer:: l,m,n
l=0
m=0
n=0
if(nigmat(4).lt.-eps) l=-1
if(nigmat(4).gt. eps) l= 1
if(nigmat(5).lt.-eps) m=-1
if(nigmat(5).gt. eps) m= 1
if(nigmat(6).lt.-eps) n=-1
if(nigmat(6).gt. eps) n= 1
end subroutine


!subroutine dist_ang2latvec(dist_ang,latvec,pi)
!!This subroutine will generate the lattice vector representation of the cell
!!from the length/angle representation
!implicit none
!real(8):: dist_ang(6),latvec(3,3),convang,pi
!convang=pi/180.d0
!
!latvec(1,1)=dist_ang(1)
!latvec(2,1)=0.d0
!latvec(3,1)=0.d0
!latvec(1,2)=dist_ang(2)*cos(convang*dist_ang(6))
!latvec(2,2)=dist_ang(2)*sin(convang*dist_ang(6))
!latvec(3,2)=0.d0
!latvec(1,3)=dist_ang(3)*cos(convang*dist_ang(5))
!latvec(2,3)=(dist_ang(2)*dist_ang(3)*cos(convang*dist_ang(4))-latvec(1,2)*latvec(1,3))/latvec(2,2)
!latvec(3,3)=(dist_ang(3)**2-latvec(1,3)**2-latvec(2,3)**2)
!if(latvec(3,3).lt.0.d0) stop "Bad angles in dist_ang2latvec"
!latvec(3,3)=sqrt(latvec(3,3))
!end subroutine

subroutine latvec2dist_ang(dist_ang,latvec,pi)
!This subroutine will generate the distance and angles of a cell starting from latvec
implicit none
real(8):: dist_ang(6),latvec(3,3),pi,convang
convang=180.d0/pi
dist_ang(1)=sqrt(dot_product(latvec(:,1),latvec(:,1)))
dist_ang(2)=sqrt(dot_product(latvec(:,2),latvec(:,2)))
dist_ang(3)=sqrt(dot_product(latvec(:,3),latvec(:,3)))
dist_ang(4)=convang*acos(dot_product(latvec(:,2),latvec(:,3))/(dist_ang(2)*dist_ang(3)))
dist_ang(5)=convang*acos(dot_product(latvec(:,3),latvec(:,1))/(dist_ang(3)*dist_ang(1)))
dist_ang(6)=convang*acos(dot_product(latvec(:,1),latvec(:,2))/(dist_ang(1)*dist_ang(2)))
end subroutine

function is_niggli_cell(nigmat,eps)
implicit none
real(8):: nigmat(6),eps
logical::is_niggli_cell,is_buerger_cell,feq,fgt
is_niggli_cell=.true.
    if (.not.is_buerger_cell(nigmat,eps)) then
        is_niggli_cell=.false.
        return
    endif
    if (feq(nigmat(4),nigmat(2),eps)) then
      if (fgt(nigmat(6),nigmat(5)+nigmat(5),eps)) then 
        is_niggli_cell=.false.
        return
      endif
    endif
    if (feq(nigmat(5),nigmat(1),eps)) then
      if (fgt(nigmat(6),nigmat(4)+nigmat(4),eps)) then 
        is_niggli_cell=.false.
        return
      endif
    endif
    if (feq(nigmat(6),nigmat(1),eps)) then
      if (fgt(nigmat(5),nigmat(4)+nigmat(4),eps)) then 
        is_niggli_cell=.false.
        return
      endif
    endif
    if (feq(nigmat(4),-nigmat(2),eps)) then
      if (.not.feq(nigmat(6),0.d0,eps)) then 
        is_niggli_cell=.false.
        return
      endif
    endif
    if (feq(nigmat(5),-nigmat(1),eps)) then
      if (.not.feq(nigmat(6),0.d0,eps)) then 
        is_niggli_cell=.false.
        return
      endif
    endif
    if (feq(nigmat(6),-nigmat(1),eps)) then
      if (.not.feq(nigmat(5),0.d0,eps)) then 
        is_niggli_cell=.false.
        return
      endif
    endif
    if (feq(nigmat(4)+nigmat(5)+nigmat(6)+nigmat(1)+nigmat(2),0.d0,eps)) then
      if (fgt(nigmat(1)+nigmat(1)+nigmat(5)+nigmat(5)+nigmat(6),0.d0,eps)) then 
        is_niggli_cell=.false.
        return
      endif
    endif
end function

function  meets_primary_conditions(nigmat,eps)
implicit none
real(8):: nigmat(6),eps
logical::meets_primary_conditions,fgt
meets_primary_conditions=.true.
    if (fgt(nigmat(1), nigmat(2),eps)) then
        meets_primary_conditions=.false.
        return
    endif
    if (fgt(nigmat(2), nigmat(3),eps)) then
        meets_primary_conditions=.false.
        return
    endif
    if (fgt(abs(nigmat(4)), nigmat(2),eps)) then
        meets_primary_conditions=.false.
        return
    endif
    if (fgt(abs(nigmat(5)), nigmat(1),eps)) then
        meets_primary_conditions=.false.
        return
    endif
    if (fgt(abs(nigmat(6)), nigmat(1),eps)) then
        meets_primary_conditions=.false.
        return
    endif
end function

function  meets_main_conditions(nigmat,eps)
implicit none
integer:: typer
real(8):: nigmat(6),eps
logical::meets_main_conditions,meets_primary_conditions,flt
meets_main_conditions=.true.
    if (.not. meets_primary_conditions(nigmat,eps)) then
      meets_main_conditions=.false.
      return 
    endif
    
    if (typer(nigmat,eps) == 0) then 
      meets_main_conditions=.false.
      return
    endif
    if (typer(nigmat,eps) == 2) then
      if (flt(nigmat(4)+nigmat(5)+nigmat(6)+nigmat(1)+nigmat(2), 0.d0,eps)) then
          meets_main_conditions=.false.
          return
      endif
    endif
end function

function is_buerger_cell(nigmat,eps)
implicit none
real(8):: nigmat(6),eps
logical::is_buerger_cell,meets_main_conditions,feq,fgt

    if (.not. meets_main_conditions(nigmat,eps)) then
        is_buerger_cell=.false.
        return
    endif
    if (feq(nigmat(1),nigmat(2),eps)) then
      if (fgt(abs(nigmat(4)),abs(nigmat(5)),eps)) then
          is_buerger_cell=.false.
          return
      endif
    endif
    if (feq(nigmat(2), nigmat(3),eps)) then
      if (fgt(abs(nigmat(5)), abs(nigmat(6)),eps))then 
          is_buerger_cell=.false.
          return
      endif
   endif
  is_buerger_cell=.true. 

end function

function typer(nigmat,eps)
implicit none
real(8):: nigmat(6),eps
integer:: nzero,npositive,i,typer
logical:: flt
typer=0
nzero=0
npositive=0
do i=4,6
!if(0.d0.lt.nigmat(i)) then
if(flt(0.d0,nigmat(i),eps)) then
  npositive=npositive+1
  elseif (.not.flt(nigmat(i),0.d0,eps)) then
  nzero=nzero+1
endif
enddo
    typer=0
    if (npositive == 3) typer=1
    if (npositive == 0) typer=2
end function

 subroutine invertmat(mat,matinv,n)
 implicit none
 real(8),intent(in) :: mat(n,n)
 integer               :: n
 real(8),allocatable   :: WORK(:)
 real(8)               :: matinv(n,n),det(3),a(n,n),div
 integer               :: IPIV(n), INFO
 integer               :: LDWORK
 !Here only for a 3*3 matrix
 if (n==3) then
 a=mat
 div=(a(1,1)*a(2,2)*a(3,3)-a(1,1)*a(2,3)*a(3,2)-a(1,2)*a(2,1)*a(3,3)+a(1,2)*a(2,3)*a(3,1)+&
 &a(1,3)*a(2,1)*a(3,2)-a(1,3)*a(2,2)*a(3,1)) 
 div=1.d0/div
      matinv(1,1) = (a(2,2)*a(3,3)-a(2,3)*a(3,2))*div
      matinv(1,2) =-(a(1,2)*a(3,3)-a(1,3)*a(3,2))*div
      matinv(1,3) = (a(1,2)*a(2,3)-a(1,3)*a(2,2))*div
      matinv(2,1) =-(a(2,1)*a(3,3)-a(2,3)*a(3,1))*div
      matinv(2,2) = (a(1,1)*a(3,3)-a(1,3)*a(3,1))*div
      matinv(2,3) =-(a(1,1)*a(2,3)-a(1,3)*a(2,1))*div
      matinv(3,1) = (a(2,1)*a(3,2)-a(2,2)*a(3,1))*div
      matinv(3,2) =-(a(1,1)*a(3,2)-a(1,2)*a(3,1))*div
      matinv(3,3) = (a(1,1)*a(2,2)-a(1,2)*a(2,1))*div
 else
          stop "ONLY 3x3 MATRIX INVERSION"
 endif
 end
!**********************************************************************************************
subroutine find_kpt(k1,k2,k3,lat,gridden)
!This code will define the KPT mesh based on the desired grid density
  implicit none
  integer       :: k1,k2,k3,i,j
  real(8)       :: lat1, lat2, lat3
  real(8)       :: angles(3),cos_arr(3)
  real(8)       :: lat(3,3),glat(3,3),crossp(3),vol,gridden,glen(3),a(3,3)
  real(8)       :: pi
  pi=acos(-1.d0)

  lat1=sqrt(lat(1,1)**2+lat(2,1)**2+lat(3,1)**2)
  lat2=sqrt(lat(1,2)**2+lat(2,2)**2+lat(3,2)**2)
  lat3=sqrt(lat(1,3)**2+lat(2,3)**2+lat(3,3)**2)
  cos_arr(1)=(lat(1,2)*lat(1,3)+lat(2,2)*lat(2,3)+lat(3,2)*lat(3,3))/lat2/lat3
  cos_arr(2)=(lat(1,1)*lat(1,3)+lat(2,1)*lat(2,3)+lat(3,1)*lat(3,3))/lat1/lat3
  cos_arr(3)=(lat(1,1)*lat(1,2)+lat(2,1)*lat(2,2)+lat(3,1)*lat(3,2))/lat1/lat2
  angles(1)=(acos(cos_arr(1))/pi)*180.d0
  angles(2)=(acos(cos_arr(2))/pi)*180.d0
  angles(3)=(acos(cos_arr(3))/pi)*180.d0

  a=lat
  vol= a(1,1)*a(2,2)*a(3,3)-a(1,1)*a(2,3)*a(3,2)-a(1,2)*a(2,1)*a(3,3)+&
       a(1,2)*a(2,3)*a(3,1)+a(1,3)*a(2,1)*a(3,2)-a(1,3)*a(2,2)*a(3,1)
  call cross_product(lat(:,2),lat(:,3),crossp(:))
  glat(:,1)=2.d0*pi*crossp(:)/vol
  call cross_product(lat(:,3),lat(:,1),crossp(:))
  glat(:,2)=2.d0*pi*crossp(:)/vol
  call cross_product(lat(:,1),lat(:,2),crossp(:))
  glat(:,3)=2.d0*pi*crossp(:)/vol
!Compute the correct kpts
  glen(1)=sqrt(glat(1,1)**2+glat(2,1)**2+glat(3,1)**2)
  glen(2)=sqrt(glat(1,2)**2+glat(2,2)**2+glat(3,2)**2)
  glen(3)=sqrt(glat(1,3)**2+glat(2,3)**2+glat(3,3)**2)
  call track_kpt(gridden,glen(1),k1)
  call track_kpt(gridden,glen(2),k2)
  call track_kpt(gridden,glen(3),k3)
end subroutine find_kpt

subroutine track_kpt(gridden,glen,kpt)
  implicit none
  integer :: kpt,j
  real(8) :: glen,gridden,d_test,pi 
  pi=acos(-1.d0)
  kpt=int(glen/(gridden*2.d0*pi))
  if (kpt == 0) kpt = 1
  d_test=glen/(kpt*2.d0*pi)
  if (d_test.ge.gridden) then
     do j=1,25
        kpt=kpt+j
        d_test=glen/(kpt*2.d0*pi)
        if (d_test.le.gridden) exit
     enddo
  endif
end subroutine track_kpt
!************************************************************************************
 subroutine cross_product(a,b,crossp)
 !a very simple implementation of the cross product
 implicit none
 real(8)::a(3),b(3)
 real(8)::crossp(3)
 crossp(1)=a(2)*b(3)-a(3)*b(2)
 crossp(2)=a(3)*b(1)-a(1)*b(3)
 crossp(3)=a(1)*b(2)-a(2)*b(1)
 return
 end subroutine

!************************************************************************************

 subroutine rxyz_cart2int(latvec,rxyzint,rxyzcart,nat)
 !This subrouine will convert the internal coordinates into cartesian coordinates
 implicit none
 real(8):: rxyzint(3,nat), rxyzcart(3,nat),latvec(3,3),latvecinv(3,3)
 integer:: nat,iat
 call invertmat(latvec,latvecinv,3)
 do iat=1,nat
  rxyzint(:,iat)=matmul(latvecinv,rxyzcart(:,iat))
 enddo
 end subroutine rxyz_cart2int

!************************************************************************************
