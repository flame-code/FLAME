!*****************************************************************************************
subroutine writexyz(filename,fn_position,nat,rat,bemoved,sat,cellvec,boundcond,comment)
    !use charges, only:nzx
    !use interactions, only: imass
    !use mod_atoms, only:
    implicit none
    integer, intent(in):: nat
    real(8), intent(in):: rat(3,nat)
    logical, intent(in):: bemoved(3,nat)
    character(5), intent(in):: sat(nat)
    real(8), intent(in):: cellvec(3,3)
    character(*), intent(in):: filename, fn_position, boundcond, comment
    !local variables
    integer:: iat, ios
    real(8):: x, y, z !, cellx, celly, cellz
    if(trim(fn_position)=='new') then
        open(unit=1358, file=filename,status='replace',iostat=ios)
    elseif(trim(fn_position)=='append') then
        open(unit=1358, file=filename,status='old',position='append',iostat=ios)
    else
        write(*,'(2a)') 'ERROR: fn_position is unknown, ',trim(fn_position)
    endif
    if(ios/=0) then;write(*,'(2a)') 'ERROR: failure openning ',trim(filename);stop;endif
    write(1358,'(i6,2x,a)')   nat,trim(comment)
    !write(1358,'(a)') '   1.66728000E+01   1.66728000E+01   1.66728000E+01'
    !cellx=cellvec(1,1)*0.529177d0
    !celly=cellvec(2,2)*0.529177d0
    !cellz=cellvec(3,3)*0.529177d0
    write(1358,'(a,2x,6(e24.15))') trim(boundcond),cellvec(1,1),cellvec(2,2),cellvec(3,3),cellvec(1,2),cellvec(1,3),cellvec(2,3)
    do iat=1,nat
        !if(nzx(imass(iat))==1 ) an=' H'
        !if(nzx(imass(iat))==14) an='Si'
        !an=' C '
        x=rat(1,iat) !*0.529177d0
        y=rat(2,iat) !*0.529177d0
        z=rat(3,iat) !*0.529177d0
        write(1358,'(a5,2x,3es24.15,2x,3l1)') trim(sat(iat)),x,y,z, &
            bemoved(1,iat),bemoved(2,iat),bemoved(3,iat)
    enddo
    close(1358)
end subroutine writexyz
!*****************************************************************************************
subroutine readxyznat(filename,nat)
    implicit none
    integer:: nat
    character(*):: filename
    !local variables
    integer:: ios
    open(unit=1358, file=trim(filename),status='old',iostat=ios)
    if(ios/=0) then;write(*,'(2a)') 'ERROR: failure openning ',trim(filename);stop;endif
    read(1358,*) nat
    close(1358)
end subroutine readxyznat
!*****************************************************************************************
subroutine readxyz(filename,nat,rat,sat,comment1,comment2,atom_motion)
    implicit none
    integer:: nat
    real(8):: rat(3,nat)
    character(5):: sat(nat)
    character(*):: filename, comment1, comment2
    logical:: atom_motion(3,nat)
    !local variables
    integer:: iat, ios
    character(256):: str1, str2
    character(3):: str_motion
    open(unit=1358, file=trim(filename),status='old',iostat=ios)
    if(ios/=0) then;write(*,'(2a)') 'ERROR: failure openning ',trim(filename);stop;endif
    read(1358,'(a)') str1
    read(1358,'(a)') str2
    comment1=str1
    comment2=str2
    do iat=1,nat
        read(1358,*) sat(iat),rat(1,iat),rat(2,iat),rat(3,iat),str_motion
        if(str_motion(1:1)=='T') then
            atom_motion(1,iat)=.true.
        elseif(str_motion(1:1)=='F') then
            atom_motion(1,iat)=.false.
        else
            stop 'ERROR: incorrect symbol for atomic motion in input configuration'
        endif
        if(str_motion(2:2)=='T') then
            atom_motion(2,iat)=.true.
        elseif(str_motion(2:2)=='F') then
            atom_motion(2,iat)=.false.
        else
            stop 'ERROR: incorrect symbol for atomic motion in input configuration'
        endif
        if(str_motion(3:3)=='T') then
            atom_motion(3,iat)=.true.
        elseif(str_motion(3:3)=='F') then
            atom_motion(3,iat)=.false.
        else
            stop 'ERROR: incorrect symbol for atomic motion in input configuration'
        endif
    enddo
    close(1358)
end subroutine readxyz
!*****************************************************************************************
