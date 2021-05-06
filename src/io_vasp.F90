!*****************************************************************************************
!This routine will write the file "filename" in ascii file format
subroutine write_poscar(filename,nat,rat,latvec,ntypat,natarr,comment,vasp5,comment2,atom_motion)
    implicit none
    integer:: nat, ntypat, natarr(128)
    character(*):: filename
    real(8):: rat(3,nat),latvec(3,3)
    character(*):: comment, comment2
    logical:: vasp5, atom_motion(3,nat)
    !local variables
    integer:: iat, ioserr, i
    !So if units==angstroem, the file will be converted to angstroem
    !   if units==bohr, the positions will not be changed
    !units="angstroem"
    !-------------------------------------------------------
    open(unit=46,file=trim(filename),status='replace',iostat=ioserr)
    if(ioserr/=0) then
        write(*,'(3a)') 'ERROR: openning file ',trim(filename),' failed.'
        stop
    endif
    write(46,'(a)') trim(comment)
    write(46,'(f20.10)') 1.d0
    !write(46,'(3es20.10)') latvec(1,1),latvec(2,1),latvec(3,1)
    !write(46,'(3es20.10)') latvec(1,2),latvec(2,2),latvec(3,2)
    !write(46,'(3es20.10)') latvec(1,3),latvec(2,3),latvec(3,3)
    write(46,'(1x,3f22.16)') latvec(1,1),latvec(2,1),latvec(3,1)
    write(46,'(1x,3f22.16)') latvec(1,2),latvec(2,2),latvec(3,2)
    write(46,'(1x,3f22.16)') latvec(1,3),latvec(2,3),latvec(3,3)
    if(vasp5) then
        write(46,'(a)') trim(comment2)
    endif
    do i=1,ntypat
        write(46,'(i4)',advance='no') natarr(i)
    enddo
    write(46,*)
    write(46,'(a)') 'selective dynamics'
    !write(46,'(a)') 'direct'
    write(46,'(a)') 'cartesian'
    do iat=1,nat
        write(46,'(3es20.10,3l3)') rat(1,iat),rat(2,iat),rat(3,iat), &
            atom_motion(1,iat),atom_motion(2,iat),atom_motion(3,iat)
    enddo
    close(46)
end subroutine write_poscar
!*****************************************************************************************
