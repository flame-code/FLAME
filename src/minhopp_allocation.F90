!*****************************************************************************************
subroutine allocate_minhopp_arrays1(nproc)
    use mod_minhopp, only: nlmin, nlminx, earr, nvisit, dtarr, ediffarr, ekinarr, &
        itagintermediate, mtagarr1, mtagarr2, do_req1, do_req2, ireqarr1, ireqarr2, nbuf
    implicit none
    integer, intent(in):: nproc
    !local variables
    integer:: istat
    if(nlmin==-1) stop 'ERROR: nlmin=-1, is not initialized.'
    if(nlminx==-1) stop 'ERROR: nlminx=-1, is not initialized.'
    allocate(earr(0:nlminx+nbuf),stat=istat)
    if(istat/=0) write(*,'(a)') 'ERROR: failure in allocating array earr of mod_minhopp.'
    allocate(nvisit(0:nlminx+nbuf),stat=istat)
    if(istat/=0) write(*,'(a)') 'ERROR: failure in allocating array nvisit of mod_minhopp.'
    allocate(dtarr(0:nproc-1),stat=istat)
    if(istat/=0) write(*,'(a)') 'ERROR: failure in allocating array dtarr of mod_minhopp.'
    allocate(ediffarr(0:nproc-1),stat=istat)
    if(istat/=0) write(*,'(a)') 'ERROR: failure in allocating array ediffarr of mod_minhopp.'
    allocate(ekinarr(0:nproc-1),stat=istat)
    if(istat/=0) write(*,'(a)') 'ERROR: failure in allocating array ekinarr of mod_minhopp.'
    allocate(itagintermediate(0:nproc-1),stat=istat)
    if(istat/=0) write(*,'(a)') 'ERROR: failure in allocating array itagintermediate of mod_minhopp.'
    allocate(mtagarr1(0:nproc-1),stat=istat)
    if(istat/=0) write(*,'(a)') 'ERROR: failure in allocating array mtagarr1 of mod_minhopp.'
    allocate(mtagarr2(0:nproc-1),stat=istat)
    if(istat/=0) write(*,'(a)') 'ERROR: failure in allocating array mtagarr2 of mod_minhopp.'
    allocate(do_req1(0:nproc-1),stat=istat)
    if(istat/=0) write(*,'(a)') 'ERROR: failure in allocating array do_req1 of mod_minhopp.'
    allocate(do_req2(0:nproc-1),stat=istat)
    if(istat/=0) write(*,'(a)') 'ERROR: failure in allocating array do_req2 of mod_minhopp.'
    allocate(ireqarr1(0:nproc-1),stat=istat)
    if(istat/=0) write(*,'(a)') 'ERROR: failure in allocating array ireqarr1 of mod_minhopp.'
    allocate(ireqarr2(0:nproc-1),stat=istat)
    if(istat/=0) write(*,'(a)') 'ERROR: failure in allocating array ireqarr2 of mod_minhopp.'
    itagintermediate(0:nproc-1)=1
    mtagarr1(0:nproc-1)=0
    mtagarr2(0:nproc-1)=1
    nvisit(0:nlminx+nbuf)=0
    do_req1(0:nproc-1)=.true.
    do_req2(0:nproc-1)=.true.
end subroutine allocate_minhopp_arrays1
!*****************************************************************************************
subroutine allocate_minhopp_arrays2(nat,nproc)
    use mod_minhopp, only: nbuf, abuf, abufall, abuf1, abuf2
    implicit none
    integer, intent(in):: nat, nproc
    !local variables
    integer:: istat
    allocate(abuf(3*nat+4),stat=istat)
    if(istat/=0) write(*,'(a)') 'ERROR: failure in allocating array abuf of mod_minhopp.'
    allocate(abufall(3*nat+4,0:nbuf-1),stat=istat)
    if(istat/=0) write(*,'(a)') 'ERROR: failure in allocating array abufall of mod_minhopp.'
    allocate(abuf1(3*nat+1,0:nproc-1),stat=istat)
    if(istat/=0) write(*,'(a)') 'ERROR: failure in allocating array abuf1 of mod_minhopp.'
    allocate(abuf2(3*nat+4,0:nproc-1),stat=istat)
    if(istat/=0) write(*,'(a)') 'ERROR: failure in allocating array abuf2 of mod_minhopp.'
end subroutine allocate_minhopp_arrays2
!*****************************************************************************************
subroutine deallocate_minhopp_arrays
    use mod_minhopp, only: earr, nvisit, abuf, abufall, dtarr, ediffarr, ekinarr, &
        itagintermediate, mtagarr1, mtagarr2, abuf1, abuf2, do_req1, do_req2, ireqarr1, ireqarr2
    implicit none
    !local variables
    integer:: istat
    deallocate(earr,stat=istat)
    if(istat/=0) write(*,'(a)') 'ERROR: failure in deallocating array earr of mod_minhopp.'
    deallocate(nvisit,stat=istat)
    if(istat/=0) write(*,'(a)') 'ERROR: failure in deallocating array nvisit of mod_minhopp.'
    deallocate(abuf,stat=istat)
    if(istat/=0) write(*,'(a)') 'ERROR: failure in deallocating array abuf of mod_minhopp.'
    deallocate(abufall,stat=istat)
    if(istat/=0) write(*,'(a)') 'ERROR: failure in deallocating array abufall of mod_minhopp.'
    deallocate(dtarr,stat=istat)
    if(istat/=0) write(*,'(a)') 'ERROR: failure in deallocating array dtarr of mod_minhopp.'
    deallocate(ediffarr,stat=istat)
    if(istat/=0) write(*,'(a)') 'ERROR: failure in deallocating array ediffarr of mod_minhopp.'
    deallocate(ekinarr,stat=istat)
    if(istat/=0) write(*,'(a)') 'ERROR: failure in deallocating array ekinarr of mod_minhopp.'
    deallocate(itagintermediate,stat=istat)
    if(istat/=0) write(*,'(a)') 'ERROR: failure in deallocating array itagintermediate of mod_minhopp.'
    deallocate(mtagarr1,stat=istat)
    if(istat/=0) write(*,'(a)') 'ERROR: failure in deallocating array mtagarr1 of mod_minhopp.'
    deallocate(mtagarr2,stat=istat)
    if(istat/=0) write(*,'(a)') 'ERROR: failure in deallocating array mtagarr2 of mod_minhopp.'
    deallocate(abuf1,stat=istat)
    if(istat/=0) write(*,'(a)') 'ERROR: failure in deallocating array abuf1 of mod_minhopp.'
    deallocate(abuf2,stat=istat)
    if(istat/=0) write(*,'(a)') 'ERROR: failure in deallocating array abuf2 of mod_minhopp.'
    deallocate(do_req1,stat=istat)
    if(istat/=0) write(*,'(a)') 'ERROR: failure in deallocating array do_req1 of mod_minhopp.'
    deallocate(do_req2,stat=istat)
    if(istat/=0) write(*,'(a)') 'ERROR: failure in deallocating array do_req2 of mod_minhopp.'
    deallocate(ireqarr1,stat=istat)
    if(istat/=0) write(*,'(a)') 'ERROR: failure in deallocating array ireqarr1 of mod_minhopp.'
    deallocate(ireqarr2,stat=istat)
    if(istat/=0) write(*,'(a)') 'ERROR: failure in deallocating array ireqarr2 of mod_minhopp.'
end subroutine deallocate_minhopp_arrays
!*****************************************************************************************
