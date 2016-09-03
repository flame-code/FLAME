!> @file
!! Include fortran file for deallocation template
!! @author
!!    Copyright (C) 2012-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS

  !      call timing(0,'AllocationProf','IR') 

  !here the size should be corrected with ndebug (or maybe not)
  ilsize=int(product(shape(array))*kind(array),kind=8)
  !fortran deallocation
  deallocate(array,stat=ierror)

  if (f_err_raise(ierror/=0,&
       'Deallocation problem, error code '//trim(yaml_toa(ierror)),ERR_DEALLOCATE)) return

  !profile address, in case of profiling activated
!  if (m%profile) then 
     !address of first element (not needed for deallocation)
     !call getaddress(array,address,len(address),ierr)
     !address of the metadata 
     call metadata_address(len(array),array,iadd)
     address=repeat(' ',len(address))
     address=trim(long_toa(iadd))

     !hopefully only address is necessary for the deallocation

     !search in the dictionaries the address
     !a error event should be raised in this case
     dict_add=>find_key(mems(ictrl)%dict_routine,trim(address))
     if (.not. associated(dict_add)) then
        dict_add=>find_key(mems(ictrl)%dict_global,trim(address))
        if (f_err_raise(.not. associated(dict_add),'address'//trim(address)//&
             'not present in dictionary',ERR_INVALID_MALLOC)) then
           return
        else
           use_global=.true.
        end if
     else
        use_global=.false.
     end if

     array_id=dict_add//arrayid
     routine_id=dict_add//routineid
     jlsize=dict_add//sizeid
     if (f_err_raise(ilsize /= jlsize,'Size of array '//trim(array_id)//&
          ' ('//trim(yaml_toa(ilsize))//') not coherent with dictionary, found='//&
          trim(yaml_toa(jlsize)),ERR_MALLOC_INTERNAL)) return

     call memocc(ierror,-int(ilsize),trim(array_id),trim(routine_id))

     if (use_global) then
        !call yaml_dict_dump(dict_global)
        call pop(mems(ictrl)%dict_global,trim(address))
     else
        call pop(mems(ictrl)%dict_routine,trim(address))
     end if
!  end if

  !      call timing(0,'AllocationProf','RS') 
