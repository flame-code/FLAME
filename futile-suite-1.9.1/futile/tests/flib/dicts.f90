!> @file
!! Test the dictionaries of flib
!! @example dicts.f90
!! Some examples about dictionaries
!! @author
!!    Copyright (C) 2013-2015 BigDFT group <br>
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Routine testing the dictionary object of flib
subroutine test_dictionaries0()
  use yaml_output
  use dictionaries
  use exception_callbacks
  implicit none
  type(dictionary), pointer :: dict1,dict2,dict3
  !local variables
  integer :: ival!,nval
!!$  character(len=2) :: val
  character(len=30) :: val2
  type(dictionary), pointer :: dict_tmp
  !finding operations
!!$  print *,' Filling a linked list'
!!$  call dict_init(list)
!!$
!!$  call add(list,'x1')
!!$  call add(list,'x2')
!!$  call add(list,'y1')
!!$  call add(list,'z1')
!!$  call add(list,'z1')
!!$
!!$  call yaml_dict_dump(list)
!!$
!!$  nval=dict_len(list)
!!$  print *,' Number of elements',dict_len(list)
!!$
!!$  do ival=0,nval-1
!!$     val=list//ival
!!$     print *,'value',ival,valget
!!$  end do
!!$
!!$  call dict_free(list)

  call yaml_comment('Now we test dictionaries inside yaml.')
  !Test a dictionary
  !alternative way of initializing a dictionary
  !call dict_init(dict1)

!!! [Creation]
  !example of dict_new usage1
  dict1=>dict_new()
  !end of example of dict_new usage1
!!! [Creation]
  call f_err_open_try()
  ival=dict1//'Toto'

  call yaml_map('ival not existing, fake value',ival)

  call yaml_map('An error has been raised',f_err_check())
  call yaml_map('Its error id is',f_get_last_error())
  !routine to retrieve the error
  call f_dump_last_error()
  call f_err_close_try()
  call yaml_map('Error pipe is still full',f_err_check())

  ! a single scalar
!!  call set(dict1//'',1)
!!$  !can be also set like that, should be avoided
!  call set(dict1,1
!  call yaml_dict_dump(dict1)
!stop
  call set(dict1//'toto',1)
!stop
  call set(dict1//'titi',1.d0)
  call set(dict1//'tutu',(/ '1', '2' /))
  call dict_init(dict2)
  call set(dict2//'a',0)

  !this had  a bug, now solved
  call set(dict1//'List',list_new((/.item. dict2,.item. '4',.item. '1.0'/)))

  !search for a list element
  call yaml_map('1.0 index',dict1//'List' .index. '1.0')

  !example of a list
  dict_tmp => list_new([.item. 'one',.item. '4',.item. '1.1'])
  !end of example of a list
  call yaml_map('1.1 index',dict_tmp .index. '1.1')
  call dict_free(dict_tmp)
  nullify(dict_tmp)


  !this works
!!$  call add(dict1//'List',dict2)
!!$  call add(dict1//'List',4)
!!$  call add(dict1//'List',1.0)

!!$  !this works too
!!$  list=>dict_new()
!!$  call add(list,dict2)
!!$  call add(list,4)
!!$  call add(list,1.0)
!!$  call set(dict1//'List',list)

!!$  !this also
!!$  list=>dict_new()
!!$  call set(list//'First',dict2)
!!$  call set(list//'Second',4)
!!$  call set(list//'Third',1.0)
!!$  call set(dict1//'List',list)

!!$  dict3=>dict1//'List'
!!$  call yaml_map('Elements of the new dictionary (elems, list)',&
!!$       (/dict3%data%nelems,dict3%data%nitems/))

  dict3=> dict1//'New key'
  call set(dict3//'Example',4)
  call yaml_dict_dump(dict1,flow=.true.)

  !test length functions of dictionaries
  call yaml_map('List length',dict_len(dict1//'List'))
  call yaml_map('Dictionary size',dict_size(dict1))
  call dict_free(dict1)

  !new test, build dictionary on-the-fly
  dict1=>dict_new(['Key1' .is. 'One',&
       'Key2' .is. 'Two','Key3' .is. 'Three'])
  !end of new test, build dictionary

  call yaml_dict_dump(dict1)
  call dict_free(dict1)

  dict1=>dict_new()

  call set(dict1//'hgrid',0.5,fmt='(1pe17.5)')
  call yaml_map('Length and size before',(/dict_len(dict1//'hgrid'),dict_size(dict1//'hgrid')/))
  !call add(dict1//'hgrid','new')
  call set(dict1//'hgrid'//0,'new')

  call yaml_mapping_open('There was a hidden problem here')
  call yaml_dict_dump(dict1)
  call yaml_mapping_close()

  call yaml_map('Value of dict1//hgrid',trim(dict_value(dict1//'hgrid')))

  !retrieve value
  val2=dict1//'hgrid' !dict_value(dict1//'hgrid')
  call yaml_map('Value retrieved with equal sign',trim(val2))

  !test value of the dictionary, explicitly
  dict_tmp=>dict1//'hgrid'
  call yaml_map('Value explicitly written in the dictionary',trim(dict_tmp%data%value))

  !test length and sizes of the dictionary
  call yaml_map('Length and size after',(/dict_len(dict_tmp),dict_size(dict_tmp)/))

  call dict_free(dict1)

!stop
  dict1=>dict_new()
  call set(dict1//'hgrid',dict_new((/'test1' .is. '1','test2' .is. '2'/)))

  !search for a dictionary item
  call yaml_map('test2 index',dict1//'hgrid' .index. 'test2')
  call yaml_map('hgrid index',dict1 .index. 'hgrid')
  call yaml_comment('Improper testing of index function',hfill='TEST')

  call yaml_map('Length and size before',(/dict_len(dict1//'hgrid'),dict_size(dict1//'hgrid')/))
  call set(dict1//'hgrid'//0,'new')

  call yaml_mapping_open('Hidden problem here')
  call yaml_dict_dump(dict1)
  call yaml_mapping_close()

  call yaml_map('Value of dict1//hgrid',trim(dict_value(dict1//'hgrid')))

  !retrieve value
  val2=dict1//'hgrid' !dict_value(dict1//'hgrid')
  call yaml_map('Value retrieved with equal sign',trim(val2))

  !test value of the dictionary, explicitly
  dict_tmp=>dict1//'hgrid'
  call yaml_map('Verify that the child is still associated',associated(dict_tmp%child))

  !test length and sizes of the dictionary
  call yaml_map('Length and size after',(/dict_len(dict_tmp),dict_size(dict_tmp)/))

  call dict_free(dict1)

  !test the values of the len after reaffectations
  call dict_init(dict1)
  call set(dict1//'key','scalar0')
  call yaml_map('Entered dict',dict1)
  call yaml_map('Length',dict_len(dict1//'key'))
  call set(dict1//'key',['one','two','thr'])
  call yaml_map('Entered dict',dict1)
  call yaml_map('Length',dict_len(dict1//'key'))
  !reaffectation
  call set(dict1//'key','scalar')
  call yaml_map('Entered dict',dict1)
  call yaml_map('Length',dict_len(dict1//'key'))
  !reaffectation again
  call set(dict1//'key',['One','Two','Thr'])
  call yaml_map('Entered dict',dict1)
  call yaml_map('Length',dict_len(dict1//'key'))

  call dict_free(dict1)

!!$
!!$  !new test, build list on-the-fly
!!$  dict1=list_new((/ .item. 'Val1', .item. 'Val2', .item. 'Val3' ,&
!!$       .item. 'Val4'/))
!!$  call yaml_dict_dump(dict1)
!!$  call dict_free(dict1)

end subroutine test_dictionaries0


subroutine test_dictionaries1()
  use f_precisions
  use yaml_output
  use yaml_strings
  use dictionaries
  use yaml_parse
  implicit none
  !local variables
   integer :: ival,i
!   integer :: j
   type(dictionary), pointer :: dict2
   type(dictionary), pointer :: dict,dictA
   type(dictionary), pointer :: dictA2,dict_tmp,zero1,zero2
   real(f_double), dimension(3) :: tmp_arr
   real(f_double), dimension(3,3) :: tmp_mat,mat
   character(len=20) :: name

   !testing add
   call dict_init(dict)
!   call set(dict//0,1)
!   call set(dict//1,2)
!   call set(dict//2,3)
   call add(dict,'1')
   call add(dict,'2')
   call add(dict,'3')
   call yaml_mapping_open('List')
   call yaml_dict_dump(dict,flow=.true.)
   call yaml_mapping_close()
!after this call the document has to finish
   call yaml_release_document()

   call yaml_new_document()


   call yaml_map('Dictionary length',dict_len(dict))
   call yaml_map('Dictionary size',dict_size(dict))

   call dict_free(dict)

   call yaml_comment('Fortran Dictionary Test',hfill='~')

   call dict_init(dict)

   !Normal filling of the dictionary
   !this fills a last level
   call set(dict//'Number of Groups',1)

   !this fills a nested level
   call set(dict//'First'//'One',1)
   call set(dict//'First'//'Two',2)

   !alternative way of filling
   dict2 => dict//'First'
   call set(dict//'First'//'Three',3)
   call set(dict2//'Threeb','3b')

   !print dictionary status
   call yaml_dict_dump(dict,flow=.true.)

   mat=1.0_f_double
   mat(2,3)=0.0_f_double

   !popping a term from the dictionary
   !only a normal pointer can be used
   !try with these examples
   call yaml_map('Size before removing',dict_size(dict2))
   call dict_remove(dict2,'One')
   call yaml_map('First removal',dict2)
   call dict_remove(dict2,'Two')
   call yaml_map('Second removal',dict2)
!   call pop(dict2,'Three')
   !a further element can be added
   call set(dict//'First'//'Four',4)
   call yaml_mapping_open('After pop')
   call yaml_dict_dump(dict)
   call yaml_mapping_close()

   !search for a key and point to it without modifying
   dict2=>find_key(dict,'Number of Gruops')
   call yaml_map('Key found',associated(dict2))
   !the key was wrong, try to find again
   dict2=>find_key(dict,'Number of Groups')
   call yaml_map('Second try, Key found',associated(dict2))
   ival=dict2
   call yaml_map('Value found',ival)
   !increase the value
   call set(dict//'Number of Groups',ival+1)
   !retrieve it
   ival=dict//'Number of Groups'
   call yaml_map('Alternative way',ival)

   call yaml_map('Search for "First" key',find_key(dict,'First'))

   !use now pop instead of remove
   !call dict_remove(dict,'First')
   !note that we do not have a garbage collector!
   ! imagine we use the dict_remove above
   !a call to this will produce a crash due to association above
   !indeed, dict2 points to dict//First, which has been removed
   !therefore dict2 is now pointing to a unallocated region
   !call set(dict2//'Five',5)

   !for the same reason, extracting the value in this way
   !ival = dict .pop. 'Number of Groups'
   !will lead to a memory leak, as the function in the right hand size
   !will not be freed

   !in fortran, the correct way to free the memory is the following
   dictA2 => dict .pop. 'Number of Groups'
   !extract value
   ival=dictA2
   call dict_free(dictA2)
   call yaml_map('Extracted value',ival)


   call yaml_map('Size after popping',dict_size(dict))
   dictA => dict .pop. 'First'
   call yaml_map('Size after popping again',dict_size(dict))
   call yaml_mapping_open('Complete pop')
   call yaml_map('Status of association',associated(dict))
   call yaml_map('Size after popping',dict_size(dict))
   call yaml_map('Present status',dict)
   call yaml_mapping_close()

   call yaml_map('DictA is associated',associated(dictA))
   call yaml_map('DictA is now of size',dict_size(dictA))
   call yaml_map('DictA is now',dictA)
   call yaml_map('DictA key',dict_key(dictA))
   call dict_free(dictA)
  !test if a complete pop will disassociate the dictionry
  call yaml_map('Dictionary associated before last pop',associated(dict))

   !now redo the pop experience with only one key in the dictionary
   call dict_init(dictA2)
   dict => dictA2 // 'head'
   call set(dict//'onlyone'//0,1.0)
   call set(dict//'onlyone'//1,1.0)
   call set(dict//'onlyone'//2,1.0)

   call yaml_map('Total dict',dictA2)

   if ('onlyone' .in. dict) then
      dictA => dict .pop. 'onlyone'
      call yaml_map('Nowdict',dictA)
      call dict_free(dictA)
   end if

   call yaml_map('Dict popped',dictA2)
   call dict_free(dictA2)
   
!  call dict_remove(dict,'Number of Groups')
!  call yaml_map('Last pop done, still associated',associated(dict))

   call dict_init(dictA)

   call dict_init(dictA2)

   call set(dictA2//'Test1'//'Toto',5)
   call set(dictA2//'Test1'//'Titi',6)

   call set(dictA//'Stack'//0,5)
   call set(dictA//'Stack'//1,4)
   call set(dictA//'Stack'//2,2)
   call set(dictA//'Stack'//3,dictA2)

   call set(dictA//'Stack2',(/'1','2','3'/))
   call set(dictA//'Stack3',(/'4 ','AQ','3g'/))
   call set(dictA//'Stack4',12)
   call set(dictA//'Matrix',mat)

   call yaml_dict_dump(dictA)

   !retrieve the value from the Stack2 key
   tmp_arr=dictA//'Stack2'

   call yaml_map('Values retrieved from the dict',tmp_arr,fmt='(1pg12.5)')

   dict2=>find_key(dictA,'Stack')
   call yaml_map('Lenght zero',dict_len(dict2))
   call yaml_map('Dict extracted',dict2)
   call dict_remove_last(dict2)
   call yaml_map('Lenght first',dict_len(dict2))


   call dict_remove_last(dict2)
   call yaml_map('Lenght second',dict_len(dict2))

   !  call push(dict2,'Element')
   !  call append(dictA,dictA2)
   call yaml_dict_dump(dictA)

   !retrieve the value from the Stack key
   tmp_arr(1:2)=dictA//'Stack'
   call yaml_map('Two values from Stack key',tmp_arr,fmt='(1pg12.5)')

   !retrieve the value from the a scalar
   tmp_arr=dictA//'Stack'//0
   call yaml_map('Array filled with a scalar',tmp_arr,fmt='(1pg12.5)')

   

   tmp_mat=45.0_f_double
   !retrieve now a rank-two array
   tmp_mat=dictA//'Matrix'

   !and verify its values
   call yaml_map('Small matrix',tmp_mat)
   call yaml_map('Original matrix',mat)

!!$   !try to see if extra information can be added after the value
!!$   call set(dictA//'Test Field',6,fmt='(i6.6)')
!!$   ival = dictA//'Test Field'
!!$   call yaml_map('Retrieving Test Field',ival)
!!$   call yaml_map('Retrieving actual value',dict_value(dictA//'Test Field'))
!!$
!!$
!!$   call set(dictA//'Test Field','6   #extra comment # extra')
!!$   ival = dictA//'Test Field'
!!$   call yaml_map('Retrieving Test Field Again',ival)
!!$   call yaml_map('Retrieving actual value',dict_value(dictA//'Test Field'))
!!$   call yaml_map('Index of comment',index(dict_value(dictA//'Test Field'),'#'))

   call yaml_comment('Prepend dictionary example',hfill='~')

   call yaml_map('Size of dict A',dict_size(dictA))
   call yaml_mapping_open('Dict A')
   call yaml_dict_dump(dictA)
   call yaml_mapping_close()


   call dict_init(dict2)
   call set(dict2//'Test1'//'Toto',5)
   call set(dict2//'Test1'//'Titi',6)
   call set(dict2//'Test2'//'Toto',4)
   call set(dict2//'Test2'//'Titi',2)


   call yaml_map('Size of dict 2',dict_size(dict2))
   call yaml_mapping_open('Dict 2')
   call yaml_dict_dump(dict2)
   call yaml_mapping_close()

   !verify the euqlity between dictionaries
   call yaml_map('Test1 and and Test2 are equal',dict2//'Test1' == dict2//'Test2')

   !now correct
   call set(dict2//'Test1'//'Toto',4)
   call set(dict2//'Test1'//'Titi',2)

   call yaml_map('Corrected version',dict2)

   !verify the equality between dictionaries
   call yaml_map('Test1 and and Test2 are equal',dict2//'Test1' == dict2//'Test2')

   !now add another element, written differently
   call set(dict2//'Test1'//'Tutu',4.d0,fmt='(1pe12.5)')
   call set(dict2//'Test2'//'Tutu','4.d0')

   call yaml_map('Added version',dict2)

   !verify the equality between dictionaries
   call yaml_map('Test1 and and Test2 are equal',dict2//'Test1' == dict2//'Test2')

   !also test the possibility two arrays filled with zeroes
   zero1=>list_new(.item. list_new(.item. '0.0000000000000000',&
        .item. '0.0000000000000000',.item. '0.0000000000000000'))
   zero2=>list_new(.item. list_new(.item. '0.',.item. '0.',.item. '0.'))

   call yaml_map('List of list of zeroes, first version',zero1)
   call yaml_map('List of list of zeroes, second version',zero2)

      !verify the equality between dictionaries
   call yaml_map('Zero1 and and Zero2 are equal',zero1==zero2)

   call dict_free(zero1)
   call dict_free(zero2)

   call yaml_map('Keys of first dict',dict_keys(dictA))
   call yaml_map('Keys of second dict',dict_keys(dict2))


   call prepend(dictA,dict2)
   call yaml_map('Size of prepended',dict_size(dictA))
   call yaml_mapping_open('Prepended')
   !call yaml_dict_dump2(dictA,verbatim=.true.)
   call yaml_dict_dump(dictA)
   call yaml_mapping_close()
   !test for scratching a dict
   !zero1=>dict_new('Test1' .is. 'scratchToto')
   !call set(dictA//'Test1','scratchToto')
   !test of the copy
   !zero1=>dict_new('Test1' .is. list_new(.item. 'scratchToto',.item. 'scratchTiti'))
   !call yaml_map('To update dict',zero1)
   !zero2 => dict_iter(zero1)
   !call yaml_map('Key found',zero2)
   !call dict_copy(dictA//'Test1',zero2)
   !zero2 => dict_next(zero2)
   !call yaml_map('Key found next',zero2)
   !!call dict_update(dictA,zero1)
   !call yaml_map('Scratched dict',dictA)
   call yaml_map('Keys of prepended dict',dict_keys(dictA))


   !perform an iterator on dictA
   dict_tmp=>dict_iter(dictA)
   do while(associated(dict_tmp))
      call yaml_map('Iterating in dictA',.true.)
      call yaml_map('Key of dictA',dict_key(dict_tmp))
      call yaml_map('Value of dictA',dict_value(dict_tmp))
      dict_tmp=>dict_next(dict_tmp)
   end do
   !end of perform an iterator on dictA

!!! [newiter]
   nullify(dict_tmp)
   do while(iterating(dict_tmp,on=dictA))
      call yaml_map('Iterating in dictA, again',.true.)
      call yaml_map('Key of dictA, again',dict_key(dict_tmp))
      call yaml_map('Value of dictA, again',dict_value(dict_tmp))
   end do
!!! [newiter]
   call dict_free(dictA)

   !test again the iterator in a compact form
   dictA=>dict_new('nolocks' .is. 'true')
   nullify(dict_tmp)
   do while(iterating(dict_tmp,on=dictA))
      call yaml_map('Iterating in dictA, compact form',.true.)
      call yaml_map('Key of dictA, again (CF)',dict_key(dict_tmp))
      call yaml_map('Value of dictA, again (CF)',dict_value(dict_tmp))
   end do
   call dict_free(dictA)

   !again the iterator on the elements of a list
   call yaml_mapping_open('Test of the iterator over a provided sequence')
   dictA=>yaml_load('[xyz,ascii,int,yaml]')
   call yaml_map('Provided list',dictA)
   nullify(dict_tmp)
   do while(iterating(dict_tmp,on=dictA))
      call yaml_map('Extension probed',trim(dict_value(dict_tmp)))
   end do
   call dict_free(dictA)
   call yaml_mapping_close()

   !fill a list and iterate over it
   dictA=>dict_new()
   do i=1,10
      call add(dictA,'Value'+yaml_toa(i))
   end do

   !perform an iterator on dict
   dict_tmp=>dict_iter(dictA)
   do while(associated(dict_tmp))
      call yaml_map('Item of dictA',dict_item(dict_tmp))
      call yaml_map('Value of dictA',dict_value(dict_tmp))
      dict_tmp=>dict_next(dict_tmp)
   end do
   call dict_free(dictA)

   !perform an iterator on a scalar, should not provide output
  dictA=>dict_new('Key' .is. 'Scalar')
   !end of perform an iterator
   dict_tmp=>dict_iter(dictA//'Key')
   do while(associated(dict_tmp))
      call yaml_map('Item of dict scalar',dict_item(dict_tmp))
      call yaml_map('Value of dict scalar',dict_value(dict_tmp))
      dict_tmp=>dict_next(dict_tmp)
   end do
   call dict_free(dictA)


   !example which has a bug
   dict_tmp => list_new((/.item.'55',.item. '66'/))
   dictA=>list_new((/.item. '5',.item. '6',.item. dict_tmp/))
!!$!call yaml_sequence_open("",flow=.false.)
!!$call yaml_sequence(advance="no")
!!$call yaml_mapping_open("SUCCESS",flow=.false.)
!!$call yaml_map("Id","0")
!!$call yaml_map("Message","Operation has succeeded")
!!$call yaml_map("Action","No action")
!!$call yaml_mapping_close()
!!$call yaml_sequence(advance="no")
!!$call yaml_mapping_open("GENERIC_ERROR",flow=.false.)
!!$call yaml_map("Id","1")
!!$call yaml_map("Message","UNSPECIFIED")
!!$call yaml_map("Action","UNKNOWN")
!!$call yaml_mapping_close()
!!$!call yaml_sequence_close()

   !what should be, also this writing has problem in the indentation
!!$    call yaml_sequence('5')
!!$    call yaml_sequence('6')
!!$    call yaml_sequence(advance='no')
!!$    call yaml_sequence_open()
!!$      call yaml_sequence('55')
!!$      call yaml_sequence('66')
!!$    call yaml_sequence_close()
!!$
   call yaml_sequence_open('List in a list')
   call yaml_dict_dump(dictA,verbatim=.true.)
   call yaml_dict_dump(dictA,flow=.false.)
   call yaml_dict_dump(dictA,flow=.true.,verbatim=.true.)
   call yaml_dict_dump(dictA,flow=.true.)
   call yaml_sequence_close()

   !perform an iterator with item on dict
   dict_tmp=>dict_next(dictA)
   do while(associated(dict_tmp))
      call yaml_map('Item of dictA',dict_item(dict_tmp))
      call yaml_map('Key of dictA',dict_key(dict_tmp))
      call yaml_map('Value of dictA',dict_value(dict_tmp))
      dict_tmp=>dict_next(dict_tmp)
   end do
   call dict_free(dictA)
   !end of perform an iterator with item
   

!!$   !try to steel a argument (does not work, should arrange routine set to be full-proof)
!!$   !fill a list and iterate over it
!!$   dictA=>dict_new()
!!$   do i=1,10
!!$      call add(dictA,trim(yaml_toa((/ (j,j=i,i+3) /))))
!!$   end do
!!$
!!$   call yaml_map('List before',dictA)
!!$
!!$   dict_tmp=>dict_new('ciao' .is. '1','hello' .is. '2')
!!$   dictA2=>dictA//3
!!$   call set(dict_tmp//'bonjour',dictA2)
!!$
!!$   call yaml_map('Thief dict',dict_tmp)
!!$
!!$   call yaml_map('List after',dictA)
!!$   call dict_free(dictA)
!!$   call dict_free(dict_tmp)

   name(1:len(name)) = 'atom'
   call dict_init(dict_tmp)
   dictA => dict_tmp // 'Test'
   call dict_init(dictA2)
   call set(dictA2 // name // 0, tmp_arr(1))
   call set(dictA2 // name // 1, tmp_arr(2))
   call set(dictA2 // name // 2, tmp_arr(3))
   call add(dictA, dictA2)

            !to be verified if one of these works
!!$            call set(fxyz // astruct%atomnames(astruct%iatype(iat)),outs%fxyz(:,iat))
!!$            call set(fxyz // astruct%atomnames(astruct%iatype(iat)), list_new( .elem. outs%fxyz(:,iat)))
!!$            fxyz => dict_new( trim(astruct%atomnames(astruct%iatype(iat))) .is. outs%fxyz(:,iat))
!!$
!!$            call f_strcpy(src=astruct%atomnames(astruct%iatype(iat)),dest=name)
!!$            call add(pos, dict_new( trim(name) .is. &
!!$                 list_new( .elem. outs%fxyz(:,iat))))
   call yaml_map('Dictionary now',dictA)
   call dict_free(dict_tmp)

   !now do it again with another mechanism
   tmp_arr = tmp_arr + 1.d0
   call dict_init(dict_tmp)
   dictA => dict_tmp // 'Test'
   call dict_init(dictA2)
   call set(dictA2 // name,[ '5', '6', '7', '8' ])
   call set(dictA2 // name,[ '3', '4' ])
   call set(dictA2 // name,[ '2', '9', '1' ])
   call set(dictA2 // name,tmp_arr)
   call add(dictA, dictA2)
   call yaml_map('Dictionary now2',dictA)
   call dict_free(dict_tmp)

   !now do it again2 with another mechanism
   tmp_arr = tmp_arr + 1.d0
   call dict_init(dict_tmp)
   dictA => dict_tmp // 'Test'
   call add(dictA, dict_new(name .is. list_new( .item. tmp_arr)))
   call yaml_map('Dictionary now3',dictA)
   call dict_free(dict_tmp)

   !now do it again3 with another mechanism (the one of choice)
   tmp_arr = tmp_arr + 1.d0
   call dict_init(dict_tmp)
   dictA => dict_tmp // 'Test'
   call add(dictA, dict_new(name .is. tmp_arr))
   call yaml_map('Dictionary now3',dictA)

   !now iterate on the list in two ways
   dictA2=> dictA // 0 // name
   call yaml_map('Dictionary',dictA2)
   do i=1,dict_len(dictA2)
      call yaml_map('i='+yaml_toa(i),dictA2//(i-1))
   end do

   dictA2=> dict_iter(dictA // 0 .get. name)
   i=0
   do while(associated(dictA2))
      i=i+1
      call yaml_map('i2='+yaml_toa(i),dictA2)
      dictA2 => dict_next(dictA2)
   end do

   call dict_free(dict_tmp)

 end subroutine test_dictionaries1

 subroutine test_copy_merge()
   use dictionaries
   use yaml_output
   use yaml_parse
   use yaml_strings
   implicit none

   type(dictionary), pointer :: dict, cpy, subd

   dict => dict_new(&
         & "__comment__" .is. 'Grid shifts', &
         & "__cond__"    .is. dict_new("__master_key__" .is. "kpt_method", "__when__" .is. list_new( .item. "MPGrid")), &
         & "__default__" .is. list_new( .item."0.", .item."0.", .item."0.") )

   call yaml_mapping_open("test dict_copy")
   call yaml_mapping_open("original")
   call yaml_dict_dump(dict)
   call yaml_mapping_close()
   nullify(cpy)
   call dict_copy(cpy, dict)
   call yaml_mapping_open("copy")
   call yaml_dict_dump(cpy)
   call yaml_mapping_close()
   call dict_free(cpy)
   call yaml_mapping_close()

   cpy => yaml_load('{__comment__ : Grid shifts, __cond__: '//&
        '   { __master_key__: kpt_method, __when__ : [ MPGrid ]},'//&
         '__default__ : ['//0.//','//0.//','// 0.//'] }')
   call yaml_mapping_open("new method")
   call yaml_dict_dump(cpy)
   call yaml_map('dicts are equal',cpy==dict)
   call dict_free(cpy)
   call yaml_mapping_close()

   !another comprehensive test
   subd => dict_new(  &
         & "__exclusive__" .is. dict_new( "123" .is. "operation 123", &
         &                                  "456" .is. "operation 456" ), &
         & "__default__"   .is. list_new(.item."1.", .item."2.", .item."3." ) )
   !end of another comprehensive test
   call yaml_mapping_open("test dict_update")
   call dict_update(dict, subd)
   call yaml_mapping_open("additional")
   call yaml_dict_dump(subd)
   call yaml_mapping_close()
   call yaml_mapping_open("after merge")
   call yaml_dict_dump(dict)
   call yaml_mapping_close()
   call yaml_mapping_close()
   call dict_free(subd)

   call dict_free(dict)
 end subroutine test_copy_merge

subroutine test_dictionary_for_atoms()
  use yaml_output
  use yaml_strings
  implicit none

!!$  character(len = 50) :: gu
  integer :: ierr
  character(len = 50) :: fmts,tmp
  double precision, dimension(3) :: cell, xred, hgrids
  double precision :: tt


  call yaml_mapping_open("Atomic structure")

  cell = 20.345752999999998
  call yaml_map('Cell', cell)

  hgrids = cell / (/ 54, 40, 40 /)

  call yaml_sequence_open('Positions')

  call yaml_sequence(advance='no')
  xred = (/ 0.2516085125D-05,  0.5826606155D-05,  20.34574212d0 /)
  call print_one_atom('Si',xred,hgrids,1)

  call yaml_sequence(advance='no')
  xred = (/ 5.094032326d0,  5.153107111d0,  0.3047989908d-01 /)
  call print_one_atom('Si',xred,hgrids,2)
!!$  call yaml_map("Si", xred, fmt="(g18.10)", advance = "no")
!!$  xred = xred / hgrids
!!$  write(gu, "('[ 'F6.2', 'F6.2', 'F6.2'] 'I4.4)") xred, 2
!!$  call yaml_comment(gu)

  call yaml_sequence(advance='no')
  xred = (/ 0.3049344014d-01,  5.153107972d0,  5.094018600d0 /)
  call print_one_atom('Si',xred,hgrids,3)
!!$  call yaml_map("Si", xred, fmt="(g18.10)", advance = "no")
!!$  xred = xred / hgrids
!!$  write(gu, "('[ 'F6.2', 'F6.2', 'F6.2'] 'I4.4)") xred, 3
!!$  call yaml_comment(gu)

  call yaml_sequence_close()

  call yaml_mapping_close()

  !now print some double precision values to understand which is the best format
  tt=real(0.5e0,kind=8) !use a conversion from float
  call yaml_map('Test clean_zeroes',clean_zeroes(yaml_toa('20')))
  call yaml_map('Retest clean_zeroes',clean_zeroes(yaml_toa('20.0e+00')))
  call yaml_map('Real without format',clean_zeroes(yaml_toa('0.2000000000000000000')))
  fmts(1:len(fmts))='(1pe25.17)'
  call yaml_map('Real with format '//trim(fmts),clean_zeroes(yaml_toa(tt,fmt=fmts)))
  fmts(1:len(fmts))='(1pe24.16)'
  call yaml_map('Real with format '//trim(fmts),tt,fmt=fmts)
  fmts(1:len(fmts))='(1pe24.16)'
  call yaml_map('Tiny with wrong format '//trim(fmts),tiny(tt),fmt=fmts)
  write(tmp,fmt=fmts,iostat=ierr)tiny(tt)
  call yaml_map('Iostat for format',ierr)
  fmts(1:len(fmts))='(1pe24.16)'
  call yaml_map('Huge with wrong format '//trim(fmts),-huge(tt),fmt=fmts)
  write(tmp,fmt=fmts,iostat=ierr)-huge(tt)
  call yaml_map('Iostat for format',ierr)
  fmts(1:len(fmts))='(i3)'
  call yaml_map('Integer with too little format '//trim(fmts),10000,fmt=fmts)
  write(tmp,fmt=fmts,iostat=ierr) 10000
  call yaml_map('Iostat for format',ierr)
  fmts(1:len(fmts))='(f3.2)'
  call yaml_map('Float with too little format '//trim(fmts),1000.4,fmt=fmts)
  fmts(1:len(fmts))='(f3.0)'
  call yaml_map('Real as integer '//trim(fmts),10.0,fmt=fmts)
  write(tmp,fmt=fmts,iostat=ierr) 1000.4
  call yaml_map('Iostat for format',ierr)
  fmts(1:len(fmts))='(f3.2)'
  call yaml_map('Double with too little format '//trim(fmts),1000.4d0,fmt=fmts)
  fmts(1:len(fmts))='(f3.2)'
  call yaml_map('Float with too little format again '//trim(fmts),10.4456,fmt=fmts)
  write(tmp,fmt=fmts,iostat=ierr)10.4456
  call yaml_map('Iostat for format',ierr)
  fmts(1:len(fmts))='(es23.16)'
  call yaml_map('Real with format '//trim(fmts),tt,fmt=fmts)
  fmts(1:len(fmts))='(es24.17)'
  call yaml_map('Real with format '//trim(fmts),tt,fmt=fmts)
  fmts(1:len(fmts))='(es25.18)'
  call yaml_map('Real with format '//trim(fmts),tt,fmt=fmts)
  fmts(1:len(fmts))='(es26.19)'
  call yaml_map('Real with format '//trim(fmts),tt+epsilon(1.d0),fmt=fmts)
  fmts(1:len(fmts))='(es27.20)'
  call yaml_map('Real with format '//trim(fmts),tt-epsilon(1.d0),fmt=fmts)
  fmts(1:len(fmts))='(es26.19)'
  call yaml_map('Real with format '//trim(fmts),epsilon(1.d0),fmt=fmts)
  fmts(1:len(fmts))='(es27.20)'
  call yaml_map('Real with format '//trim(fmts),-epsilon(1.d0),fmt=fmts)

  !Cite a missing citation
  call yaml_cite('C61')
  !call a valid paper
  call yaml_cite('C60')

  contains

    subroutine print_one_atom(atomname,rxyz,hgrids,id)
      implicit none
      integer, intent(in) :: id
      character(len=*), intent(in) :: atomname
      double precision, dimension(3), intent(in) :: rxyz,hgrids
      !local variables
      character(len=*), parameter :: fmtat='(g18.10)',fmtg='(F6.2)'
      integer :: i

      call yaml_sequence_open(atomname,flow=.true.)
      do i=1,3
         call yaml_sequence(yaml_toa(rxyz(i),fmt=fmtat))
      end do
      call yaml_sequence_close(advance='no')
      call yaml_comment(trim(yaml_toa(rxyz/hgrids,fmt=fmtg))//trim(yaml_toa(id))) !we can also put tabbing=

    end subroutine print_one_atom

    !> when str represents a real number, clean it if there are lot of zeroes after the decimal point
    !! pure
    function clean_zeroes(str)
      implicit none
      integer, parameter:: max_value_length=95
      character(len=*), intent(in) :: str
      character(len=max_value_length) :: clean_zeroes
      !local variables
      integer :: idot,iexpo,i

      !first fill with all the values up to the dot if it exist
      idot=scan(str,'.')
      if (idot==0) then
         !no dot, nothing to clean
         clean_zeroes(1:max_value_length)=str
      else
         !then search for the position of the exponent or of the space if present
         iexpo=scan(str(idot+2:),'eE ')+idot+1
         !print *,'there',iexpo,'str',str(idot+2:)
         if (iexpo==idot+1) iexpo=len(str)+1
         i=iexpo
         find_last_zero: do while(i > idot+1) !first digit after comma always stays
            i=i-1
            if (str(i:i) /= '0') exit find_last_zero
         end do find_last_zero
         clean_zeroes(1:i)=str(1:i)
         !print *,'here',i,clean_zeroes(1:i),'iexpo',iexpo,str(iexpo:)
         !then copy the exponent
         if (str(iexpo:) /= 'E+00' .and. str(iexpo:) /= 'e+00' .and. str(iexpo:) /= 'E+000' .and. &
              str(iexpo:) /= 'e+000') then
            call f_strcpy(src=str(iexpo:),dest=clean_zeroes(i+1:))
            !clean_zeroes(i+1:max_value_length)=str(iexpo:)
         else
            clean_zeroes(i+1:max_value_length)=' '
         end if
      end if
    end function clean_zeroes

end subroutine test_dictionary_for_atoms


!> Test the usage of the new f_trees structure
subroutine test_f_trees()
  use f_trees
  use yaml_output
  implicit none
  type(f_tree) :: dict1

  !initialization
  dict1=f_tree_new()

  call f_tree_push(dict1//'Key1','val1')
  !to be updated with proper method for trees
  call yaml_map('Initialized tree',dict1%d)

  call yaml_mapping_open('Test of dump')
  call f_tree_dump(dict1)
  call yaml_mapping_close()

  call f_tree_free(dict1)
end subroutine test_f_trees


!> This routine consider the usage of dictionaries for intensive data storage (of course to be avoided)
!! and compares it to the usage of an array for doing similar things
subroutine profile_dictionary_usage()
  use dictionaries
  use f_utils
  use yaml_output
  implicit none
  !local variables
  integer :: nprof,ntry,nstep,iprof,jprof,itry,ival
  integer(kind=8) :: t0,t1
  double precision :: tel,tot
  type(dictionary), pointer :: dict
  integer, dimension(:), allocatable :: itest !< used to simulate search with an array
  type(f_progress_bar) :: bar


!!$!$  !profiling
  nprof=100001
  ntry=1000
  nstep=10000
  allocate(itest(nprof))
  itest=0
  bar=f_progress_bar_new(nstep=nprof/nstep)
  do iprof=1,nprof,nstep

     !call system_clock(ncount0,ncount_rate,ncount_max)
     t0=f_time()
     do itry=1,ntry
        do jprof=1,nprof
           itest(jprof)=itest(jprof)+iprof+itry
        end do
     end do
     t1=f_time()
     !call system_clock(ncount1,ncount_rate,ncount_max)
     !tel=dble(ncount1-ncount0)/dble(ncount_rate)*(1d6/dble(ntry))
     tot=dble(ntry)*dble(nprof)
     tel = dble(t1-t0)/tot
     !here all the steps are identical
     call dump_progress_bar(bar,step=iprof/nstep)
!!$     call yaml_mapping_open('Timings for search',flow=.true.)
!!$     call yaml_map('No. of items',iprof)
!!$     call yaml_map('Elapsed time (ns)',tel,fmt='(f12.2)')
!!$     call yaml_mapping_close()
  end do
  call yaml_map('Some value',itest(1)+itest(ntry))
  deallocate(itest)

  !profiling
  nprof=20001
  ntry=100
  nstep=5000
  call dict_init(dict)
  bar=f_progress_bar_new(nstep=nprof/nstep)
  do iprof=1,nprof,nstep
     do jprof=0,nstep-1
        call set(dict//'Test'//(jprof+iprof-1),jprof+iprof-1)
     end do

     !call system_clock(ncount0,ncount_rate,ncount_max)
     tot=0.d0
     t0=f_time()
     do itry=1,ntry
        ival=dict//'Test'//(iprof-1)
        tot=tot+dble(ival)
     end do
     t1=f_time()
     !call system_clock(ncount1,ncount_rate,ncount_max)
     !tel=dble(ncount1-ncount0)/dble(ncount_rate)*(1d6/dble(ntry))
     tel=dble(t1-t0)/dble(ntry)*1.d-3
     !here each step costs more
     call dump_progress_bar(bar,step=iprof/nstep)

!!$     call yaml_mapping_open('Timings for search',flow=.true.)
!!$     call yaml_map('No. of items',iprof)
!!$     call yaml_map('Elapsed time (mus)',tel,fmt='(f12.2)')
!!$     call yaml_mapping_close()
  end do
  call yaml_map('Other value',tot)
  call dict_free(dict)


end subroutine profile_dictionary_usage
