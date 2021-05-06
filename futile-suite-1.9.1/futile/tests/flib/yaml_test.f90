!> @file
!! Test yaml_output module
!! @example yaml_test.f90
!! Extensive tests about yaml output generations
!! @author
!!    Copyright (C) 2012-2015 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Test yaml output module
program yaml_test
   use yaml_output
   use dictionaries, dict_char_len=> max_field_length
   use dynamic_memory
   use yaml_parse
   use f_utils
   use f_bibliography
   implicit none
   integer, parameter :: isize=16 !< size of functionality
   character(len=isize), parameter :: YAML               ='yaml'
   character(len=isize), parameter :: YAML_EXTRAS        ='yaml_extras'
   character(len=isize), parameter :: YAML_PARSER        ='yaml_parse'
   character(len=isize), parameter :: EXCEPTIONS         ='exceptions'
   character(len=isize), parameter :: TREES              ='trees'
   character(len=isize), parameter :: TREES_EXTRAS       ='trees_extras'
   character(len=isize), parameter :: ALLOCATIONS        ='allocations'
   character(len=isize), parameter :: UTILS              ='utils'
   character(len=isize), dimension(8), parameter :: FUNCTIONALITIES=&
        [ YAML               ,&
          YAML_EXTRAS        ,&
          YAML_PARSER        ,&
          EXCEPTIONS         ,&
          TREES              ,&
          TREES_EXTRAS       ,&
          UTILS              ,&
          ALLOCATIONS        ]
   
   type(dictionary), pointer :: dict_tmp,run,dict_mp,biblio
   type(yaml_cl_parse) :: parser
   external :: get_database

   call f_lib_initialize()
!   call yaml_set_stream(record_length=92)
   !test output level
   call f_malloc_set_status(output_level=2,logfile_name='memstatus.yaml')

   parser=yaml_cl_parse_null()
   !set valid options
   call yaml_cl_parse_option(parser,'test','All',&
        'Decide what to test (--help for more info)','t',&
        dict_new('Usage' .is. &
        'Specify, as a yaml list, the functionalities which have to be tested',&
        'Allowed values' .is. &
        list_new(.item. FUNCTIONALITIES)))
   !other options to test the parser
   call yaml_cl_parse_option(parser,'test2','1',&
        'this is another test option','s',first_option=.true.)
   call yaml_cl_parse_option(parser,'test3','None',&
        'this is a test option','T',&
        dict_new('Usage' .is. &
        'Unused option3, just for testing the command line parser, '//&
        'also the long help lines have to be tested in order to understand if it works'),&
        conflicts='[test2,test]')
   call yaml_cl_parse_option(parser,'test_mp','None',&
        'sandbox to test the I/O of point multipoles','m',&
        dict_new('Usage' .is. &
        'Just to test the format of the multipoles'))

   !verify the parsing
   call yaml_cl_parse_cmd_line(parser)

   call yaml_map('Parsed options',parser%options)
   call yaml_map('Parsed info',parser%args)

   !construct the runs
   nullify(run)
   dict_tmp = parser%args .get. 'test'
   if (trim(dict_value(dict_tmp)) == 'All') then
      run => list_new(.item. FUNCTIONALITIES)
   else if (dict_len(dict_tmp) > 0) then
      call dict_copy(run,dict_tmp)
   end if
   nullify(dict_tmp)

   dict_mp = parser%args .get. 'test_mp'
   call yaml_map('Multipole list found',associated(dict_mp))
   if (associated(dict_mp)) then
       call check_multipoles(parser)
   end if

!    call dict_init(dict_tmp)
!   !call set(dict_tmp//'Ciao','1')
!    call set(dict_tmp//0,'1')
!    call set(dict_tmp//1,'5')
!    dict_tmp=>dict_new() !d={}
!    dict_tmp=>dict_new('Ciao' .is. '1')
!    dict_tmp=>list_new() !d=[]
!    dict_tmp=>list_new([.item. '1',.item. '5']) !d=[1,5]
!    dict_tmp=>list_new(.item. ['1','ciao']) !d=[1,ciao]

   call yaml_cl_parse_free(parser)
!!$
!!$
!!$   !call profile_dictionary_usage()

   !take database information
   call yaml_parse_database(biblio,get_database)
   call f_bib_update(biblio//0) !only first document
   call dict_free(biblio)

   if (YAML .in. run) then
      !test to solve the bug of yaml_comment
      call yaml_comment('AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA')
      call yaml_comment('BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB')
      call yaml_comment('CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC')
      call yaml_comment('DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD')

      !First document  
      call yaml_new_document()
      call yaml_comment('Yaml Output Module Test',hfill='~')
      call test_yaml_output1()
      call yaml_release_document()
      
      !Second document
      call yaml_new_document()
      call test_yaml_output2()
      call yaml_release_document()
   end if

   if (YAML_EXTRAS .in. run) then
      !Third document
      call yaml_new_document()
      !Check calling twice yaml_new_document
      call yaml_new_document()
      call test_yaml_output_sequences1()
      call yaml_release_document()

   !Fourth document
      call yaml_new_document()
      call test_yaml_output_sequences2()
      call yaml_release_document()

   !Fourth-A document
      call yaml_new_document()
      call yaml_invoice_example()
      call yaml_release_document()

   !Fourth-A2 document
      call yaml_new_document()
      call yaml_invoice_example_with_dictionaries()
      call yaml_release_document()
   end if

   if (EXCEPTIONS .in. run) then
      !Fourth-B document, to be moved
      call yaml_new_document()
      call test_error_handling()
      call yaml_release_document()
   end if

   if (TREES .in. run) then
      !Fifth document
      call yaml_new_document()
      call test_dictionaries0()
      call yaml_release_document()
      
      !Sixth document: test dictionaries
      call yaml_new_document()
      call test_dictionaries1()
      call yaml_release_document()
   end if

   if (ALLOCATIONS .in. run) then
      !Seventh document: Test dynamic memory allocation
      call yaml_new_document()
      call test_dynamic_memory()
      call yaml_release_document()
      dict_tmp=>dict_new()
      call f_malloc_dump_status(dict_summary=dict_tmp)
      call yaml_map('Summary',dict_tmp)
      call dict_free(dict_tmp)
      call yaml_map('Test for pid',f_getpid())
   end if

   if (TREES .in. run) then
      call yaml_new_document()
      call test_copy_merge()
      call yaml_release_document()


      call yaml_new_document()
      call test_f_trees()
      call yaml_release_document()

   end if

   if (YAML_PARSER .in. run) then
      !test the yaml parsing
      call yaml_parse_file_and_string()
   end if
   
   if (YAML_EXTRAS .in. run) then
      call yaml_new_document()
      call test_dictionary_for_atoms()
      call yaml_release_document()
   end if
   if (TREES_EXTRAS .in. run) then
      call profile_dictionary_usage()
   end if
!!$   if (ALLOCATIONS .in. run) then
!!$      call verify_heap_allocation_status()
!!$   end if
   if (UTILS .in. run) then
      call f_utils_test()
      call f_inputfile_test()
   end if

   call dict_free(run)

   !prepare the finalization of the library
   call f_lib_finalize()

  call f_lib_initialize()
  !Finalize without report
  call f_lib_finalize_noreport()


end program yaml_test

subroutine yaml_parse_file_and_string()
  use dictionaries
  use yaml_parse
  use yaml_output
  use f_precisions, only: cr=>f_cr
  implicit none
  type(dictionary), pointer :: dict_parse
  character(len=*), parameter ::stream=&
       "---"//cr//&
       "Key1: field1"//cr//&
       "Key2: "//cr//&
       " Subkey2-1: ciao"//cr//&
       " Subkey2-2: [0, hello, 2]"//cr//&
       "Key3:"//cr//&
       " - One"//cr//&
       " - Two"//cr//&
       " - Three"//cr//&
       "---"//cr//&
       "#example of the input variables"//cr//&
       "inputvars: #name of the variable as declared in the code"//cr//&
       " COMMENT: 'This is the description of the variable as will appear in the logfile'"//cr//&
       " RANGE: ['from','to'] #always put two numbers (also .inf) can be put"//cr//&
       " EXCLUSIVE: #here follows a list of allowed values (either RANGE or EXCLUSIVE)"//cr//&
       "  - Value1:  'comments of value1'"//cr//&
       "  - Value2:  'comment of value2'"//cr//&
       " CONDITION: #here the conditions for which the variable makes sense are written "//cr//&
       "   MASTER_KEY: foo #this means that inputvars makes sense only if foo is specified "//cr//&
       "   WHEN: #provide a list of allowed values of foo for which inputvars is meaningful "//cr//&
       "     - fooval1 "//cr//&
       "     - fooval2   "//cr//&
       "#then the profiles follows, which gives to the variables the allowed name "//cr//&
       " default: 'value of the default, written as a string' "//cr//&
       " profile1: 'value1' # if the user specifies inputvars: profile1 then inputvars will be value1 "//cr//&
       " profile2: 'value2'"//cr

  
  call yaml_comment('Yaml parsing',hfill='-')

!  call yaml_parse_from_file(dict_parse,'testfile.yaml')
!  call yaml_dict_dump_all(dict_parse)
!  call dict_free(dict_parse)

  call yaml_parse_from_string(dict_parse,stream)
  call yaml_dict_dump_all(dict_parse)
  call dict_free(dict_parse)

  
end subroutine yaml_parse_file_and_string

  subroutine help_screen()
    write(*,*)' Usage of the command line instruction'
    write(*,*)' --taskgroup-size=<mpi_groupsize>'
    write(*,*)' --runs-file=<list_posinp filename>'
    write(*,*)' --run-id=<name of the run>: it can be also specified as unique argument'
    write(*,*)' --help : prints this help screen'
  end subroutine help_screen



subroutine check_multipoles(parser)
  use yaml_output
  use yaml_strings
   use dictionaries
   !use dynamic_memory
   use yaml_parse
   !use f_utils
  implicit none
  ! Calling arguments
  type(yaml_cl_parse) :: parser

  ! Local variables
  type(dictionary), pointer :: dict_mp, iter
  integer :: ilist, imp, nmplist
  character(len=2) :: key
  real(kind=8),dimension(3) :: rxyz
  real(kind=8),dimension(7) :: mp
  character(len=20) :: sym

   !!!before freeing the options just test a putative way of inserting multipoles
   !!!first retrieve the dictionary if it has been entered
   dict_mp = parser%args .get. 'test_mp'
   !!call yaml_map('Multipole list found',associated(dict_mp))

   if  (associated(dict_mp)) then
       nmplist = dict_len(dict_mp)
       call yaml_map('Size of the mp list',nmplist)
       call yaml_sequence_open('Values')
       do ilist=0,nmplist-1
          call yaml_sequence()
          call yaml_map('Size of element'//trim(yaml_toa(ilist)),dict_size(dict_mp//ilist))
          iter => dict_mp//ilist
          ! Retrieve the description of the multipole center (typically atom type)
          sym = iter//'name'
          !retrieve atomic positions, compulsory
          if ('r' .notin. iter) call f_err_throw('For the item .. the r should  be present')
          rxyz = iter//'r'
          call yaml_map('rxyz',rxyz)
          do imp=0,3
               key='q'+imp
               if (key .in. iter) then
                 if (dict_len(iter//key)/=2*imp+1) then
                     call f_err_throw('Wrong len ('//dict_len(iter//key)//') of the mutipole')
                 end if
                 mp(1:2*imp+1) = iter//key
                 !call yaml_map('dict_size of '//key,dict_size(iter//key))  
                 !call yaml_map('dict_len of '//key,dict_len(iter//key))  
                 call yaml_map(key,mp(1:2*imp+1))
                 !call yaml_map('key',key)
                 !call yaml_map('key',iter//key)
                 !call yaml_dict_dump(iter//key)  
               end if
          end do 
       end do
       call yaml_sequence_close()
   end if

end subroutine check_multipoles
