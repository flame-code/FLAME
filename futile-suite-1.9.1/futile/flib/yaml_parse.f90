!> @file
!! Module to parse the yaml (flib library)
!! @author
!!    Copyright (C) 2013-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Module containing routines to parse a yaml input
module yaml_parse
  use dictionaries, only: dictionary,max_field_length
  implicit none

  private

  integer :: STREAM_START, STREAM_END
  integer :: DOCUMENT_START, DOCUMENT_END
  integer :: SEQUENCE_START, SEQUENCE_END
  integer :: MAPPING_START, MAPPING_END
  integer :: ALIAS, SCALAR, ERROR

  !> Contains some keys to better explain the yaml parse errors
  type(dictionary), pointer :: dict_yaml_errs=>null()

  !> Two possible errors when parsing
  integer, public :: YAML_PARSE_ERROR       = 0
  integer, public :: YAML_PARSE_UNSUPPORTED = 0
  integer, public :: ERROR_YAML_COMMAND_LINE_PARSER = 0

  character(len=*), parameter :: OPTDEFAULT = 'default'
  character(len=*), parameter :: OPTSNAME   = 'shortname'
  character(len=*), parameter :: OPTSHELP   = 'shorthelp'
  character(len=*), parameter :: OPTHELP    = 'help'
  character(len=*), parameter :: OPTCONFL    = 'conflicts'

  type, public :: yaml_cl_parse
     !> command line parser to determine options
     !>default value of the key for the first command line options
     character(len=max_field_length) :: first_command_key
     !> dictionary specifying the valid options
     type(dictionary), pointer :: options
     !> parsed dictionary, filled by options which have valid default value
     type(dictionary), pointer :: args
  end type yaml_cl_parse

  interface yaml_cl_parse_option
     module procedure yaml_cl_parse_option_from_string,yaml_cl_parse_option,yaml_cl_parse_option_from_dict
  end interface

  public :: yaml_parse_from_file
  public :: yaml_parse_from_char_array
  public :: yaml_parse_from_string,yaml_parse_database
  public :: yaml_cl_parse_null,yaml_cl_parse_free
  public :: yaml_cl_parse_option,yaml_argparse

  !for internal f_lib usage
  public :: yaml_parse_errors
  public :: yaml_parse_errors_finalize
  public :: yaml_load,yaml_cl_parse_cmd_line

contains


  !creators
  pure function yaml_cl_parse_null() result(parser)
    implicit none
    type(yaml_cl_parse) :: parser
    call nullify_yaml_cl_parse(parser)
  end function yaml_cl_parse_null
  pure subroutine nullify_yaml_cl_parse(parser)
    use yaml_strings, only: f_strcpy
    implicit none
    type(yaml_cl_parse), intent(out) :: parser

    call f_strcpy(src=' ',dest=parser%first_command_key)
    nullify(parser%options)
    nullify(parser%args)
  end subroutine nullify_yaml_cl_parse

  subroutine yaml_cl_parse_free(parser)
    use dictionaries_base, only: dict_free
    implicit none
    type(yaml_cl_parse), intent(inout) :: parser
    call dict_free(parser%options)
    call dict_free(parser%args)
    parser = yaml_cl_parse_null()
  end subroutine yaml_cl_parse_free

  !> used to set up input parameters
  !! similar behaviour as python optparse options
  subroutine yaml_cl_parse_option(parser,name,default,help_string,&
       shortname,help_dict,first_option,conflicts)
    use dictionaries
    use yaml_strings
    implicit none
    !> the parser which has to be updated
    type(yaml_cl_parse), intent(inout) :: parser
    !> name of the parser option, should correspond to the name of the long key invoked with --name=value
    !! the only protected name is "help" as the latter dumps out the help screen
    character(len=*), intent(in) :: name
    !> value of the default option. If "None" is given, after parsing the parser will not contain this key
    character(len=*), intent(in) :: default
    !> help string of the key, should be interesting for the user as a short help
    character(len=*), intent(in) :: help_string
    !> if true, it states that this is the first option of the parser, which can be given without key
    !! if given in multiple options, the last given overrides the others
    logical, intent(in), optional :: first_option
    !> shortkey name. in this case the option invoked is as -k value
    !! the only protected key if h as it stands for shorthelp
    character(len=1), intent(in), optional :: shortname
    !> character string containing list (in yaml format) of
    !! options that conflicts with the one given
    !! the options may not have been defined yet.
    !!The check is performed at the end of the parsing procedure
    character(len=*), intent(in), optional :: conflicts

    !> help dict intended as a more eextended help, to be invoked when help on command line is given
    !! this dictionary is stolen by the parser and nullified at exit
    type(dictionary), pointer, optional :: help_dict
    !local variables
    logical :: first
    character(len=1) :: shname
    type(dictionary), pointer :: option,iter

    if (trim(name) == 'help') then
       call f_err_throw('The "help" option is reserved.',err_id=ERROR_YAML_COMMAND_LINE_PARSER)
       return
    end if

    if (trim(name) .in. parser%options) then
       call f_err_throw('Option "'//trim(name)//&
            '" already specified',err_id=ERROR_YAML_COMMAND_LINE_PARSER)
       return
    end if

    option => dict_new(OPTDEFAULT .is. default, OPTSHELP .is. help_string)

    if (present(shortname)) then
       if (shortname == 'h') then
          call f_err_throw('Option shortname "h" is reserved for shorthelp.',&
               err_id=ERROR_YAML_COMMAND_LINE_PARSER)
          return
       end if
       if (len_trim(shortname) /=0) then
          !search if another shortname has already been specified
          shname='h' !protected value
          iter => dict_iter(parser%options)
          do while (associated(iter))
             shname = iter .get. OPTSNAME
             if (shname == shortname) then
                call f_err_throw('Option "'//trim(name)//&
                     '" error; shortname "'//trim(shortname)//&
                     '" already specified in option "'//&
                     trim(dict_key(iter))//'"',&
                     err_id=ERROR_YAML_COMMAND_LINE_PARSER)
                return
             end if
             iter => dict_next(iter)
          end do

          call set(option//OPTSNAME,shortname)
       end if
    end if
    if (present(help_dict)) then
       if (associated(help_dict)) call set(option//OPTHELP,help_dict)
    end if

    if (present(conflicts)) then
       if (trim(default) /= "None" .and. len_trim(conflicts) /=0) then
          call f_err_throw('Error in defining "'//&
               trim(name)//&
               '"; a conflicting option must have "None" as its default value.',&
               err_id=ERROR_YAML_COMMAND_LINE_PARSER)
          return
       end if

       if (len_trim(conflicts) /=0) then
          iter=>yaml_load(trim(conflicts),OPTCONFL)
          call dict_update(option,iter)
          call dict_free(iter)
       end if
    end if

    !then insert it in the dict of options
    if (.not. associated(parser%options)) call dict_init(parser%options)
    !then add the option
    call set(parser%options//trim(name),option)

    first=.false.
    if (present(first_option)) first=first_option
    if (first) then
       !check if this has already been initialized
       if (len_trim(parser%first_command_key) /= 0) then
          call f_err_throw('Option "'//trim(parser%first_command_key)//&
               '" has already been set as first_option.',&
               err_id=ERROR_YAML_COMMAND_LINE_PARSER)
          return
       end if
       call f_strcpy(src=name,dest=parser%first_command_key)
    end if

  end subroutine yaml_cl_parse_option

  !> accept the definition of the input variable from a yaml_string
  subroutine yaml_cl_parse_option_from_string(parser,string)
    use dictionaries
    implicit none
    !> the parser which has to be updated
    type(yaml_cl_parse), intent(inout) :: parser
    !>definition of the input variable, contains the informations
    !in a dictionary
    character(len=*), intent(in) :: string
    !local variables
    type(dictionary), pointer :: dict,iter

    !load the input variable definition
    dict=>yaml_load(string)

    if (dict_islist(dict)) then
       nullify(iter)
       do while(iterating(iter,on=dict))
          call yaml_cl_parse_option_from_dict(parser,iter)
       end do
    else
       call yaml_cl_parse_option_from_dict(parser,dict)
    end if
    call dict_free(dict)

  end subroutine yaml_cl_parse_option_from_string

  subroutine yaml_cl_parse_option_from_dict(parser,dict)
    use dictionaries
    use f_utils, only: f_zero
    implicit none
    !> the parser which has to be updated
    type(yaml_cl_parse), intent(inout) :: parser
    type(dictionary), pointer :: dict
    !local variables
    logical :: first_option
    character(len=1) :: shortname
    character(len=max_field_length) :: name,default,help_string,conflicts
    type(dictionary), pointer :: help_dict

    !then check the options
    call f_zero(name)
    call f_zero(default)
    call f_zero(help_string)
    call f_zero(conflicts)
    call f_zero(shortname)
    call f_zero(first_option)
    nullify(help_dict)

    !compulsory
    name=dict .get. 'name'
    default=dict .get. 'default'
    help_string=dict .get. 'help_string'
    !optionals
    conflicts=dict .get. 'conflicts'
    shortname=dict .get. 'shortname'
    first_option=dict .get. 'first_option'
    if ('help_dict' .in. dict) call dict_copy(src=dict // 'help_dict',dest=help_dict)

    call yaml_cl_parse_option(parser,name,default,help_string,&
         shortname,help_dict,first_option,conflicts)

  end subroutine yaml_cl_parse_option_from_dict

  subroutine yaml_argparse(options,string)
    use dictionaries
    use f_utils, only: f_zero
    implicit none
    !> the dictionary of the options, should be nullified as input
    type(dictionary), pointer :: options
    !>definition of the input variables, given with a single string
    character(len=*), intent(in) :: string
    !local variables
    type(yaml_cl_parse) :: parser !< command line parser

    !define command-line options
    parser=yaml_cl_parse_null()
    call yaml_cl_parse_option(parser,string)
    !parse command line, and retrieve arguments
    call yaml_cl_parse_cmd_line(parser,args=options)
    !free command line parser information
    call yaml_cl_parse_free(parser)

  end subroutine yaml_argparse


  !> Fill the parsed dictionary with default values
  subroutine yaml_cl_parse_init(parser)
    use dictionaries
    implicit none
    !> the parser which has to be updated
    type(yaml_cl_parse), intent(inout) :: parser
    !local variable
    type(dictionary), pointer :: opt_iter
    character(len=max_field_length) :: default

    call dict_free(parser%args)
    call dict_init(parser%args)

    opt_iter=> dict_iter(parser%options)
    do while(associated(opt_iter))
       default=opt_iter//OPTDEFAULT
       if (trim(default) /= 'None') then
          call set(parser%args//dict_key(opt_iter),default)
       end if
       opt_iter=> dict_next(opt_iter)
    end do
  end subroutine yaml_cl_parse_init


  !> dump on stdout the long help for each of the options
  subroutine parser_help(parser,short)
    use dictionaries
    use yaml_output
    implicit none
    logical, intent(in) :: short
    type(yaml_cl_parse), intent(in) :: parser
    !local variables
    character(len=1) :: shortopt
    character(len=max_field_length) :: default
    type(dictionary), pointer :: iter

    call yaml_comment('Command-line arguments help',hfill='-',unit=6)
    !iterate on all the options
    iter => dict_iter(parser%options)
    do while(associated(iter))
       !find if this is the first option
       if (trim(parser%first_command_key) /= dict_key(iter)) then
          call option_help(iter)
       end if
       iter => dict_next(iter)
    end do
    call yaml_map('--help (-h)',&
         'Print this screen message (-h would produce a compact form)',&
         unit=6)
    if (len_trim(parser%first_command_key) /=0) then
       call option_help(parser%options//trim(parser%first_command_key))
       call yaml_comment('The above option can also be specified as a first argument')
    end if

    contains

      subroutine option_help(opt)
        implicit none
        type(dictionary), pointer :: opt
        !find if this is the first option
        if (OPTSNAME .in. opt) then
           shortopt=opt//OPTSNAME
           call yaml_mapping_open('--'//trim(adjustl(dict_key(opt)))//&
                ' (-'//shortopt//')',unit=6)
        else
           call yaml_mapping_open('--'//trim(adjustl(dict_key(opt))),&
                unit=6)
        end if
        if ((OPTHELP .in. opt) .and. .not. short) then
           call yaml_dict_dump(opt//OPTHELP,unit=6)
        else
           call yaml_dict_dump(opt//OPTSHELP,unit=6)
        end if
        default=opt//OPTDEFAULT
        if (trim(default) /= 'None') then
           call yaml_map('Default value (when not specified)',default,&
                unit=6)
        end if
        call yaml_mapping_close(unit=6)

      end subroutine option_help

  end subroutine parser_help


  !> routine for parsing the command line
  subroutine yaml_cl_parse_cmd_line(parser,args)
    use dictionaries
    use yaml_strings, only:f_strcpy
    use yaml_output
    !use yaml_output
    implicit none
    !> the parser which has to be updated
    !! resulting command line arguments are written in parser%args
    type(yaml_cl_parse), intent(inout) :: parser
    !> dictionary of the arguments.
    !! if present, the arguments are copied into this dictionary
    type(dictionary), pointer, intent(out), optional :: args
    !local variables
    integer :: icommands,ncommands
    type(dictionary), pointer :: dict,conf
    !fill the parser with default values
    call yaml_cl_parse_init(parser)

    !first see how many arguments are present
    ncommands=COMMAND_ARGUMENT_COUNT()

    icommands=1
    do while(icommands <= ncommands)
       call parse_command(dict)
       !cycle if command is unreadable
       if (f_err_check(err_id=ERROR_YAML_COMMAND_LINE_PARSER)) cycle
       !fill the parser with the parsed dictionary
       call dict_update(parser%args,dict)
       call dict_free(dict)
    end do

    !coherence checks.
    dict => dict_iter(parser%args)
    do while(associated(dict))
       !See if some option conflicts with others
       conf => dict_iter(parser%options//dict_key(dict) .get. OPTCONFL)
       do while(associated(conf))
          if (trim(dict_value(conf)) .in. parser%args) then
             call f_err_throw('The option "'//trim(dict_key(dict))//&
                  '" conflicts with the presence of "'//&
                  trim(dict_value(conf))//'".',&
                  err_id=ERROR_YAML_COMMAND_LINE_PARSER)
          end if
          conf => dict_next(conf)
       end do

       dict => dict_next(dict)
    end do

    if (present(args)) then
       nullify(args)
       call dict_copy(src=parser%args,dest=args)
    end if

    contains

      !> Parse the input command and returns the dictionary which is associated to it
      subroutine parse_command(dict)
        implicit none
        !> value of the key of the parser as it has been defined in the parser dictionary
        !! it can be 'help' or 'h' for dumping the help dictionary (in long or short version respectively)
        !> dictionary associated to the value (nullified on output if key absent)
        type(dictionary), pointer, intent(out) :: dict
        !local variables
        logical :: found
        integer :: ipos,jpos
        character(len=max_field_length) :: command,test,key
        character(len=1) :: short_key
        type(dictionary), pointer :: opt_iter

        nullify(dict)
        found=.false.
        call get_cmd(icommands,command)
        icommands=icommands+1
        !search for the long key value
        ipos=index(command,'--')
        if (ipos > 0 .and. len_trim(command(1:ipos-1)) == 0) then
           ipos=ipos+2
           !the help key is a particular case
           if (trim(command(ipos:)) == 'help') then
              if (icommands-1 /=1 .or. ncommands > 1) then
                 call f_err_throw('"help" keyword is only accepted as unique command.',&
                      err_id=ERROR_YAML_COMMAND_LINE_PARSER)
                 return
              end if
              call parser_help(parser,.false.)
              !usually, after the help screen is called, the code should exit
              call yaml_cl_parse_free(parser)
              call f_lib_finalize_noreport()
              stop
           end if
           !a long key has always to be specified as --key=yaml_dict
           jpos=index(command(ipos+1:),'=')
           if (jpos == 0) then
              call f_err_throw('A long key must be followed by its value. Write (--'//&
                   trim(command(ipos:))//'=<val>)',&
                   err_id=ERROR_YAML_COMMAND_LINE_PARSER)
              return
           end if
           !insert the key, it has to have the value associated to the name
           call f_strcpy(src=trim(command(ipos:ipos+jpos-1)),dest=key)

           !check if the key is among the known values
           if (key .notin. parser%options) then
              call f_err_throw('The key "'//trim(key)//&
                   '" is not among the known values.',err_id=ERROR_YAML_COMMAND_LINE_PARSER)
           end if

           !then parse the value as a yaml_string
           jpos=jpos+1
           if (len_trim(command(ipos+jpos:)) > 0) then
              dict=>yaml_load(trim(command(ipos+jpos:)),key)
           else
              call f_err_throw('Empty value for key "'//trim(key)//&
                   '"',err_id=ERROR_YAML_COMMAND_LINE_PARSER)
           end if
        else
           !search for short string format (one letter only)
           !example -k yaml_dict
           ipos=index(command,'-')
           if (ipos > 0 .and. len_trim(command(1:ipos-1)) == 0) then
              ipos=ipos+1
              short_key=command(ipos:ipos)
              if (len_trim(command(ipos+1:)) > 0) then
                 call f_err_throw('Short command line option '//&
                      trim(command(ipos-1:))//&
                      ' has to be on one letter only.',&
                      err_id=ERROR_YAML_COMMAND_LINE_PARSER)
                 return
              end if
              if (short_key == 'h') then
                 if (icommands-1 /=1 .or. ncommands > 1) then
                    call f_err_throw('"-h" option is only accepted as unique one.',&
                         err_id=ERROR_YAML_COMMAND_LINE_PARSER)
                    return
                 end if
                 call parser_help(parser,.true.)
                 !usually, after the help screen is called, the code should exit
                 call yaml_cl_parse_free(parser)
                 call f_lib_finalize_noreport()
                 stop
              end if
              !then fill the key according to the set values
              opt_iter=>dict_iter(parser%options)
              find_shkey: do while(associated(opt_iter))
                 if (OPTSNAME .in. opt_iter) then
                    test=opt_iter//OPTSNAME
                    if (trim(test)==short_key) then
                       found=.true.
                       call f_strcpy(src=dict_key(opt_iter),dest=key)
                       exit find_shkey
                    end if
                 end if
                 opt_iter=>dict_next(opt_iter)
              end do find_shkey
              if (.not. found) then
                 call f_err_throw('Unrecognized short command line option '//trim(command(ipos-1:)),&
                      err_id=ERROR_YAML_COMMAND_LINE_PARSER)
                 return
              end if
              if (trim(key) /= 'help') then
                 !then the command has to be increased and parsed
                 call get_cmd(icommands,command)
                 icommands=icommands+1
                 !then parse the value as a yaml_string
                 dict=>yaml_load(trim(command),key)
              end if
           else
              !only the first command can be given without option
              if (icommands-1 == 1 .and. len_trim(parser%first_command_key) > 0) then
                 key=parser%first_command_key
                 dict=>yaml_load(trim(command),key)
              else
                 call f_err_throw('Unrecognized option "'//&
                      trim(command)//'"',&
                      err_id=ERROR_YAML_COMMAND_LINE_PARSER)
                 return
              end if
           end if
        end if
      end subroutine parse_command

  end subroutine yaml_cl_parse_cmd_line


  !> Retrieve the command as a string
  subroutine get_cmd(icommands,command)
    use yaml_strings
    use dictionaries, only: f_err_throw
    implicit none
    integer, intent(in) :: icommands
    character(len=*), intent(out) :: command
    !local variables
    integer :: ierr

    command=repeat(' ',len(command))
    call get_command_argument(icommands,value=command,status=ierr)
    if (ierr/=0) then
       call f_err_throw('Error in get_command_argument, ierr='//trim(yaml_toa(ierr))//&
            ', command no.='//trim(yaml_toa(icommands)),err_id=ERROR_YAML_COMMAND_LINE_PARSER)
       command=repeat(' ',len(command))
    end if
  end subroutine get_cmd


  !> Define the errors related to the parsing of the yaml files
  subroutine yaml_parse_errors()
    use dictionaries
    implicit none

    call f_err_define(err_name='YAML_PARSE_ERROR',&
         err_msg='YAML parse error.',&
         err_action='Modify your inputs.',&
         err_id=YAML_PARSE_ERROR)

    !Define a dictionary to have a more verbosity of yaml_parse_error
    dict_yaml_errs => dict_new("<document start>" .is. &
         & "The first indentation is different at this line in front of the given key.", &
         &                     "mapping values" .is. &
         & "The indentation is different at this line." )

    call f_err_define(err_name='YAML_PARSE_UNSUPPORTED',&
         err_msg='YAML standard not supported.',&
         err_action='kindly ask developers to finish implementation.',&
         err_id=YAML_PARSE_UNSUPPORTED)

    call f_err_define(err_name='ERROR_YAML_COMMAND_LINE_PARSER',&
         err_msg='Error in yaml parsing of the command line',&
         err_action='Check the allowed options and their values (--help).',&
         err_id=ERROR_YAML_COMMAND_LINE_PARSER)

  end subroutine yaml_parse_errors


  !> Create a dict from a file (fname is the buffer containing all the file)
  subroutine yaml_parse_from_file(dict, fname)
    use dictionaries
    implicit none
    type(dictionary), pointer :: dict
    character(len = *), intent(in) :: fname

    integer(kind = 8) :: parser

    call yaml_parser_c_init(parser, fname, len(fname))
    dict => yaml_parse_(parser)
  end subroutine yaml_parse_from_file


  !> parse a database from the corresponding external routine
  !! such a routine is created from the Makefile macro
  !! yaml_db.mk which has to be included in the Makefile
  subroutine yaml_parse_database(dict,symbol)
    use f_precisions
    implicit none
    type(dictionary), pointer :: dict
    external :: symbol
    !local variables
    integer(f_integer) :: db_size
    integer(f_address) :: parser
    character, dimension(:), allocatable :: f_string !<not profiled in f_malloc

    !get reference to the database and size of the data
    !first get the size
    db_size=int(0,f_integer)
    call symbol(parser,db_size) !<parser here is ignored
    !allocate array
    allocate(f_string(db_size))
    call symbol(f_string,db_size) !<we fill the array
    call yaml_parse_from_char_array(dict,f_string)
    deallocate(f_string)
  end subroutine yaml_parse_database


  subroutine yaml_parse_from_char_array(dict, carr)
    use dictionaries
    implicit none
    type(dictionary), pointer :: dict
    character, dimension(:), intent(in) :: carr

    integer(kind = 8) :: parser

    call yaml_parser_c_init_from_buf(parser, carr(1), size(carr))
    dict => yaml_parse_(parser)
  end subroutine yaml_parse_from_char_array


  subroutine yaml_parse_from_string(dict, str)
    use dictionaries
    implicit none
    type(dictionary), pointer :: dict
    character(len = *), intent(in) :: str

    integer(kind = 8) :: parser

    call yaml_parser_c_init_from_buf(parser, str, len_trim(str))
    dict => yaml_parse_(parser)
  end subroutine yaml_parse_from_string

  function yaml_parse_(parser) result(output)
    use dictionaries
    implicit none
    integer(kind = 8), intent(in) :: parser

    type(dictionary), pointer :: dict, doc, output
    integer :: event, errid
    character(max_field_length) :: val

    ! Get event values from C.
    call yaml_parser_c_get_stream_start(STREAM_START)
    call yaml_parser_c_get_stream_end(STREAM_END)
    call yaml_parser_c_get_document_start(DOCUMENT_START)
    call yaml_parser_c_get_document_end(DOCUMENT_END)
    call yaml_parser_c_get_sequence_start(SEQUENCE_START)
    call yaml_parser_c_get_sequence_end(SEQUENCE_END)
    call yaml_parser_c_get_mapping_start(MAPPING_START)
    call yaml_parser_c_get_mapping_end(MAPPING_END)
    call yaml_parser_c_get_alias(ALIAS)
    call yaml_parser_c_get_scalar(SCALAR)
    call yaml_parser_c_get_error(ERROR)

    ! Initialise error if required.
    if (YAML_PARSE_ERROR == 0) then
       call f_err_define(err_name='YAML_PARSE_ERROR',&
            err_msg='YAML parse error.',&
            err_action='correct your YAML stream.',&
            err_id=YAML_PARSE_ERROR)
    end if
    if (YAML_PARSE_UNSUPPORTED == 0) then
       call f_err_define(err_name='YAML_PARSE_UNSUPPORTED',&
            err_msg='YAML standard not supported.',&
            err_action='kindly ask developers to finish implementation.',&
            err_id=YAML_PARSE_UNSUPPORTED)
    end if

    call f_err_open_try()
    call dict_init(dict)
    event = 0
    nullify(doc)
    val=repeat(' ',len(val))
    do while (event /= STREAM_END)
       call yaml_parser_c_next(parser, event, val, max_field_length)
       !print *,'event',event_toa(event),event,trim(val),'end'
       if (event == ERROR) then
          !search
          call f_err_throw(err_id = YAML_PARSE_ERROR, err_msg = trim(val))
          exit
       end if

       if (event == DOCUMENT_END) then
          if (.not.associated(doc)) call dict_init(doc) ! empty document case
          call add(dict, doc)
          nullify(doc)
       else if (event == MAPPING_START) then
          doc => build_map(parser)                      ! dictionary document case
       else if (event == SEQUENCE_START) then
          doc => build_seq(parser)                      ! list document case
       else if (event == SCALAR) then
          call dict_init(doc)                           ! scalar document case
          call set(doc, val)
       end if

       if (f_err_check(YAML_PARSE_ERROR)) exit
    end do
    errid = f_get_last_error(val)
    if (f_err_check(YAML_PARSE_ERROR)) then
       if (associated(doc)) call add(dict, doc)
    end if
    call f_err_close_try()

    if (event /= STREAM_END) call yaml_parser_c_finalize(parser)

    output => dict

    !if (errid == YAML_PARSE_ERROR) call f_err_throw(err_id = errid, err_msg = trim(val))
    if (errid == YAML_PARSE_ERROR) call yaml_parse_error_throw(val)

  contains

    !> Determine which is the event that has been recognized, to be used mostly for debugging purposes
    function event_toa(event) result(toa)
      implicit none
      integer, intent(in) :: event
      character(len=32) :: toa
      if(event==STREAM_START) then
         toa(1:len(toa))='STREAM_START'
      else if(event==STREAM_END) then
         toa(1:len(toa))='STREAM_END'
      else if(event==DOCUMENT_START) then
         toa(1:len(toa))='DOCUMENT_START'
      else if(event==DOCUMENT_END) then
         toa(1:len(toa))='DOCUMENT_END'
      else if(event==SEQUENCE_START) then
         toa(1:len(toa))='SEQUENCE_START'
      else if(event==SEQUENCE_END) then
         toa(1:len(toa))='SEQUENCE_END'
      else if(event==MAPPING_START) then
         toa(1:len(toa))='MAPPING_START'
      else if(event==MAPPING_END) then
         toa(1:len(toa))='MAPPING_END'
      else if(event==ALIAS) then
         toa(1:len(toa))='ALIAS'
      else if(event==SCALAR) then
         toa(1:len(toa))='SCALAR'
      else if(event==ERROR) then
         toa(1:len(toa))='ERROR'
      else
         toa(1:len(toa))='UNKNOWN'
      end if
    end function event_toa

  end function yaml_parse_


  recursive function build_map(parser) result(map)
    use dictionaries
    !use yaml_output
    implicit none
    integer(kind = 8), intent(in) :: parser

    type(dictionary), pointer :: m, sub, map
    integer :: event
    character(max_field_length) :: val, key

    call dict_init(m)
    map => m

    event = 0
    key(1:max_field_length) = " "
    do while (event /= STREAM_END)
       call yaml_parser_c_next(parser, event, val, max_field_length)
       !write(*,*) "map", event

       if (event == ERROR) then
          call f_err_throw(err_id = YAML_PARSE_ERROR, err_msg = trim(val))
          !call yaml_parse_error_throw(val)
          return
       end if

       if (event == MAPPING_END) then
          exit
       else if (event == MAPPING_START) then
          sub => build_map(parser)
          if (len_trim(key) == 0) stop "no key"
          call set(m // key, sub)
          key(1:max_field_length) = " "
       else if (event == SEQUENCE_START) then
          sub => build_seq(parser)
          if (len_trim(key) == 0) stop "no key"
          call set(m // key, sub)
          key(1:max_field_length) = " "
       else if (event == SCALAR) then
          if (len_trim(key) > 0) then
             ! This is a simple key / val entry.
             call set(m // key, val)
             key(1:max_field_length) = " "
          else
             ! We store a key for later usage.
             key = val
             !write(*,*) "set ", key
          end if
       else if (event == ALIAS) then
          call f_err_throw(err_id = YAML_PARSE_UNSUPPORTED, &
               & err_msg = "unsupported alias to " // trim(val))
          ! Fallback to stringified alias.
          if (len_trim(key) == 0) stop "no key"
          call set(m // key, "*" // trim(val))
          key(1:max_field_length) = " "
       end if

       if (f_err_check(YAML_PARSE_ERROR)) return

    end do
  end function build_map


  recursive function build_seq(parser) result(seq)
    use dictionaries
    implicit none
    integer(kind = 8), intent(in) :: parser

    type(dictionary), pointer :: s, sub, seq
    integer :: event
    character(max_field_length) :: val

    call dict_init(s)
    seq => s
    nullify(s)

    event = 0
    do while (event /= STREAM_END)
       call yaml_parser_c_next(parser, event, val, max_field_length)
       !write(*,*) "seq", event

       if (event == ERROR) then
          call f_err_throw(err_id = YAML_PARSE_ERROR, err_msg = trim(val))
          !call yaml_parse_error_throw(val)
          return
       end if

       if (event == SEQUENCE_END) then
          exit
       else if (event == MAPPING_START) then
          sub => build_map(parser)
          call add(seq, sub, s)
       else if (event == SEQUENCE_START) then
          sub => build_seq(parser)
          call add(seq, sub, s)
       else if (event == SCALAR) then
          call add(seq, val, s)
       else if (event == ALIAS) then
          call f_err_throw(err_id = YAML_PARSE_UNSUPPORTED, &
               & err_msg = "unsupported alias to " // trim(val))
          ! Fallback to stringified alias.
          call add(seq, "*" // trim(val), s)
       end if

       if (f_err_check(YAML_PARSE_ERROR)) return

    end do

  end function build_seq

  function yaml_load(string,key) result(dict)
    use dictionaries
    use yaml_strings, only: f_strcpy
    !use yaml_output !to be removed
    implicit none
    character(len=*), intent(in) :: string
    !> key of the parsed string. the dictionary will have this key
    !! if this variable is given
    character(len=*), intent(in), optional :: key
    type(dictionary), pointer :: dict
    !local variables
    type(dictionary), pointer :: loaded_string,test
    !parse from the given string
    call yaml_parse_from_string(loaded_string,string)

    !extract the first document
    dict => loaded_string .pop. 0

    call dict_free(loaded_string)

    if (present(key)) then !to be defined better
       select case(dict_value(dict))
       case(TYPE_DICT,TYPE_LIST)
          call dict_init(test)
          call dict_copy(test//trim(key),dict)
          call dict_free(dict)
          dict=>test
       case default
          test=>dict_new(trim(key) .is. dict_value(dict_iter(dict)))
          call dict_free(dict)
          dict=>test
       end select
    end if

  end function yaml_load

  !> Throw an error with YAML_PARSE_ERROR trying to give a better understandable message
  subroutine yaml_parse_error_throw(val)
    use dictionaries
    implicit none
    !Argument
    character(len=max_field_length) :: val
    !Local variables
    type(dictionary), pointer :: error
    !type(dictionary), pointer :: dict_error
    character(max_field_length) :: key,message
    !integer :: pp
    !dict_error=>dict_new()
    !pp = index(val,':')
    !key = val(1:pp-1)
    !message = val(pp+1:)
    !call set(dict_error//trim(key),trim(message))
    error=>dict_iter(dict_yaml_errs)
    message = trim(val)
    do while(associated(error))
      key = dict_key(error)
      if (index(val,trim(key)) /= 0) then
        !call set(dict_error//'Info',trim(dict_value(error)))
        message = trim(message) // '. ' // trim(dict_value(error))
        exit
      end if
      error=>dict_next(error)
    end do
    !call f_err_throw(err_id = YAML_PARSE_ERROR, err_dict = dict_error)
    !We add quite to have a yaml standard output.
    call f_err_throw(err_msg = '"'//trim(message)//'"',err_id = YAML_PARSE_ERROR)
  end subroutine yaml_parse_error_throw


  !> Nullify the dictionary dict_yaml_errs
  subroutine yaml_parse_errors_finalize()
     use dictionaries_base, only: dict_free
     implicit none
     call dict_free(dict_yaml_errs)

  end subroutine yaml_parse_errors_finalize


end module yaml_parse
