!> @file
!! Test the argument parsing of the library.e
!! @example yaml_test.f90
!! Extensive tests command line parsing from yaml
!! @author
!!    Copyright (C) 2015-2016 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


program yaml_argparse_main
  use futile
  implicit none
  character(len=*), parameter :: input1=&
       "  {name: ndim, shortname: n, default: 30,"//&
       "  help_string: Size of the simulation domain,"//&
       "  help_dict: {Allowed values: list of integers}}"
  character(len=*), parameter :: input2=&
       "  {name: geocode,shortname: g, default: P,"//&
       "  help_string: Boundary conditions,"//&
       "  help_dict: {Usage: set the boundary conditions of the run,"//&
       "              Allowed values: [F, S , W, P]}}"
  character(len=*), parameter :: input3=&       
       "  {name: angdeg, shortname: d, default: 90.0,"//&
       "  help_string: Degrees of the angles between the directions,"//&
       "  help_dict: {Allowed values: arrays of floats}}"
  character(len=*), parameter :: input4=&
       "  {name: input, shortname: i, default: None,"//&
       "  help_string: Inpufile of Poisson Solver,"//&
       "  help_dict: {Allowed values: dictionary in yaml format (mapping)}}"
  character(len=*), parameter :: input5=&
       "  {name: boldify, shortname: b, default: None,"//&
       "  help_string: Boldify the string as a test,"//&
       "  help_dict: {Allowed values: string scalar}}"

  character(len=*), parameter :: input6=&
       "  {name: blinkify, shortname: l, default: None,"//&
       "  help_string: Make the string blinking,"//&
       "  help_dict: {Allowed values: string scalar}}"


  character(len=*), parameter :: inputs=&
       '-'//input1//f_cr//&
       '-'//input2//f_cr//&
       '-'//input3//f_cr//&
       '-'//input4//f_cr//&
       '-'//input5//f_cr//&
       '-'//input6
  character(len=1) :: geocode
  character(len=32) :: bold,blink
  integer, dimension(3) :: ndims
  real(f_double), dimension(3) :: angdeg
  type(dictionary), pointer :: dict,options


  call f_lib_initialize()
  call yaml_new_document()

  nullify(dict)
  call f_zero(bold)
  call f_zero(blink)
  call yaml_argparse(options,inputs)
  call yaml_map('Commandline options provided',options)
  ndims=options//'ndim'
  geocode=options//'geocode'
  angdeg=options//'angdeg'
  bold=options .get. 'boldify'
  blink=options .get. 'blinkify'
!  input=options .get. 'input'
!  call dict_copy(dict,input) !null if absent

  if (len_trim(bold)>0) call yaml_map('Boldify test',yaml_bold(bold))
  if (len_trim(blink)>0) call yaml_map('Blinkify test',yaml_bold(yaml_blink(blink)))

  call dict_free(options)
!  call dict_free(dict)
  call f_lib_finalize()

end program yaml_argparse_main
