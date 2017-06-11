!> @file
!! Define operations to insert and cite the bibliography
!! @author
!!    Copyright (C) 2014-2016 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
module f_bibliography
  use dictionaries
  implicit none

  private

  type(dictionary), pointer :: bibliography=>null()

  public :: f_bib_initialize,f_bib_finalize,f_bib_update,f_bib_item_exists
  public :: f_bib_get_item,get_bib_filename

  contains

    subroutine f_bib_initialize()
      implicit none
      call dict_init(bibliography)
    end subroutine f_bib_initialize

    subroutine f_bib_finalize()
      call dict_free(bibliography)
    end subroutine f_bib_finalize

    !> update bibliography database from an array of characters
    !! the char array might 
    subroutine f_bib_update(dict)
      implicit none
      type(dictionary), pointer :: dict
      call dict_copy(bibliography,dict)
    end subroutine f_bib_update

    function f_bib_item_exists(item) result(yes)
      implicit none
     character(len=*), intent(in) :: item
     logical :: yes

     yes=item .in. bibliography
   end function f_bib_item_exists

   function f_bib_get_item(item) result(dict)
     implicit none
     character(len=*), intent(in) :: item
     type(dictionary), pointer :: dict

     dict = bibliography .get. item

   end function f_bib_get_item

   pure subroutine get_bib_filename(stream_out,filename)
     use yaml_strings
     implicit none
     character(len=*), intent(in) :: stream_out
     character(len=*), intent(out) :: filename

     if (trim(stream_out) == 'stdout') then
        filename='citations.bib'
     else
        filename=stream_out
        call rstrip(filename,'.yaml')
        filename=trim(filename)//'.bib'
     end if
     
   end subroutine get_bib_filename

end module f_bibliography
