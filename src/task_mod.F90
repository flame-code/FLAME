!*****************************************************************************************
module mod_task
    implicit none
    !logical, allocatable:: comment_line(:)
    real(8):: time_start=0.d0
    real(8):: time_end=0.d0
    logical:: time_exceeded=.false.
    !-------------------------------------------------------
    type typ_file_ini
        integer:: nline_max=1000
        integer:: nline=-1
        character(256), allocatable:: file_lines(:)
        logical, allocatable:: stat_line_is_read(:)
        integer:: iline
        integer:: iline_header
        integer:: iline_next_header
        character(256):: str_par !first entry in the line (before equal sign)
        character(256):: str_val !next entries in the line (after equal sign)
        character(256):: str_tmp
    end type typ_file_ini
end module mod_task
!*****************************************************************************************
