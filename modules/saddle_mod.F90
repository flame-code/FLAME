!*****************************************************************************************
module mod_saddle
    implicit none
    !private
    real(8):: ampl=-1.d0
    real(8):: dimsep=-1.d0
    !real(8):: alphax=-1.d0
    !real(8):: fnrmtol
    real(8):: beta=-1.d0
    real(8):: epotprime
    integer:: maxitec
    integer:: maxitsd
    integer:: maxitcg
    integer:: nmatr=-1
    integer:: moving_atoms_rand(10)
    integer:: nit !total number of iteration.
    logical:: ecconverged
    logical:: sdconverged
    logical:: cgconverged
    logical:: dmconverged
    logical:: do_elim_trans
    logical:: do_elim_rot
    character(256):: str_moving_atoms_rand
end module mod_saddle
!*****************************************************************************************
