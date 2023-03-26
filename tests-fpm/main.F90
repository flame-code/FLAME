!*****************************************************************************************
program tests
    !use iso_fortran_env, only: error_unit, output_unit
    use mod_flm_futile
    implicit none
    integer:: iunit, ierr
    call f_lib_initialize()
    iunit=f_get_free_unit(10**5)
    call yaml_set_stream(unit=iunit,filename='tests-fpm/flame_log.yaml',&
         record_length=92,istat=ierr,setdefault=.false.,tabbing=0,position='rewind')
     call yaml_release_document(iunit)
    call yaml_set_default_stream(iunit,ierr)
    call yaml_new_document()

    !---------------------------------------------------------------------------
    call test_cal_pot_gauss_s()
    !---------------------------------------------------------------------------
    call test_cal_pot_gauss_p()
    !---------------------------------------------------------------------------
    call test_cal_pot_r2gauss_s()
    !---------------------------------------------------------------------------
    call test_cent2_analytic()
    !---------------------------------------------------------------------------
    call test_optimize_basis_functions()
    !---------------------------------------------------------------------------
    !call test_get_basis_functions_cent2()
    !---------------------------------------------------------------------------
    call test_psolver_p3d()
    !---------------------------------------------------------------------------
    call test_symfunc_atom_behler()
    !---------------------------------------------------------------------------

    call f_lib_finalize()
end program tests
!*****************************************************************************************
