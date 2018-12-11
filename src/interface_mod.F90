!*****************************************************************************************
module mod_interface
    implicit none
interface
! ./src/acceleration.F90 :
subroutine acceleration(pressure,accpos,acclat,accvol,vpos,vlat,vvol,strten,fcart,latvec,amass,latmass,f0inv,md_type,nat) 
integer:: iat,i,j,md_type,nat
real(8),dimension(3,nat):: accpos,vpos,fcart,fpos
real(8),dimension(3,3)  :: acclat,vlat,latvec,tmplat,pressure,a,velmat,sigma,lattrans,latdottrans,gdot,g,ginv,gtot,str_matrix
real(8),dimension(3,3)  :: term1,term2,term3,term4,term5,term5_1,term5_2,sigmatrans,f0inv
real(8):: amass(nat),latmass,crossp(3),strten(6),vol,vpostmp(3),volvel,trace3
real(8):: accvol,vvol,vol_1_3
end subroutine acceleration
! ./src/ann_best_symfunc.F90 :
subroutine ann_best_symfunc(parini)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann
    use mod_symfunc, only: typ_symfunc_arr
    use mod_atoms, only: typ_atoms_arr
    type(typ_parini), intent(in):: parini
end subroutine ann_best_symfunc
subroutine cal_symfunc_diversity(n_tot,his,ann,disparity)
    use mod_ann, only: typ_ann
    integer, intent(in):: n_tot
    real(8), intent(in):: his(1000,n_tot)
    type(typ_ann), intent(inout):: ann
    real(8), intent(inout):: disparity
end subroutine cal_symfunc_diversity
subroutine gbounds_distro(ann,atoms_arr,strmess)
    use mod_atoms, only: typ_atoms_arr
    use mod_ann, only: typ_ann
    type(typ_ann), intent(inout):: ann
    type(typ_atoms_arr), intent(inout):: atoms_arr
    character(*), intent(in):: strmess
end subroutine gbounds_distro
! ./src/ann_check_symmetry_function.F90 :
subroutine ann_check_symmetry_function(parini)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_symfunc, only: typ_symfunc, typ_symfunc_arr
    use mod_atoms, only: typ_atoms_arr
    type(typ_parini), intent(in):: parini
end subroutine ann_check_symmetry_function
! ./src/ann_evaluate.F90 :
subroutine ann_evaluate_subtask(parini)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann
    use mod_symfunc, only: typ_symfunc_arr
    use mod_atoms, only: typ_atoms_arr
    type(typ_parini), intent(in):: parini
end subroutine ann_evaluate_subtask
! ./src/ann_gen_symmetry_function.F90 :
subroutine ann_gen_symmetry_function(parini)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr, ann_arr_deallocate
    use mod_symfunc, only: typ_symfunc, typ_symfunc_arr
    use mod_atoms, only: typ_atoms_arr
    type(typ_parini), intent(in):: parini
end subroutine ann_gen_symmetry_function
! ./src/ann_io.F90 :
subroutine read_input_ann(parini,iproc,ann_arr)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    type(typ_parini), intent(in):: parini
    integer, intent(in):: iproc
    type(typ_ann_arr), intent(inout):: ann_arr
end subroutine read_input_ann
subroutine read_symmetry_functions(parini,iproc,ifile,ann,rcut)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann
    type(typ_parini), intent(in):: parini
    integer, intent(in):: iproc, ifile
    type(typ_ann), intent(inout):: ann
    real(8), intent(out):: rcut
end subroutine read_symmetry_functions
subroutine set_radial_atomtype(parini,sat1,ityp)
    use mod_parini, only: typ_parini
    type(typ_parini), intent(in):: parini
    character(*), intent(in):: sat1
    integer, intent(out):: ityp(1)
end subroutine set_radial_atomtype
subroutine set_angular_atomtype(parini,sat1,sat2,ityp)
    use mod_parini, only: typ_parini
    type(typ_parini), intent(in):: parini
    character(*), intent(in):: sat1, sat2
    integer, intent(out):: ityp(2)
end subroutine set_angular_atomtype
subroutine write_ann_all(parini,ann_arr,iter)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(in):: ann_arr
    integer, intent(in):: iter
end subroutine write_ann_all
subroutine write_ann(parini,filename,ann)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann
    type(typ_parini), intent(in):: parini
    character(*):: filename
    type(typ_ann), intent(in):: ann
end subroutine write_ann
subroutine read_ann(parini,ann_arr)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
end subroutine read_ann
subroutine read_data_old(parini,filename_list,atoms_arr)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms_arr, atom_copy_old, set_rat_atoms
    type(typ_parini), intent(in):: parini
    character(*), intent(in):: filename_list
    type(typ_atoms_arr), intent(inout):: atoms_arr
end subroutine read_data_old
! ./src/ann_io_yaml.F90 :
subroutine read_input_ann_yaml(parini,iproc,ann_arr)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    type(typ_parini), intent(in):: parini
    integer, intent(in):: iproc
    type(typ_ann_arr), intent(inout):: ann_arr
end subroutine read_input_ann_yaml
subroutine get_symfunc_parameters_yaml(parini,iproc,fname,ann,rcut)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann
    type(typ_parini), intent(in):: parini
    type(typ_ann), intent(inout):: ann
    integer, intent(in):: iproc
    real(8), intent(out):: rcut
    character(50):: fname, sat1, sat2
end subroutine get_symfunc_parameters_yaml
subroutine write_ann_all_yaml(parini,ann_arr,iter)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(in):: ann_arr
    integer, intent(in):: iter
end subroutine write_ann_all_yaml
subroutine write_ann_yaml(parini,filename,ann,rcut)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann
    type(typ_parini), intent(in):: parini
    character(*):: filename
    type(typ_ann), intent(in):: ann
    real(8), intent(in):: rcut
end subroutine write_ann_yaml
subroutine read_ann_yaml(parini,ann_arr)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
end subroutine read_ann_yaml
subroutine set_dict_ann(ann,fname,stypat)
    use mod_ann, only: typ_ann
    type(typ_ann), intent(inout):: ann
    character(5):: stypat
    character(len=*):: fname 
end subroutine set_dict_ann
subroutine read_data_yaml(parini,filename_list,atoms_arr)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms_arr, atom_allocate_old, atom_deallocate, atom_copy_old
    type(typ_parini), intent(in):: parini
    character(*), intent(in):: filename_list
    type(typ_atoms_arr), intent(inout):: atoms_arr
end subroutine read_data_yaml
! ./src/ann_pot_atom.F90 :
subroutine cal_ann_atombased(parini,atoms,symfunc,ann_arr)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_ann, only: typ_ann_arr, convert_ann_epotd
    use mod_symfunc, only: typ_symfunc
    use mod_linked_lists, only: typ_pia_arr
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_symfunc), intent(inout):: symfunc
end subroutine cal_ann_atombased
! ./src/ann_pot_cent1.F90 :
subroutine cal_ann_cent1(parini,atoms,symfunc,ann_arr)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_ann, only: typ_ann_arr, convert_ann_epotd
    use mod_symfunc, only: typ_symfunc
    use mod_electrostatics, only: typ_poisson
    use mod_linked_lists, only: typ_pia_arr
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_symfunc), intent(inout):: symfunc
end subroutine cal_ann_cent1
subroutine get_qat_from_chi_cent1(parini,ann_arr,atoms,poisson,a)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_atoms, only: typ_atoms
    use mod_electrostatics, only: typ_poisson
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_atoms), intent(inout):: atoms
    type(typ_poisson), intent(inout):: poisson
    real(8), intent(inout):: a(atoms%nat+1,atoms%nat+1)
end subroutine get_qat_from_chi_cent1
subroutine get_qat_from_chi_dir(parini,ann_arr,atoms,a)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_atoms, only: typ_atoms
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_atoms), intent(inout):: atoms
    real(8), intent(inout):: a(atoms%nat+1,atoms%nat+1)
end subroutine get_qat_from_chi_dir
subroutine init_electrostatic_cent1(parini,atoms,ann_arr,a,poisson)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, update_ratp
    use mod_ann, only: typ_ann_arr
    use mod_electrostatics, only: typ_poisson
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
    type(typ_ann_arr), intent(inout):: ann_arr
    real(8), intent(inout):: a(atoms%nat+1,atoms%nat+1)
    type(typ_poisson), intent(inout):: poisson
end subroutine init_electrostatic_cent1
subroutine get_amat_cent1(atoms,ann_arr,a)
    use mod_atoms, only: typ_atoms, update_ratp
    use mod_ann, only: typ_ann_arr
    type(typ_atoms), intent(inout):: atoms
    type(typ_ann_arr), intent(inout):: ann_arr
    real(8), intent(inout):: a(atoms%nat+1,atoms%nat+1)
end subroutine get_amat_cent1
subroutine fini_electrostatic_cent1(parini,atoms,poisson)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_electrostatics, only: typ_poisson
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
    type(typ_poisson), intent(inout):: poisson
end subroutine fini_electrostatic_cent1
subroutine get_electrostatic_cent1(parini,atoms,ann_arr,epot_c,a,poisson)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_ann, only: typ_ann_arr
    use mod_electrostatics, only: typ_poisson
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
    type(typ_ann_arr), intent(inout):: ann_arr
    real(8), intent(out):: epot_c
    real(8), intent(inout):: a(atoms%nat+1,atoms%nat+1)
    type(typ_poisson), intent(inout):: poisson
end subroutine get_electrostatic_cent1
subroutine cal_electrostatic_ann(parini,atoms,ann_arr,a,poisson)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, update_ratp
    use mod_ann, only: typ_ann_arr
    use mod_electrostatics, only: typ_poisson
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
    type(typ_ann_arr), intent(inout):: ann_arr
    real(8), intent(in):: a(atoms%nat+1,atoms%nat+1)
    type(typ_poisson), intent(inout):: poisson
end subroutine cal_electrostatic_ann
subroutine charge_analysis(parini,atoms,ann_arr)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_ann, only: typ_ann_arr
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(in):: atoms
    type(typ_ann_arr), intent(inout):: ann_arr
end subroutine charge_analysis
subroutine get_qat_from_chi_iter(parini,ann_arr,atoms,a)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_atoms, only: typ_atoms
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_atoms), intent(inout):: atoms
    real(8), intent(inout):: a(atoms%nat+1,atoms%nat+1)
end subroutine get_qat_from_chi_iter
subroutine get_ener_gradient_cent1(parini,poisson,ann_arr,atoms,g,qtot)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_atoms, only: typ_atoms, update_ratp
    use mod_electrostatics, only: typ_poisson
    type(typ_parini), intent(in):: parini
    type(typ_poisson),intent(inout):: poisson
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_atoms), intent(inout):: atoms
    real(8), intent(out):: g(atoms%nat), qtot
end subroutine get_ener_gradient_cent1
subroutine get_qat_from_chi_operator(parini,poisson,ann_arr,atoms)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_atoms, only: typ_atoms, set_qat, update_ratp
    use mod_electrostatics, only: typ_poisson
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_atoms), intent(inout):: atoms
    type(typ_poisson),intent(inout):: poisson
end subroutine get_qat_from_chi_operator
! ./src/ann_pot_cent2.F90 :
subroutine cal_ann_cent2(parini,atoms,symfunc,ann_arr)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, update_ratp
    use mod_ann, only: typ_ann_arr, typ_cent, convert_ann_epotd
    use mod_symfunc, only: typ_symfunc
    use mod_linked_lists, only: typ_pia_arr
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_symfunc), intent(inout):: symfunc
end subroutine cal_ann_cent2
subroutine get_qat_from_chi_cent2(parini,ann_arr,atoms,cent)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr, typ_cent
    use mod_atoms, only: typ_atoms
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_atoms), intent(inout):: atoms
    type(typ_cent), intent(inout):: cent
end subroutine get_qat_from_chi_cent2
subroutine init_cent2(parini,ann_arr,atoms,cent)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr, typ_cent
    use mod_atoms, only: typ_atoms
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_atoms), intent(inout):: atoms
    type(typ_cent), intent(inout):: cent
end subroutine init_cent2
subroutine final_cent2(cent)
    use mod_ann, only: typ_cent
    type(typ_cent), intent(inout):: cent
end subroutine final_cent2
subroutine cent2_force(parini,ann_arr,atoms,cent)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr, typ_cent
    use mod_atoms, only: typ_atoms
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_atoms), intent(inout):: atoms
    type(typ_cent), intent(inout):: cent
end subroutine cent2_force
subroutine cal_potential_cent2(parini,ann_arr,atoms,cent)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr, typ_cent
    use mod_atoms, only: typ_atoms
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_atoms), intent(inout):: atoms
    type(typ_cent), intent(inout):: cent
end subroutine cal_potential_cent2
subroutine cal_cent2_pot_pairsum(parini,ann_arr,atoms,cent,epot_es)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_ann, only: typ_ann_arr, typ_cent
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_atoms), intent(inout):: atoms
    type(typ_cent), intent(inout):: cent
    real(8), intent(inout):: epot_es
end subroutine cal_cent2_pot_pairsum
subroutine cal_cent2_pairsum_force(parini,ann_arr,atoms,cent)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_ann, only: typ_ann_arr, typ_cent
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_atoms), intent(inout):: atoms
    type(typ_cent), intent(inout):: cent
end subroutine cal_cent2_pairsum_force
subroutine cal_cent2_pot_bps(parini,ann_arr,atoms,cent,epot_es)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_ann, only: typ_ann_arr, typ_cent
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_atoms), intent(inout):: atoms
    type(typ_cent), intent(inout):: cent
    real(8), intent(inout):: epot_es
end subroutine cal_cent2_pot_bps
subroutine put_cent2_gauss_to_grid(parini,atoms,cent)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_ann, only: typ_cent
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(in):: atoms
    type(typ_cent), intent(inout):: cent
end subroutine put_cent2_gauss_to_grid
subroutine cal_cent2_shortrange_ewald(parini,ann_arr,atoms,cent,epot_es)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_atoms, only: typ_atoms
    use mod_ann, only: typ_ann_arr, typ_cent
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(in):: ann_arr
    type(typ_atoms), intent(in):: atoms
    type(typ_cent), intent(inout):: cent
    real(8), intent(inout):: epot_es
end subroutine cal_cent2_shortrange_ewald
subroutine cal_shortrange_ewald_force_cent2(parini,ann_arr,atoms,cent)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr, typ_cent
    use mod_atoms, only: typ_atoms
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(in):: ann_arr
    type(typ_atoms), intent(inout):: atoms
    type(typ_cent), intent(inout):: cent
end subroutine cal_shortrange_ewald_force_cent2
subroutine erf_over_r_taylor(r,funcval,funcval_der)
    real(8), intent(in):: r
    real(8), intent(out):: funcval, funcval_der
end subroutine erf_over_r_taylor
subroutine calc_multipoles_cent2(parini,atoms,poisson,rel)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_electrostatics, only: typ_poisson
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(in):: atoms
    real(8), intent(in):: rel(3,atoms%nat)
    type(typ_poisson), intent(inout):: poisson
end subroutine calc_multipoles_cent2
subroutine calc_multipoles_grid_cent2(parini,atoms,poisson)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_electrostatics, only: typ_poisson
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(in):: atoms
    type(typ_poisson), intent(inout):: poisson
end subroutine calc_multipoles_grid_cent2
! ./src/ann_pot_cent3.F90 :
subroutine cal_ann_cent3(parini,atoms,symfunc,ann_arr)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, update_ratp
    use mod_ann, only: typ_ann_arr, typ_cent, convert_ann_epotd
    use mod_symfunc, only: typ_symfunc
    use mod_linked_lists, only: typ_pia_arr
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_symfunc), intent(inout):: symfunc
end subroutine cal_ann_cent3
subroutine get_dqat_from_chi_dir_cent3(parini,ann_arr,atoms,cent,a)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr, typ_cent
    use mod_atoms, only: typ_atoms
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_atoms), intent(inout):: atoms
    type(typ_cent), intent(in):: cent
    real(8), intent(inout):: a(atoms%nat+1,atoms%nat+1)
end subroutine get_dqat_from_chi_dir_cent3
subroutine get_qat_from_chi_cent3(parini,ann_arr,atoms,cent)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr, typ_cent
    use mod_atoms, only: typ_atoms
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_atoms), intent(inout):: atoms
    type(typ_cent), intent(inout):: cent
end subroutine get_qat_from_chi_cent3
subroutine get_qat_from_chi_dir_cent3(parini,ann_arr,atoms,cent,a)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr, typ_cent
    use mod_atoms, only: typ_atoms
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_atoms), intent(inout):: atoms
    type(typ_cent), intent(in):: cent
    real(8), intent(inout):: a(atoms%nat+1,atoms%nat+1)
end subroutine get_qat_from_chi_dir_cent3
subroutine get_qat_from_chi_operator_cent3(parini,ann_arr,atoms,cent)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr, typ_cent
    use mod_atoms, only: typ_atoms
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_atoms), intent(inout):: atoms
    type(typ_cent), intent(inout):: cent
end subroutine get_qat_from_chi_operator_cent3
subroutine init_cent3(parini,ann_arr,atoms,cent)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr, typ_cent
    use mod_atoms, only: typ_atoms
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_atoms), intent(inout):: atoms
    type(typ_cent), intent(inout):: cent
end subroutine init_cent3
subroutine final_cent3(cent)
    use mod_ann, only: typ_cent
    type(typ_cent), intent(inout):: cent
end subroutine final_cent3
subroutine cent3_force(parini,ann_arr,atoms,cent)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr, typ_cent
    use mod_atoms, only: typ_atoms
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_atoms), intent(inout):: atoms
    type(typ_cent), intent(inout):: cent
end subroutine cent3_force
subroutine cal_potential_cent3(parini,ann_arr,atoms,cent)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr, typ_cent
    use mod_atoms, only: typ_atoms
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_atoms), intent(inout):: atoms
    type(typ_cent), intent(inout):: cent
end subroutine cal_potential_cent3
subroutine cal_cent3_pot_pairsum(parini,ann_arr,atoms,cent,epot_es)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_ann, only: typ_ann_arr, typ_cent
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_atoms), intent(inout):: atoms
    type(typ_cent), intent(inout):: cent
    real(8), intent(inout):: epot_es
end subroutine cal_cent3_pot_pairsum
subroutine cal_cent3_pairsum_force(parini,ann_arr,atoms,cent)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_ann, only: typ_ann_arr, typ_cent
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_atoms), intent(inout):: atoms
    type(typ_cent), intent(inout):: cent
end subroutine cal_cent3_pairsum_force
subroutine cal_cent3_pot_bps(parini,ann_arr,atoms,cent,epot_es)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_ann, only: typ_ann_arr, typ_cent
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_atoms), intent(inout):: atoms
    type(typ_cent), intent(inout):: cent
    real(8), intent(inout):: epot_es
end subroutine cal_cent3_pot_bps
subroutine put_cent3_gauss_to_grid(parini,atoms,cent)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_ann, only: typ_cent
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(in):: atoms
    type(typ_cent), intent(inout):: cent
end subroutine put_cent3_gauss_to_grid
subroutine cal_cent3_shortrange_ewald(parini,ann_arr,atoms,cent,epot_es)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_atoms, only: typ_atoms
    use mod_ann, only: typ_ann_arr, typ_cent
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(in):: ann_arr
    type(typ_atoms), intent(in):: atoms
    type(typ_cent), intent(inout):: cent
    real(8), intent(inout):: epot_es
end subroutine cal_cent3_shortrange_ewald
subroutine cal_shortrange_ewald_force_cent3(parini,ann_arr,atoms,cent)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr, typ_cent
    use mod_atoms, only: typ_atoms
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(in):: ann_arr
    type(typ_atoms), intent(inout):: atoms
    type(typ_cent), intent(inout):: cent
end subroutine cal_shortrange_ewald_force_cent3
! ./src/ann_pot_cent_common.F90 :
subroutine cal_force_chi_part1(parini,symfunc,iat,atoms,out_ann,ann_arr)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_symfunc, only: typ_symfunc
    use mod_atoms, only: typ_atoms
    type(typ_parini), intent(in):: parini
    type(typ_symfunc), intent(in):: symfunc
    integer, intent(in):: iat
    type(typ_atoms), intent(in):: atoms
    real(8), intent(in):: out_ann
    type(typ_ann_arr), intent(inout):: ann_arr
end subroutine cal_force_chi_part1
subroutine cal_force_chi_part2(parini,symfunc,atoms,ann_arr)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_symfunc, only: typ_symfunc
    use mod_atoms, only: typ_atoms
    type(typ_parini), intent(in):: parini
    type(typ_symfunc), intent(in):: symfunc
    type(typ_atoms), intent(inout):: atoms
    type(typ_ann_arr), intent(inout):: ann_arr
end subroutine cal_force_chi_part2
subroutine repulsive_potential_cent(parini,atoms,ann_arr)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, update_ratp
    use mod_ann, only: typ_ann_arr
    use mod_linked_lists, only: typ_linked_lists
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
    type(typ_ann_arr), intent(inout):: ann_arr
end subroutine repulsive_potential_cent
! ./src/ann_pot_main.F90 :
subroutine get_fcn_ann(parini,idp,str_dataset,ann_arr,opt_ann,fcn_ann,fcn_ref)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_opt_ann, only: typ_opt_ann, set_opt_ann_grad
    use mod_atoms, only: typ_atoms, atom_copy_old, update_ratp
    type(typ_parini), intent(in):: parini
    integer, intent(in):: idp
    character(*), intent(in):: str_dataset
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_opt_ann), intent(inout):: opt_ann
    real(8), intent(out):: fcn_ann
    real(8), intent(out):: fcn_ref
end subroutine get_fcn_ann
subroutine cal_ann_main(parini,atoms,symfunc,ann_arr,opt_ann)
    use mod_tightbinding, only: typ_partb
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_ann, only: typ_ann_arr
    use mod_symfunc, only: typ_symfunc
    use mod_opt_ann, only: typ_opt_ann
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_symfunc), intent(inout):: symfunc
    type(typ_opt_ann), intent(inout):: opt_ann
end subroutine cal_ann_main
subroutine prefit_cent_ener_ref(parini,ann_arr,symfunc_train,symfunc_valid,atoms_train,atoms_valid,opt_ann)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_symfunc, only: typ_symfunc_arr
    use mod_opt_ann, only: typ_opt_ann, convert_opt_x_ann_arr
    use mod_atoms, only: typ_atoms, typ_atoms_arr, atom_copy_old
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_symfunc_arr), intent(inout):: symfunc_train, symfunc_valid
    type(typ_atoms_arr), intent(inout):: atoms_train
    type(typ_atoms_arr), intent(inout):: atoms_valid
    type(typ_opt_ann), intent(inout):: opt_ann
end subroutine prefit_cent_ener_ref
subroutine prefit_cent(parini,ann_arr,symfunc_train,symfunc_valid,atoms_train,atoms_valid,opt_ann)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_symfunc, only: typ_symfunc_arr
    use mod_opt_ann, only: typ_opt_ann, convert_opt_x_ann_arr
    use mod_atoms, only: typ_atoms, typ_atoms_arr, atom_copy_old
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_symfunc_arr), intent(inout):: symfunc_train, symfunc_valid
    type(typ_atoms_arr), intent(inout):: atoms_train
    type(typ_atoms_arr), intent(inout):: atoms_valid
    type(typ_opt_ann), intent(inout):: opt_ann
end subroutine prefit_cent
! ./src/ann_pot_tb.F90 :
subroutine cal_ann_tb(parini,partb,atoms,ann_arr,symfunc,opt_ann)
    use mod_parini, only: typ_parini
    use mod_tightbinding, only: typ_partb
    use mod_potl, only: potl_typ
    use mod_atoms, only: typ_atoms, update_ratp
    use mod_ann, only: typ_ann_arr, convert_ann_epotd
    use mod_symfunc, only: typ_symfunc
    use mod_opt_ann, only: typ_opt_ann, set_opt_ann_grad
    use mod_linked_lists, only: typ_pia_arr, typ_linked_lists
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
    type(typ_opt_ann), intent(inout):: opt_ann
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_symfunc), intent(inout):: symfunc
    type(typ_partb), intent(inout):: partb
    real(8):: hgen_der(4,1:atoms%nat,1:atoms%nat)  , ttxyz !derivative of 
end subroutine cal_ann_tb
subroutine lenoskytb_ann(parini,ann_arr,pia_arr,linked_lists,partb,atoms,natsi,count_md)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_atoms, only: typ_atoms
    use mod_potl, only: potl_typ
    use mod_linked_lists, only: typ_pia_arr, typ_linked_lists
    use mod_tightbinding, only: typ_partb, lenosky
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_linked_lists), intent(in):: linked_lists
    type(typ_pia_arr), intent(in):: pia_arr
    type(typ_partb), intent(inout):: partb
    type(typ_atoms), intent(inout):: atoms
    integer, intent(in):: natsi
    real(8), intent(inout):: count_md
end subroutine lenoskytb_ann
subroutine fit_hgen(parini,ann_arr,opt_ann)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, atom_allocate_old, set_rat_iat, update_ratp
    use mod_ann, only: typ_ann_arr, convert_x_ann
    use mod_symfunc, only: typ_symfunc
    use mod_opt_ann, only: typ_opt_ann, get_opt_ann_x, set_opt_ann_x
    use mod_parlm, only: typ_parlm
    type(typ_parini), intent(in):: parini
    type(typ_opt_ann), intent(inout):: opt_ann
    type(typ_ann_arr), intent(inout):: ann_arr
end subroutine fit_hgen
subroutine fcn_hgen(m,n,x,fvec,fjac,ldfjac,iflag,iann,ann_arr,hgen_ltb,yall)
    use mod_ann, only: typ_ann_arr, convert_ann_epotd
    integer, intent(in):: m, n, ldfjac, iflag, iann
    type(typ_ann_arr), intent(inout):: ann_arr
    real(8), intent(in):: x(n), hgen_ltb(4,325), yall(ann_arr%ann(1)%nn(0),1000,325)
    real(8), intent(inout):: fvec(m), fjac(m,n)
end subroutine fcn_hgen
! ./src/ann_process.F90 :
subroutine cal_architecture(ann,epot)
    use mod_ann, only: typ_ann
    type(typ_ann), intent(inout):: ann
    real(8), intent(inout):: epot
end subroutine cal_architecture
subroutine cal_architecture_1hiddenlayer(ann,epot)
    use mod_ann, only: typ_ann
    type(typ_ann), intent(inout):: ann
    real(8), intent(inout):: epot
end subroutine cal_architecture_1hiddenlayer
subroutine cal_architecture_2hiddenlayer(ann,epot)
    use mod_ann, only: typ_ann
    type(typ_ann), intent(inout):: ann
    real(8), intent(inout):: epot
end subroutine cal_architecture_2hiddenlayer
subroutine cal_architecture_der(ann,epot)
    use mod_ann, only: typ_ann
    type(typ_ann), intent(inout):: ann
    real(8), intent(inout):: epot
end subroutine cal_architecture_der
subroutine cal_architecture_der_1hiddenlayer(ann,epot)
    use mod_ann, only: typ_ann
    type(typ_ann), intent(inout):: ann
    real(8), intent(inout):: epot
end subroutine cal_architecture_der_1hiddenlayer
subroutine cal_architecture_der_2hiddenlayer(ann,epot)
    use mod_ann, only: typ_ann
    type(typ_ann), intent(inout):: ann
    real(8), intent(inout):: epot
end subroutine cal_architecture_der_2hiddenlayer
! ./src/ann_symfunc_atom_behler.F90 :
subroutine symmetry_functions_driver(parini,ann_arr,atoms,symfunc)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_symfunc, only: typ_symfunc
    use mod_atoms, only: typ_atoms
    use mod_linked_lists, only: typ_pia_arr !,typ_linked_lists
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_atoms), intent(in):: atoms
    type(typ_symfunc), intent(inout):: symfunc
end subroutine symmetry_functions_driver
subroutine symmetry_functions_g02_atom(ann_arr,pia,ib,iat,isat,jsat,symfunc)
    use mod_ann, only: typ_ann_arr
    use mod_symfunc, only: typ_symfunc
    use mod_linked_lists, only: typ_pia
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_pia), intent(in):: pia
    integer, intent(in):: ib, iat, isat, jsat
    type(typ_symfunc), intent(inout):: symfunc
end subroutine symmetry_functions_g02_atom
subroutine symmetry_functions_g04_atom(ann_arr,isat,iat,jsat,jat_maincell,ksat,kat_maincell,rij,rik,rjk,drij,drik,drjk,fcij,fcdij,fcik,fcdik,fcjk,fcdjk,symfunc)
    use mod_ann, only: typ_ann_arr
    use mod_symfunc, only: typ_symfunc
    type(typ_ann_arr), intent(inout):: ann_arr
    integer, intent(in):: isat, iat, jsat, jat_maincell, ksat, kat_maincell
    real(8), intent(in):: rij, rik, rjk, drij(3), drik(3), drjk(3), fcij, fcdij, fcik, fcdik, fcjk, fcdjk
    type(typ_symfunc), intent(inout):: symfunc
end subroutine symmetry_functions_g04_atom
subroutine symmetry_functions_g05_atom(ann_arr,piaij,piaik,ibij,ibik,iat,isat,jsat,ksat,symfunc)
    use mod_ann, only: typ_ann_arr
    use mod_symfunc, only: typ_symfunc
    use mod_linked_lists, only: typ_pia
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_pia), intent(in):: piaij, piaik
    integer, intent(in):: ibij, ibik, isat, iat, jsat, ksat
    type(typ_symfunc), intent(inout):: symfunc
end subroutine symmetry_functions_g05_atom
subroutine symmetry_functions_g06_atom(ann,iat,jat_maincell,r,dr,fc,fcd)
    use mod_ann, only: typ_ann
    type(typ_ann), intent(inout):: ann
    integer, intent(in):: iat, jat_maincell
    real(8), intent(in):: r, dr(3), fc, fcd
end subroutine symmetry_functions_g06_atom
function cutoff_function(r, rc) result(fc)
    real(8), intent(in):: r, rc
    real(8):: fc, pi
end function cutoff_function
function cutoff_function_der(r, rc) result(fcd)
    real(8), intent(in):: r, rc
    real(8):: fcd, pi
end function cutoff_function_der
subroutine symmetry_functions_g05_atom2(ann_arr,piaij,piaik,ibij,ibik,iat,isat,jsat,ksat,symfunc)
    use mod_ann, only: typ_ann_arr
    use mod_symfunc, only: typ_symfunc
    use mod_linked_lists, only: typ_pia
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_pia), intent(in):: piaij, piaik
    integer, intent(in):: ibij, ibik, isat, iat, jsat, ksat
    type(typ_symfunc), intent(inout):: symfunc
end subroutine symmetry_functions_g05_atom2
! ./src/ann_symfunc_atom_stefan.F90 :
subroutine symmetry_functions_driver_stefan(parini,ann_arr,atoms,symfunc)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_symfunc, only: typ_symfunc
    use mod_atoms, only: typ_atoms, set_rcov, get_rat
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_atoms), intent(inout):: atoms
    type(typ_symfunc), intent(inout):: symfunc
end subroutine symmetry_functions_driver_stefan
subroutine fingerprint_power(nat, npl, alat, rxyz, rcov, fpall)
  integer, parameter:: natx_sphere=500,lseg=4
  integer, parameter:: nwork=100
  real(8):: fpall(4,npl,nat),amplitude(natx_sphere)
  real(8):: rxyz(3,nat),rcov(nat)
  real(8):: alat(3, 3),alatalat(3,3),eigalat(3),aa(4,4),aaev(4)
end subroutine fingerprint_power
subroutine fingerprint_periodic(nat, natx_sphere, lseg, alat, rxyz, rcov, fpall)
end subroutine fingerprint_periodic
          subroutine Wblock(nat_sphere,lseg,om)
end subroutine wblock
    subroutine mltampl_4_alborz(nat,amplitude,om)
end subroutine mltampl_4_alborz
    subroutine mltampl_1_alborz(nat,amplitude,om)
end subroutine mltampl_1_alborz
subroutine frac2cart(nat, alat, xyzred, rxyz, convert)
end subroutine frac2cart
         subroutine create_om_1_alborz(nat,rxyz,rcov,om)
end subroutine create_om_1_alborz
         subroutine create_om_4_alborz(nat,rxyz,rcov,om)
end subroutine create_om_4_alborz
         subroutine OLDcreate_om_4(nat,rxyz,rcov,om)
end subroutine oldcreate_om_4
         subroutine create_om(nat,rxyz,rcov,nid,om)
end subroutine create_om
! ./src/ann_symfunc_main.F90 :
subroutine symmetry_functions(parini,ann_arr,atoms,symfunc,apply_gbounds)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_symfunc, only: typ_symfunc
    use mod_atoms, only: typ_atoms
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_atoms), intent(inout):: atoms
    type(typ_symfunc), intent(inout):: symfunc
    logical, intent(in):: apply_gbounds
end subroutine symmetry_functions
! ./src/ann_symfunc_pair_behler.F90 :
subroutine symmetry_functions_driver_bond(parini,ann_arr,atoms,symfunc)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_symfunc, only: typ_symfunc
    use mod_atoms, only: typ_atoms
    use mod_linked_lists, only: typ_pia_arr !,typ_linked_lists
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_atoms), intent(in):: atoms
    type(typ_symfunc), intent(inout):: symfunc
end subroutine symmetry_functions_driver_bond
subroutine symmetry_functions_driver_bond_tmp(ann_arr,atoms)
    use mod_ann, only: typ_ann_arr
    use mod_atoms, only: typ_atoms, get_rat
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_atoms), intent(in):: atoms
    integer, parameter::lwork=100
end subroutine symmetry_functions_driver_bond_tmp
subroutine symmetry_functions_g01_bond(ann_arr,ib,pia,symfunc)
    use mod_linked_lists, only: typ_pia
    use mod_ann, only: typ_ann_arr
    use mod_symfunc, only: typ_symfunc
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_pia), intent(in):: pia
    integer, intent(in):: ib
    type(typ_symfunc), intent(inout):: symfunc
end subroutine symmetry_functions_g01_bond
subroutine symmetry_functions_g02_bond(ann_arr,iat,jat,rij,drij,fcij,fcdij,rik,drik,fcik,fcdik,rjk,drjk,fcjk,fcdjk)
    use mod_ann, only: typ_ann_arr
    type(typ_ann_arr), intent(inout):: ann_arr
    integer, intent(in):: iat, jat
    real(8), intent(in):: fcij, fcik, fcjk, fcdij, fcdik, fcdjk
    real(8), intent(in):: rij, rik, rjk
    real(8), intent(in):: drij(3), drik(3), drjk(3)
end subroutine symmetry_functions_g02_bond
subroutine symmetry_functions_g04_bond(ann_arr,iat,jat,rij,drij,fcij,fcdij,rik,drik,fcik,fcdik,rjk,drjk,fcjk,fcdjk)
    use mod_ann, only: typ_ann_arr
    type(typ_ann_arr), intent(inout):: ann_arr
    integer, intent(in):: iat, jat
    real(8), intent(in):: fcij, fcik, fcjk, fcdij, fcdik, fcdjk
    real(8), intent(in):: rij, rik, rjk
    real(8), intent(in):: drij(3), drik(3), drjk(3)
end subroutine symmetry_functions_g04_bond
! ./src/atoms_minhocao.F90 :
subroutine atmdata(amu,rcov,symbol,znucl)
 real(8),intent(in) :: znucl
 real(8),intent(out) :: amu,rcov
 character(len=2),intent(out) :: symbol
end subroutine atmdata
subroutine mlj_atmdata(amu,sigma,eps,rcov,symbol,znucl)
 real(8),intent(in) :: znucl
 real(8),intent(out) :: amu,rcov,sigma,eps
 character(len=2),intent(out) :: symbol
end subroutine mlj_atmdata
subroutine symbol2znucl(amu,rcov,symbol,znucl)
 real(8),intent(out) :: znucl
 real(8),intent(out) :: amu,rcov
 character(len=2),intent(in) :: symbol
end subroutine symbol2znucl
subroutine mlj_symbol2znucl(amu,sigma,eps,rcov,symbol,znucl)
 real(8),intent(out) :: znucl
 real(8),intent(out) :: amu,rcov,sigma,eps
 character(len=2),intent(in) :: symbol
end subroutine mlj_symbol2znucl
! ./src/bader_mod.F90 :
! ./src/bader_neargrid.F90 :
subroutine bader_neargrid(parini)
    use mod_parini, only: typ_parini
    use mod_poisson_neargrid, only: typ_poisson
    type(typ_parini), intent(in):: parini
    integer,parameter::matl=1e5
end subroutine bader_neargrid
  subroutine ongrid_neargrid(poisson,d,i_dist)
    use mod_poisson_neargrid, only: typ_poisson
    type(typ_poisson), intent(inout):: poisson
    integer, intent(inout):: d(3)
    real(8),intent(in)::i_dist(-1:1,-1:1,-1:1)
end subroutine ongrid_neargrid
  function max_point(rho_cube)
    real(8),intent(in) :: rho_cube(3,3,3)
    logical :: max_point
end function max_point
  subroutine bounds(p1,p2,p3,p_max1,p_max2,p_max3)
    integer,intent(in)::p_max1,p_max2,p_max3
    integer,intent(inout)::p1,p2,p3
end subroutine bounds
  subroutine mat_vec(m,v,vp)
    real(8),intent(in),dimension(3,3) :: m
    real(8),intent(in),dimension(3) :: v
    real(8),intent(out),dimension(3) :: vp
end subroutine mat_vec
  subroutine vec_mat(v,m,vp)
    real(8),intent(in),dimension(3,3) :: m
    real(8),intent(in),dimension(3) :: v
    real(8),intent(out),dimension(3) :: vp
end subroutine vec_mat
  subroutine transposee(lattice,lat_car)
    real(8),intent(in),dimension(3,3) :: lattice
    real(8),intent(inout),dimension(3,3) :: lat_car
end subroutine transposee
  subroutine inverse(a,b)
    real(8),intent(in),dimension(3,3) :: a
    real(8),intent(out),dimension(3,3) :: b
end subroutine inverse
  function mat_vol(h)
    real(8),intent(in),dimension(3,3) :: h
    real(8) :: mat_vol
end function mat_vol
subroutine edag_refinement (last_iter,iter,poisson,i_dist,car_lat)
    use mod_poisson_neargrid, only: typ_poisson
    logical, intent(out):: last_iter
    integer, intent(in):: iter
    real(8), intent(in):: i_dist(-1:1,-1:1,-1:1),car_lat(3,3)
    type(typ_poisson), intent(inout):: poisson
end subroutine edag_refinement
  function is_edge_neargrid (poisson,d) result(is_edge)
    use mod_poisson_neargrid, only: typ_poisson
    type(typ_poisson):: poisson
    logical :: is_edge
    integer,intent(in) :: d(3)
end function is_edge_neargrid
  function m_point (poisson,d)
    use mod_poisson_neargrid, only: typ_poisson
    type(typ_poisson):: poisson
    logical :: m_point
    integer,intent(in) :: d(3)
end function m_point
subroutine cube_read_neargrid(filename,nat,rat,qat,poisson,vol)
    use mod_poisson_neargrid, only: typ_poisson
    character(*), intent(in):: filename
    real(8), intent(out):: rat(3,nat), qat(nat)
    real(8), intent(in):: vol
    type(typ_poisson), intent(out):: poisson
    integer:: nat, iat, igpx, igpy, igpz, ind, ios, iline, iatom
end subroutine cube_read_neargrid
subroutine m_neargrad (poisson,d,i_dist,car_lat)
    use mod_poisson_neargrid, only: typ_poisson
    type(typ_poisson), intent(inout):: poisson 
    integer, intent(inout):: d(3)
    real(8), intent(in):: i_dist(-1:1,-1:1,-1:1),car_lat(3,3)
end subroutine m_neargrad
subroutine path_reallocate(poisson, allocat2)
    use mod_poisson_neargrid, only: typ_poisson
    type(typ_poisson), intent(inout):: poisson 
    integer :: allocat2
end subroutine path_reallocate
subroutine vol_reallocate(poisson, allocat2)
    use mod_poisson_neargrid, only: typ_poisson
    type(typ_poisson), intent(inout):: poisson 
    integer :: allocat2
end subroutine vol_reallocate
subroutine near_grad (poisson,d,i_dist,car_lat)
    use mod_poisson_neargrid, only: typ_poisson
    type(typ_poisson), intent(inout):: poisson
    integer, intent(inout):: d(3)
    real(8), intent(in):: i_dist(-1:1,-1:1,-1:1),car_lat(3,3)
end subroutine near_grad
! ./src/bader_ongrid.F90 :
subroutine bader_ongrid(parini)
    use mod_parini, only: typ_parini
    use mod_poisson_ongrid, only: typ_poisson
    type(typ_parini), intent(in):: parini
end subroutine bader_ongrid
  subroutine ongrid_ongrid(i_dist,rho_cube,nx_now_grid,nx_cube_up,nx_cube_down,ny_now_grid,ny_cube_up,ny_cube_down,nz_now_grid,nz_cube_up,nz_cube_down,nx_now,ny_now,nz_now)
    real(8),intent(in)::i_dist(-1:1,-1:1,-1:1)
    integer,intent(in)::nx_cube_up,nx_cube_down,ny_cube_up,ny_cube_down,nz_cube_up,nz_cube_down
    real(8),intent(in)::rho_cube(3,3,3)
    INTEGER,INTENT(in) ::nx_now,ny_now,nz_now
    integer,intent(out)::nx_now_grid,ny_now_grid,nz_now_grid
end subroutine ongrid_ongrid
subroutine cube_read_ongrid(filename,nat,rat,qat,poisson)
    use mod_poisson_ongrid, only: typ_poisson
    character(*), intent(in):: filename
    real(8), intent(out):: rat(3,nat), qat(nat)
    type(typ_poisson), intent(out):: poisson
    integer:: nat, iat, igpx, igpy, igpz, ind, ios, iline, iatom
end subroutine cube_read_ongrid
! ./src/bader_weight.F90 :
subroutine bader_weight(parini)
    use mod_parini, only: typ_parini
    use mod_poisson_weight, only: typ_poisson
    type(typ_parini), intent(in):: parini
    integer,parameter::matl=1e5
end subroutine bader_weight
  subroutine calc_weight(poisson, p,nat)
    use mod_poisson_weight, only: typ_poisson
    type(typ_poisson), intent(inout):: poisson
    integer,intent(in) :: nat
    integer :: p(3), pn(3)
end subroutine calc_weight
 function is_neighbor(poisson, p, vol)
    use mod_poisson_weight, only: typ_poisson
    type(typ_poisson):: poisson
    logical :: is_neighbor
    integer,dimension(3),intent(in) :: p
    integer,intent(in) :: vol
end function is_neighbor
  subroutine ongrid_weight(poisson,d,i_dist)
    use mod_poisson_weight, only: typ_poisson
    type(typ_poisson), intent(inout):: poisson
    integer, intent(inout):: d(3)
    real(8),intent(in)::i_dist(-1:1,-1:1,-1:1)
end subroutine ongrid_weight
subroutine edag_refinement_weight (last_iter,iter,poisson,i_dist,car_lat)
    use mod_poisson_weight, only: typ_poisson
    logical, intent(out):: last_iter
    integer, intent(in):: iter
    real(8), intent(in):: i_dist(-1:1,-1:1,-1:1),car_lat(3,3)
    type(typ_poisson), intent(inout):: poisson
end subroutine edag_refinement_weight
  function is_edge_weight (poisson,d) result(is_edge)
    use mod_poisson_weight, only: typ_poisson
    type(typ_poisson):: poisson
    logical :: is_edge
    integer,intent(in) :: d(3)
end function is_edge_weight
  function m_point_weight (poisson,d) result(m_point)
    use mod_poisson_weight, only: typ_poisson
    type(typ_poisson):: poisson
    logical :: m_point
    integer,intent(in) :: d(3)
end function m_point_weight
subroutine cube_read_weight(filename,nat,rat,qat,poisson,vol)
    use mod_poisson_weight, only: typ_poisson
    character(*), intent(in):: filename
    real(8), intent(out):: rat(3,nat), qat(nat)
    real(8), intent(in):: vol
    type(typ_poisson), intent(out):: poisson
    integer:: nat, iat, igpx, igpy, igpz, ind, ios, iline, iatom
end subroutine cube_read_weight
subroutine m_neargrad_weight (poisson,d,i_dist,car_lat)
    use mod_poisson_weight, only: typ_poisson
    type(typ_poisson), intent(inout):: poisson 
    integer, intent(inout):: d(3)
    real(8), intent(in):: i_dist(-1:1,-1:1,-1:1),car_lat(3,3)
end subroutine m_neargrad_weight
subroutine path_reallocate_weight(poisson, allocat2)
    use mod_poisson_weight, only: typ_poisson
    type(typ_poisson), intent(inout):: poisson 
    integer :: allocat2
end subroutine path_reallocate_weight
subroutine vol_reallocate_weight(poisson, allocat2)
    use mod_poisson_weight, only: typ_poisson
    type(typ_poisson), intent(inout):: poisson 
    integer :: allocat2
end subroutine vol_reallocate_weight
subroutine near_grad_weight (poisson,d,i_dist,car_lat)
    use mod_poisson_weight, only: typ_poisson
    type(typ_poisson), intent(inout):: poisson
    integer, intent(inout):: d(3)
    real(8), intent(in):: i_dist(-1:1,-1:1,-1:1),car_lat(3,3)
end subroutine near_grad_weight
! ./src/basic.F90 :
subroutine elim_white_space(string)
    character(256), intent(inout):: string
end subroutine elim_white_space
function delta_kronecker(i,j) result(delta)
    integer, intent(in):: i, j
    real(8):: delta
end function delta_kronecker
subroutine backtocell_alborz(nat,latvec,rxyz_red)
    integer, intent(in):: nat
    real(8), intent(in):: latvec(3,3)
    real(8), intent(inout):: rxyz_red(3,nat)
end subroutine backtocell_alborz
subroutine getvol_alborz(cellvec,vol)
    real(8), intent(in):: cellvec(3,3)
    real(8), intent(out):: vol
end subroutine getvol_alborz
subroutine pbc_distance1_alborz(cellvec,xred_1,xred_2,distance2,dxyz)
    real(8), intent(in):: cellvec(3,3), xred_1(3), xred_2(3)
    real(8), intent(inout):: distance2, dxyz(3)
end subroutine pbc_distance1_alborz
subroutine n_rep_dim_alborz(cellvec,rcut,nec1,nec2,nec3)
    real(8), intent(in):: cellvec(3,3), rcut
    integer, intent(out):: nec1, nec2, nec3
end subroutine n_rep_dim_alborz
subroutine nveclatvec_alborz(cellvec,vn)
    real(8), intent(in) :: cellvec(3,3)
    real(8), intent(out):: vn(3,3)
end subroutine nveclatvec_alborz
subroutine dist2plane_alborz(r1,vn,r0,dist)
    real(8), intent(in):: r1(3), vn(3), r0(3)
    real(8), intent(out):: dist
end subroutine dist2plane_alborz
subroutine write_atomic_file_ascii_alborz(filename,nat,xred,latvec0,energy,pressure,printval1,printval2,kinds)
integer:: nat,natin,iat
character(40):: filename,units
real(8):: pos(3,nat),xred(3,nat),latvec(3,3),latvec0(3,3),dproj(6),rotmat(3,3),v(3,3),ucvol
real(8):: energy, etotal, enthalpy, enthalpy_at,pressure,printval1,printval2
integer:: Kinds(nat)
end subroutine write_atomic_file_ascii_alborz
subroutine dproj2latvec_alborz(dproj,cellvec)
    real(8), intent(in):: dproj(6)
    real(8), intent(out):: cellvec(3,3)
end subroutine dproj2latvec_alborz
subroutine latvec2dproj_alborz(dproj,latvec,rotmat,rxyz,nat)
    integer, intent(in):: nat
    real(8),intent(inout):: dproj(6), latvec(3,3), rotmat(3,3), rxyz(3,nat)
end subroutine latvec2dproj_alborz
subroutine cross_product_alborz(a,b,c)
    real(8), intent(in):: a(3), b(3)
    real(8), intent(out):: c(3)
end subroutine cross_product_alborz
subroutine rotation_alborz(angle,axe,rotmat)
    real(8), intent(in):: angle
    real(8), intent(in):: axe(3)
    real(8), intent(out):: rotmat(3,3)
end subroutine rotation_alborz
subroutine fxyz_cart2int_alborz(nat,v_cart,cv,v_int)
    integer, intent(in):: nat
    real(8), intent(in):: v_cart(3,nat), cv(3,3)
    real(8), intent(out):: v_int(3,nat)
end subroutine fxyz_cart2int_alborz
subroutine fxyz_red2cart(nat,fint,cv,fcart)
    integer, intent(in):: nat
    real(8), intent(in):: fint(3,nat), cv(3,3)
    real(8), intent(out):: fcart(3,nat)
end subroutine fxyz_red2cart
subroutine count_words(str,n)
    character(*), intent(in):: str
    integer, intent(out):: n
end subroutine count_words
subroutine count_substring(str1,str2,n)
    character(*), intent(in) :: str1, str2
    integer, intent(out):: n
end subroutine count_substring
! ./src/basic_minhocao.F90 :
 subroutine cross_product(a,b,crossp)
 real(8)::a(3),b(3)
 real(8)::crossp(3)
end subroutine cross_product
 subroutine dot_p(a,b,dotp)
 real(8)::a(3),b(3)
 real(8)::dotp(3)
end subroutine dot_p
subroutine bin_write(filename,array,n)
integer:: n
real(8):: array(n)
character(40):: filename
end subroutine bin_write
subroutine bin_read(filename,array,n)
integer:: n
real(8):: array(n)
character(40):: filename
end subroutine bin_read
function round(enerd,accur)
  real*8 enerd,accur,round
end function round
 subroutine rotation(rotmat,angle,axe)
 real(8),INTENT(IN) :: angle
 real(8),INTENT(IN) :: axe(3)
 real(8):: rotmat(3,3),cosang,sinang
end subroutine rotation
 subroutine invertmat(mat,matinv,n)
 real(8),intent(in) :: mat(n,n)
 integer               :: n
 real(8)               :: matinv(n,n),det(3),a(n,n),div
 integer               :: IPIV(n), INFO
end subroutine invertmat
subroutine hunt(xx,n,x,jlo)
  integer :: jlo,n
  real(kind=8) :: x,xx(n)
end subroutine hunt
! ./src/basic_utilities.F90 :
subroutine elim_moment_alborz(nat,atomic_vector)
  integer, intent(in):: nat
  real(8), intent(inout):: atomic_vector(3,nat)
end subroutine elim_moment_alborz
subroutine elim_moment_mass(nat,atomic_vector,atomic_mass)
  integer, intent(in):: nat
  real(8), intent(inout):: atomic_vector(3,nat)
  real(8), intent(in):: atomic_mass(nat)
end subroutine elim_moment_mass
subroutine calc_rotation_eigenvectors(nat,rat0,vrot)
  integer, intent(in) :: nat
  real(8), dimension(3*nat), intent(in) :: rat0
  real(8), dimension(3*nat,3), intent(out) :: vrot
  character(len=*), parameter :: subname='elim_torque_reza_alborz'
  real(8), dimension(3*nat) :: rat
end subroutine calc_rotation_eigenvectors
subroutine elim_torque_reza_alborz(nat,rat0,fat)
  integer, intent(in) :: nat
  real(8), dimension(3*nat), intent(in) :: rat0
  real(8), dimension(3*nat), intent(inout) :: fat
  character(len=*), parameter :: subname='elim_torque_reza_alborz'
  real(8), dimension(3*nat) :: rat
  real(8), dimension(3*nat,3) :: vrot
end subroutine elim_torque_reza_alborz
subroutine mycross(a,b,c)
  real(8), dimension(3), intent(in):: a,b
  real(8), dimension(3), intent(out):: c
end subroutine mycross
subroutine moment_of_inertia_alborz(nat,rat,teneria,evaleria)
  integer, intent(in) :: nat
  real(8), dimension(3,nat), intent(in) :: rat
  real(8), dimension(3), intent(out) :: evaleria
  real(8), dimension(3,3), intent(out) :: teneria
  character(len=*), parameter :: subname='moment_of_inertia_alborz'
  integer, parameter::lwork=100
end subroutine moment_of_inertia_alborz
subroutine normalizevector_alborz(n,v)
    integer, intent(in):: n
    real(8), intent(inout):: v(n)
end subroutine normalizevector_alborz
subroutine calnorm(n,v,vnrm)
    integer, intent(in):: n
    real(8), intent(in):: v(n)
    real(8), intent(out):: vnrm
end subroutine calnorm
subroutine calmaxforcecomponent(n,v,vmax)
    integer, intent(in):: n
    real(8), intent(in):: v(n)
    real(8), intent(out):: vmax
end subroutine calmaxforcecomponent
subroutine rxyz_cart2int_alborz(nat,latvec,rxyzcart,rxyzint)
    integer, intent(in):: nat
    real(8), intent(in):: rxyzcart(3,nat), latvec(3,3)
    real(8), intent(out):: rxyzint(3,nat)
end subroutine rxyz_cart2int_alborz
subroutine rxyz_int2cart_alborz(nat,cellvec,rat_int,rat_cart)
    integer, intent(in):: nat
    real(8), intent(in):: cellvec(3,3), rat_int(3,nat)
    real(8), intent(inout):: rat_cart(3,nat)
end subroutine rxyz_int2cart_alborz
subroutine invertmat_alborz(a,ainv)
    real(8),intent(in):: a(3,3)
    real(8),intent(out):: ainv(3,3)
end subroutine invertmat_alborz
subroutine invertmat_alborz_qp(a,ainv)
    real(16),intent(in):: a(3,3)
    real(16),intent(out):: ainv(3,3)
end subroutine invertmat_alborz_qp
subroutine convertupper(str)
    character(*), intent(inout):: str
end subroutine convertupper
subroutine convertlower(str)
    character(*), intent(inout) :: str
end subroutine convertlower
subroutine check_whether_time_exceeded
end subroutine check_whether_time_exceeded
subroutine expdist(n,x)
    integer, intent(in):: n
    real(8), intent(out):: x(n)
    real(8), parameter::eps=1.d-8
end subroutine expdist
subroutine gausdist_alborz(n,x)
    integer, intent(in) ::n
    real(8), intent(out) :: x(n)
    real(8), parameter:: eps=1.d-8
end subroutine gausdist_alborz
subroutine randdist(a,n,x)
    real(8), intent(in) :: a
    integer, intent(in):: n
    real(8), intent(out) :: x(n)
end subroutine randdist
subroutine hunt2(n,x,p,ip)
    integer, intent(in):: n
    real(8), intent(in):: x(n), p
    integer, intent(out):: ip
end subroutine hunt2
subroutine hpsort(n,ra)
    real*8 ::ra(n)
end subroutine hpsort
function flm_index(str1,str2) result(ind)
    character(*), intent(in):: str1, str2
    integer:: ind, indp, len_str
end function flm_index
subroutine elim_fixed_at(parini,nat,x)
use mod_parini, only: typ_parini
type(typ_parini), intent(in):: parini
integer:: iat,nat
real(8):: x(3,nat)
end subroutine elim_fixed_at
subroutine elim_fixed_lat(parini,latvec,x)
use mod_parini, only: typ_parini
type(typ_parini), intent(in):: parini
real(8):: x(3,3),latvec(3,3),lenlat,tmpvec(3)
end subroutine elim_fixed_lat
        subroutine elim_moment(nat,vxyz,atmass)
        real(8):: vxyz(3,nat),sx,sz,sy,atmass(nat)
        integer:: iat,nat       
end subroutine elim_moment
subroutine elim_torque_cell(latvec0,vlat)
real(8), intent(in)    :: latvec0(3,3)
real(8), intent(inout) :: vlat(3,3)
end subroutine elim_torque_cell
subroutine diagcomp(latvec,x)
real(8):: latvec(3,3),x(3,3),xnrm,latvect(3,3),latvecinv(3,3),sigma(3,3)
end subroutine diagcomp
 subroutine backtocell(nat,latvec,rxyz_red)
 integer:: nat,i,iat,j
 real(8) :: latvec(3,3), rxyz_red(3,nat), rxyz(3,nat), crossp(3),a(3),b(3), nvec(3,3), dist(6),eps,count
end subroutine backtocell
 subroutine backtocell_cart(nat,latvec,rxyz)
 integer:: nat,i,iat,j
 real(8) :: latvec(3,3), rxyz(3,nat), crossp(3),a(3),b(3), nvec(3,3), dist(6),eps,count
 real(8) :: v(3,3),vol,rxyz_red(3,nat)
end subroutine backtocell_cart
! ./src/buckingham.F90 :
subroutine set_buckingham(atoms,tosifumi)
    use mod_atoms, only: typ_atoms
    use mod_shortrange, only: typ_tosifumi
    type(typ_atoms), intent(inout):: atoms
    type(typ_tosifumi), intent(inout):: tosifumi
end subroutine set_buckingham
! ./src/cell_linkedlists.F90 :
subroutine linkedlists_init(parini,atoms,cell,linked_lists)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, atom_allocate_old
    use mod_electrostatics, only: typ_linked_lists
    type(typ_parini), intent(in):: parini 
    type(typ_atoms), intent(in):: atoms
    real(8), intent(out):: cell(3)
    type(typ_linked_lists), intent(inout):: linked_lists
end subroutine linkedlists_init
subroutine linkedlists_final(linked_lists)
    use mod_electrostatics, only: typ_linked_lists
    type(typ_linked_lists), intent(inout):: linked_lists
end subroutine linkedlists_final
subroutine prepprimelast(atoms,linked_lists)
    use mod_atoms, only: typ_atoms, get_rat, update_rat
    use mod_electrostatics, only: typ_linked_lists
    type(typ_atoms), intent(in):: atoms
    type(typ_linked_lists), intent(inout):: linked_lists
end subroutine prepprimelast
subroutine make_list_new(parini,atoms,linked_lists,cell)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, get_rat
    use mod_electrostatics, only: typ_linked_lists
    type(typ_parini), intent(in):: parini 
    type(typ_atoms), intent(in):: atoms
    type(typ_linked_lists), intent(inout):: linked_lists
    real(8), intent(in):: cell(3)
end subroutine make_list_new
subroutine determine_sclimitsphere(linked_lists)
    use mod_electrostatics, only: typ_linked_lists
    type(typ_linked_lists), intent(inout):: linked_lists
end subroutine determine_sclimitsphere
subroutine call_linkedlist(parini,atoms,dbl_count,linked_lists,pia_arr)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, type_pairs, update_ratp
    use mod_linked_lists, only: typ_linked_lists, typ_pia_arr
    type(typ_parini), intent(in):: parini 
    type(typ_atoms), intent(in):: atoms 
    logical, intent(in):: dbl_count
    type(typ_linked_lists), intent(inout):: linked_lists
    type(typ_pia_arr), intent(inout):: pia_arr
end subroutine call_linkedlist
! ./src/cell_niggli.F90 :
subroutine fixcell_niggli(nat,latvec,xred)
integer:: nat,iat
real(8):: latvec(3,3),xred(3,nat),latvec_out(3,3),eps,transmat(3,3),imat(3,3),epsvol,a(3,3),vol
end subroutine fixcell_niggli
subroutine niggli(latvec_in,latvec_out,transmat,eps)
real(8):: latvec_in(3,3),latvec_out(3,3),transmat(3,3),eps,tmpmat(3,3)
end subroutine niggli
subroutine  a1_action(nigmat,tmpmat,eps)
real(8):: nigmat(6),tmpmat(3,3),eps
end subroutine a1_action
subroutine  a2_action(nigmat,tmpmat,eps)
real(8):: nigmat(6),tmpmat(3,3),eps
end subroutine a2_action
subroutine  a3_action(nigmat,tmpmat,eps)
real(8):: nigmat(6),tmpmat(3,3),eps
end subroutine a3_action
subroutine  a4_action(nigmat,tmpmat,eps)
real(8):: nigmat(6),tmpmat(3,3),eps
end subroutine a4_action
subroutine  a5_action(nigmat,tmpmat,eps)
real(8):: nigmat(6),tmpmat(3,3),eps
end subroutine a5_action
subroutine  a6_action(nigmat,tmpmat,eps)
real(8):: nigmat(6),tmpmat(3,3),eps
end subroutine a6_action
subroutine  a7_action(nigmat,tmpmat,eps)
real(8):: nigmat(6),tmpmat(3,3),eps
end subroutine a7_action
subroutine  a8_action(nigmat,tmpmat,eps)
real(8):: nigmat(6),tmpmat(3,3),eps
end subroutine a8_action
subroutine  n1_action(nigmat,tmpmat,eps)
real(8):: nigmat(6),tmpmat(3,3),eps,tmp
end subroutine n1_action
subroutine  n2_action(nigmat,tmpmat,eps)
real(8):: nigmat(6),tmpmat(3,3),eps,tmp
end subroutine n2_action
subroutine  n3_true_action(nigmat,tmpmat,eps)
real(8):: nigmat(6),tmpmat(3,3),eps,tmp
end subroutine n3_true_action
subroutine  n3_false_action(nigmat,tmpmat,eps)
real(8):: nigmat(6),tmpmat(3,3),eps,tmp
end subroutine n3_false_action
function def_gt_0(nigmat,eps)
real(8):: nigmat(6),eps
logical:: flt,def_gt_0
end function def_gt_0
function fpositive(nigmat,eps)
real(8):: nigmat(6),eps
logical:: fpositive,flt
end function fpositive
function flt(a,b,eps) 
real(8)::a,b,eps
logical:: flt
end function flt
function fgt(a,b,eps) 
real(8)::a,b,eps
logical:: fgt,flt
end function fgt
function fle(a,b,eps) 
real(8)::a,b,eps
logical:: fle
end function fle
function fge(a,b,eps) 
real(8)::a,b,eps
logical:: fge
end function fge
function feq(a,b,eps) 
real(8)::a,b,eps
logical:: feq,flt
end function feq
subroutine init_nigmat(latvec,nigmat,l,m,n,eps,pi)
real(8):: dist_ang(6),nigmat(6),pi,convang,eps,latvec(3,3)
integer:: l,m,n
end subroutine init_nigmat
subroutine init_lmn(nigmat,l,m,n,eps)
real(8):: dist_ang(6),nigmat(6),eps
integer:: l,m,n
end subroutine init_lmn
subroutine latvec2dist_ang(dist_ang,latvec,pi)
real(8):: dist_ang(6),latvec(3,3),pi,convang
end subroutine latvec2dist_ang
function is_niggli_cell(nigmat,eps) 
real(8):: nigmat(6),eps
logical::is_niggli_cell,is_buerger_cell,feq,fgt
end function is_niggli_cell
function  meets_primary_conditions(nigmat,eps) 
real(8):: nigmat(6),eps
logical::meets_primary_conditions,fgt
end function meets_primary_conditions
function  meets_main_conditions(nigmat,eps) 
real(8):: nigmat(6),eps
logical::meets_main_conditions,meets_primary_conditions,flt
end function meets_main_conditions
function is_buerger_cell(nigmat,eps) 
real(8):: nigmat(6),eps
logical::is_buerger_cell,meets_main_conditions,feq,fgt
end function is_buerger_cell
function typer(nigmat,eps) 
real(8):: nigmat(6),eps
integer:: nzero,npositive,i,typer
end function typer
! ./src/cell_oganov.F90 :
 subroutine correct_latvec_oganov(latvec,pos_red,nat,iproc)
 real(8)              :: latvec(3,3),rxyz(3,nat),pos_red(3,nat)  
 integer              :: i,nat,counter,iproc
end subroutine correct_latvec_oganov
! ./src/compare_lammps.F90 :
subroutine compare_lammps(parini,parres)
use mod_parini, only: typ_parini
type(typ_parini), intent(in):: parini
type(typ_parini), intent(inout):: parres
real(8):: latvec(3,3),xred(3,parini%nat),xcart(3,parini%nat),f_lammps(3,parini%nat),f(3,parini%nat),e_lammps,e,tmp_r,tmp_i,tilts(6),latvec_in(3,3),strten(6),latvec_box(3,3)
end subroutine compare_lammps
! ./src/constants_mod.F90 :
! ./src/convcheck.F90 :
subroutine convcheck(parini,nat,latvec_in,fcart_in,strten_in,target_pressure_habohr,strfact,fmax,fmax_at,fmax_lat,tolmxf,iexit)
use mod_parini, only: typ_parini
type(typ_parini), intent(in):: parini
integer:: nat, iexit,iat,istr,i
real(8):: latvec_in(3,3),fcart_in(3,nat),strten_in(6),target_pressure_habohr,fmax,dstr(6)
real(8):: tolmxf,strtarget(6),strfact,fmax_at,fmax_lat
end subroutine convcheck
! ./src/correct_latvec.F90 :
subroutine correct_latvec(latvec,pos_red,nat,correctalg,iout)
integer:: correctalg,nat,iproc,iout
real(8):: latvec(3,3),pos_red(3,nat),latvec0(3,3),diff(9)
end subroutine correct_latvec
! ./src/dynamics_md_fixlat.F90 :
subroutine MD_fixlat(parini,parres,latvec_in,xred_in,fcart_in,strten_in,vel_in,etot_in,iprec,counter,folder)
 use mod_parini, only: typ_parini
    type(typ_parini), intent(in):: parini
    type(typ_parini), intent(inout):: parres
    integer:: iat,iprec,istep
    real(8):: latvec_in(3,3), xred_in(3,parini%nat),fcart_in(3,parini%nat),vel_in(3,parini%nat), strten_in(6), etot_in, counter
    real(8):: rxyz(3,parini%nat),fxyz(3,parini%nat),fxyz_old(3,parini%nat),vxyz(3,parini%nat),amass(parini%nat)
    character(40):: filename,folder
end subroutine md_fixlat
! ./src/dynamics_mod.F90 :
! ./src/electrostatics_mod.F90 :
! ./src/enthalpy.F90 :
subroutine get_enthalpy(latvec,energy,pressure,enthalpy)
real(8):: acell(3),v(3,3),ucvol,pressure,latvec(3,3),energy,enthalpy
end subroutine get_enthalpy
! ./src/envelope.F90 :
SUBROUTINE envelope(x, y, n, vertex, nvert, iwk)
INTEGER :: n, vertex(n), nvert, iwk(n)
REAL(8) :: x(n), y(n)
end subroutine envelope
! ./src/es_coulomb_p3d_bias.F90 :
subroutine bias_potener_forces(parini,poisson,atoms,epotplane)
    use mod_electrostatics, only: typ_poisson
    use mod_atoms, only: typ_atoms
    use mod_parini, only: typ_parini
    type(typ_poisson), intent(inout):: poisson
    type(typ_atoms), intent(inout):: atoms
    type(typ_parini), intent(in):: parini
    real(8):: epotlong, epotplane 
end subroutine bias_potener_forces
subroutine erfc_surface_zero(parini,atoms,poisson,nlayer)
    use mod_electrostatics, only: typ_poisson
    use mod_atoms, only: typ_atoms, update_ratp
    use mod_electrostatics, only: typ_linked_lists
    use mod_parini, only: typ_parini
    type(typ_parini), intent(in):: parini
    type(typ_poisson), intent(inout):: poisson
    type(typ_atoms), intent(inout):: atoms
    integer:: nimat !number of image atoms.
    integer:: nlayer, igpx,igpy,igpz,mx,my,mz, mlimnlayer
end subroutine erfc_surface_zero
subroutine sollaplaceq(poisson,hz,cell,vl,vu)
    use mod_electrostatics, only: typ_poisson
    type(typ_poisson), intent(inout):: poisson
    real(8):: vl, vu , zlmzu , sinhzlmzu, zlmzuinv
    real(8):: cell(3)
    real(8):: hz , vlmvu, vlzumvuzl 
end subroutine sollaplaceq
 subroutine calculate_force_ener_plane(atoms,poisson,epot,nbgpz)
    use mod_electrostatics, only: typ_poisson
    use mod_atoms, only: typ_atoms, update_ratp
    type(typ_poisson), intent(inout):: poisson
    type(typ_atoms), intent(inout):: atoms
    integer, intent(in):: nbgpz
    real(8):: x,y,z ,t,tl ,epot ,t1,t2
    real(8):: fatp(3,atoms%nat) 
end subroutine calculate_force_ener_plane
subroutine LG_weight(nlx,nly,nlz,hx,hy,hz,wx,wy,wz)
    integer:: nlx ,nly, nlz !number of point for Lagrange interpolation
    real(8):: hx ,hy ,hz , hxp , hyp, hzp
    real(8):: wx(nlx), wy(nly), wz(nlz) 
end subroutine lg_weight
subroutine LGW(n, w, h, x, LGx, DLGx, ix1, nbgp)
    integer:: ix, jx,kx, ix1,ixo ,n ,n2 ,nbgp
    real(8):: w(n), q(n),LGx(n),DLGx(n),h ,x ,diffx,x1,protot,t
end subroutine lgw
subroutine LGW4(n, w, h, x, LGx, DLGx, ix1, nbgp)
    integer:: ix, jx, ix1,ixo ,n ,n2 ,nbgp
    real(8):: w(n), q(n), qinv(n), LGx(n), DLGx(n), h ,x ,diffx,x1,protot
end subroutine lgw4
subroutine surface_charge(parini,poisson,pot_short,vl,vu)
    use mod_electrostatics, only: typ_poisson
    use mod_parini, only: typ_parini
    type(typ_parini), intent(in):: parini
    type(typ_poisson), intent(inout):: poisson
    real(8):: t, tt ,density(poisson%ngpx,poisson%ngpy,2),vl,vu
    real(8)::hgzinv,pi,pot_layerl,pot_layeru,pot_short(poisson%ngpx,poisson%ngpy,2,5)
    real(8), parameter, dimension(5) :: cf5 = [-25.d0/12.d0,4.d0,-3.d0,+4.d0/3.d0,-0.25d0]
end subroutine surface_charge
subroutine determine_limitsphere(poisson,mboundg,mboundgy,nbgpx,nbgpy,nbgpz)
    use mod_electrostatics, only: typ_poisson
    type(typ_poisson), intent(inout):: poisson
    integer:: ix, iy, iz, mboundg(2,-nbgpy:nbgpy,-nbgpz:nbgpz), mboundgy(2,-nbgpz:nbgpz)
    integer:: nbgpx, nbgpy, nbgpz
end subroutine determine_limitsphere
subroutine bias_field_potener_forces(parini,poisson,atoms,epotplane)
    use mod_electrostatics, only: typ_poisson
    use mod_atoms, only: typ_atoms, update_ratp
    use mod_parini, only: typ_parini
    type(typ_poisson), intent(inout):: poisson
    type(typ_atoms), intent(inout):: atoms
    type(typ_parini), intent(in):: parini
    real(8):: epotlong, epotplane !, epotshort
end subroutine bias_field_potener_forces
! ./src/es_coulomb_p3d_dielec.F90 :
subroutine dielec_potener_forces(parini,poisson,atoms,epot_dielec)
    use mod_electrostatics, only: typ_poisson
    use mod_atoms, only: typ_atoms
    use mod_parini, only: typ_parini
    type(typ_poisson), intent(inout):: poisson
    type(typ_atoms), intent(inout):: atoms
    type(typ_parini), intent(in):: parini
    real(8), intent(out):: epot_dielec
    REAL(8), PARAMETER :: PI = 3.14159265358979312
end subroutine dielec_potener_forces
subroutine sollaplaceq_dielctric(parini,poisson,hz,cell)
    use mod_electrostatics, only: typ_poisson
    use mod_parini, only: typ_parini
    type(typ_parini), intent(in):: parini
    type(typ_poisson), intent(inout):: poisson
    real(8):: cell(3)
    real(8):: hz , vlmvu, vlzumvuzl 
end subroutine sollaplaceq_dielctric
subroutine diff_pot_pp(parini,poisson,pot_short,vl,vu,nlayer)
    use mod_electrostatics, only: typ_poisson
    use mod_parini, only: typ_parini
    type(typ_parini), intent(in):: parini
    type(typ_poisson), intent(inout):: poisson
    integer::ix, iy, iz, npl,npu,nlayer
    real(8):: vl,vu
    real(8) :: pot_short(poisson%ngpx,poisson%ngpy,2,nlayer)
    real(8), parameter, dimension(7) :: cf7 = [-1.d0/60.d0,3.d0/20.d0,-3.d0/4.d0,0.d0,3.d0/4.d0,-3.d0/20.d0,1.d0/60.d0]
end subroutine diff_pot_pp
! ./src/es_coulomb_p3d.F90 :
subroutine calculate_forces_energy(parini,poisson,atoms)
    use mod_electrostatics, only: typ_poisson
    use mod_atoms, only: typ_atoms, get_rat
    use mod_parini, only: typ_parini
    type(typ_poisson), intent(inout):: poisson
    type(typ_atoms), intent(inout):: atoms
    type(typ_parini), intent(in):: parini
end subroutine calculate_forces_energy
! ./src/es_coulomb_spline.F90 :
subroutine build_shortrange_spline(shortrange,spline,rcut,a)
    use mod_shortrange, only: typ_shortrange
    use mod_spline, only: typ_spline
    type(typ_shortrange), intent(in):: shortrange
    type(typ_spline), intent(inout):: spline
    real(8), intent(in):: rcut !first and second cutoff for the function
    real(8), intent(in):: a !prefacor in exponent of exponential function
    real(16):: fdspq(0:3,0:spline%nsp)
    real(16):: fspq(0:4,0:spline%nsp-1)
    real(16):: fdspq_1(0:3,0:spline%nsp)
    real(16):: fspq_1(0:4,0:spline%nsp-1)
    real(16):: fdspq_2(0:3,0:spline%nsp)
    real(16):: fspq_2(0:4,0:spline%nsp-1)
    real(16):: fdspq_3(0:3,0:spline%nsp)
    real(16):: fspq_3(0:4,0:spline%nsp-1)
    real(16):: fdspq_4(0:3,0:spline%nsp)
    real(16):: fspq_4(0:4,0:spline%nsp-1)
end subroutine build_shortrange_spline
subroutine build_spline(cal_f_fd_fdd,rcutq,hspq,aq,nsp,fspq,fdspq)
    external:: cal_f_fd_fdd
    real(16), intent(in):: rcutq !first and second cutoff for the function
    real(16), intent(in):: hspq
    real(16), intent(in):: aq !prefacor in exponent of exponential function
    integer, intent(in):: nsp
    real(16), intent(out):: fdspq(0:3,0:nsp)
    real(16), intent(out):: fspq(0:4,0:nsp-1)
end subroutine build_spline
subroutine erf_over_r(r,a,hsp,func,funcder,funcsecder)
    real(16), intent(in):: r
    real(16), intent(in):: a
    real(16), intent(in):: hsp
    real(16), intent(out):: func
    real(16), intent(out):: funcder
    real(16), intent(out):: funcsecder
end subroutine erf_over_r
subroutine one_over_r6(r,a,hsp,func,funcder,funcsecder)
    real(16), intent(in):: r
    real(16), intent(in):: a
    real(16), intent(in):: hsp
    real(16), intent(out):: func
    real(16), intent(out):: funcder
    real(16), intent(out):: funcsecder
end subroutine one_over_r6
subroutine one_over_r8(r,a,hsp,func,funcder,funcsecder)
    real(16), intent(in):: r
    real(16), intent(in):: a
    real(16), intent(in):: hsp
    real(16), intent(out):: func
    real(16), intent(out):: funcder
    real(16), intent(out):: funcsecder
end subroutine one_over_r8
subroutine exp_ar(r,a,hsp,func,funcder,funcsecder)
    real(16), intent(in):: r
    real(16), intent(in):: a
    real(16), intent(in):: hsp
    real(16), intent(out):: func
    real(16), intent(out):: funcder
    real(16), intent(out):: funcsecder
end subroutine exp_ar
! ./src/es_hartree_bps.F90 :
subroutine get_psolver_bps(poisson,atoms,ehartree)
    use mod_atoms, only: typ_atoms
    use mod_electrostatics, only: typ_poisson
    type(typ_poisson),intent(inout):: poisson
    type(typ_atoms), intent(inout):: atoms
    real(8), intent(out):: ehartree
end subroutine get_psolver_bps
subroutine init_psolver_bps(parini,atoms,poisson)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_electrostatics, only: typ_poisson
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(in):: atoms
    type(typ_poisson), intent(inout):: poisson
end subroutine init_psolver_bps
subroutine fini_psolver_bps(poisson)
    use mod_electrostatics, only: typ_poisson
    type(typ_poisson), intent(inout):: poisson
end subroutine fini_psolver_bps
subroutine set_ngp_bps(parini,atoms,poisson_rough,poisson)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_electrostatics, only: typ_poisson
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(in):: atoms
    type(typ_poisson), intent(in):: poisson_rough
    type(typ_poisson), intent(inout):: poisson
end subroutine set_ngp_bps
! ./src/es_hartree_fourier.F90 :
subroutine get_psolver_fourier(parini,poisson,atoms,gausswidth,ehartree,g)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, update_ratp
    use mod_electrostatics, only: typ_poisson
    type(typ_parini), intent(in):: parini
    type(typ_poisson),intent(inout):: poisson
    type(typ_atoms), intent(inout):: atoms
    real(8), intent(in):: gausswidth(atoms%nat)
    real(8), intent(out):: ehartree, g(atoms%nat)
end subroutine get_psolver_fourier
subroutine get_psolver_fourier_various(iverbose,nat,rat,ratred,qat,cv,gwsq,ecut,ehartree,fat,eqd,stress,celldv)
    integer, intent(in):: iverbose, nat
    real(8), intent(in):: rat(3,nat), qat(nat)
    real(8), intent(in):: cv(3,3), gwsq(nat), ecut
    real(8), intent(out):: ratred(3,nat), fat(3,nat), eqd(nat), ehartree, stress(3,3), celldv(3,3)
end subroutine get_psolver_fourier_various
subroutine get_psolver_fourier_identical(iverbose,nat,rat,ratred,qat,cv,alphasq,ecut,ehartree,fat,eqd,stress,celldv)
    integer, intent(in):: iverbose, nat
    real(8), intent(in):: rat(3,nat), qat(nat)
    real(8), intent(in):: cv(3,3), alphasq, ecut
    real(8), intent(out):: ratred(3,nat), fat(3,nat), eqd(nat), ehartree, stress(3,3), celldv(3,3)
end subroutine get_psolver_fourier_identical
! ./src/es_hartree_main.F90 :
subroutine init_hartree(parini,atoms,poisson,gausswidth)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_electrostatics, only: typ_poisson
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(in):: atoms
    type(typ_poisson), intent(inout):: poisson
    real(8), intent(in):: gausswidth(atoms%nat)
end subroutine init_hartree
subroutine fini_hartree(parini,atoms,poisson)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_electrostatics, only: typ_poisson
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(in):: atoms
    type(typ_poisson), intent(inout):: poisson
end subroutine fini_hartree
subroutine init_hartree_bps(parini,atoms,poisson)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_electrostatics, only: typ_poisson
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(in):: atoms
    type(typ_poisson), intent(inout):: poisson
end subroutine init_hartree_bps
subroutine init_hartree_p3d(parini,atoms,poisson)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_electrostatics, only: typ_poisson
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(in):: atoms
    type(typ_poisson), intent(inout):: poisson
end subroutine init_hartree_p3d
subroutine put_charge_density(parini,poisson)
    use mod_parini, only: typ_parini
    use mod_electrostatics, only: typ_poisson
    type(typ_parini), intent(in):: parini
    type(typ_poisson),intent(inout):: poisson
end subroutine put_charge_density
subroutine get_psolver(parini,poisson,atoms,gausswidth,ehartree)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_electrostatics, only: typ_poisson
    type(typ_parini), intent(in):: parini
    type(typ_poisson),intent(inout):: poisson
    type(typ_atoms), intent(inout):: atoms
    real(8), intent(in):: gausswidth(atoms%nat)
    real(8), intent(out):: ehartree
end subroutine get_psolver
subroutine get_hartree_grad_rho(parini,poisson,atoms,ehartree)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_electrostatics, only: typ_poisson
    type(typ_parini), intent(in):: parini
    type(typ_poisson),intent(inout):: poisson
    type(typ_atoms), intent(inout):: atoms
    real(8), intent(out):: ehartree
end subroutine get_hartree_grad_rho
subroutine get_hartree_force(parini,poisson,atoms)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_electrostatics, only: typ_poisson
    type(typ_parini), intent(in):: parini
    type(typ_poisson),intent(inout):: poisson
    type(typ_atoms), intent(inout):: atoms
end subroutine get_hartree_force
subroutine get_hartree(parini,poisson,atoms,gausswidth,ehartree)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_electrostatics, only: typ_poisson
    type(typ_parini), intent(in):: parini
    type(typ_poisson),intent(inout):: poisson
    type(typ_atoms), intent(inout):: atoms
    real(8), intent(in):: gausswidth(atoms%nat)
    real(8), intent(out):: ehartree
end subroutine get_hartree
subroutine apply_external_field(parini,atoms,poisson,ehartree,g,flag)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, update_ratp
    use mod_electrostatics, only: typ_poisson
    type(typ_parini), intent(in):: parini
    type(typ_poisson),intent(inout):: poisson
    type(typ_atoms), intent(inout):: atoms
    real(8), intent(inout):: ehartree, g(atoms%nat)
    character(5)::flag
end subroutine apply_external_field
subroutine real_part(parini,atoms,gausswidth,alpha,epotreal,gg,stress)
    use mod_parini, only: typ_parini
    use mod_linked_lists, only: typ_linked_lists
    use mod_atoms, only: typ_atoms, update_ratp
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
    real(8), intent(in):: gausswidth(atoms%nat),alpha
    real(8):: gg(atoms%nat),rr
    real(8):: cell(3) , vol, stress(3,3)
    real(8)::epotreal, alphatwoinv, rbetainv, alphasq, betainv
end subroutine real_part
! ./src/es_hartree_p3d.F90 :
subroutine init_psolver_p3d(poisson)
    use mod_electrostatics, only: typ_poisson
    type(typ_poisson), intent(inout):: poisson
end subroutine init_psolver_p3d
subroutine fini_psolver_p3d(poisson)
    use mod_electrostatics, only: typ_poisson
    type(typ_poisson), intent(inout):: poisson
end subroutine fini_psolver_p3d
subroutine get_psolver_p3d(parini,poisson,cell,hx,hy,hz,epot)
    use mod_parini, only: typ_parini
    use mod_electrostatics, only: typ_poisson
    type(typ_parini), intent(in):: parini
    type(typ_poisson), intent(inout):: poisson
    real(8):: cell(3) !cell array contains size of the simulation box.
    real(8):: hx, hy, hz
    real(8):: epot
end subroutine get_psolver_p3d
subroutine solve_syslinequ_p3d(poisson,hz,cell)
    use mod_electrostatics, only: typ_poisson
    type(typ_poisson), intent(inout):: poisson
    real(8):: hz, cell(3)
    integer, parameter:: nem=8 
    real(8):: d(poisson%ngpz+2*8) !nem was replaced by 8 to be able to compile interface_mod.F90
    real(8):: e1(poisson%ngpz), e2(poisson%ngpz-1), c(poisson%ngpz)
end subroutine solve_syslinequ_p3d
subroutine fdcoeff(ngpz,e1,e2,g,gsq,hz,hzsq)
    integer::ngpz,ngpzm1
    real(8)::e1(ngpz) !Diagonal elements of the matrix
    real(8)::e2(ngpz-1) !Offdiagonal elements of the matrix
    real(8)::gsqhzsq,a,hz,hzsq,gsq,g,diagonal_fl,diagonal,offdiagonal
end subroutine fdcoeff
subroutine prepare00(ngpz,nem,f,c,hz)
    integer::ngpz,i,nem,n
    real(8)::f(ngpz+2*nem),c(ngpz),eta(6),hz 
end subroutine prepare00
subroutine prepare(ngpz,nem,f,c,gsq,hz,hzsq)
    integer::ngpz,i,nem,n
    real(8)::f(ngpz+2*nem),c(ngpz),eta(6),gsq,hz,hzsq,a
end subroutine prepare
subroutine prepcoeff(hz,eta,coefftot1,coefftoti,coefftotn)
    real(8)::hz,coeff(16,8),coefftot1(17),coefftoti(17),coefftotn(17),eta(6)
end subroutine prepcoeff
subroutine get_beta_grid(hzsq,ngpz,analc00,beta_grid)
    real(8), intent(in):: hzsq, analc00(ngpz)
    integer, intent(in):: ngpz
    real(8), intent(out):: beta_grid
end subroutine get_beta_grid
! ./src/fingerprint_atorb.F90 :
subroutine get_fp_malypso(nat,rxyz,rcov,latvec,r_cut_in,kinds,nkinds,fp_dim,nl,fp)
integer:: nl !Number of l components, here only even ones 
integer:: fp_dim !Number of AB interactions, doublecounting eliminated
integer:: nat,nkinds
integer:: nbond(fp_dim,nat),kinds(nat)
integer:: llist(nl)
real(8):: latvec(3,3),rxyz(3,nat),fp(nl,fp_dim,nat),fp_ri(2,nl,nat),rcov(nkinds)
real(8):: sigma,r_cut_in(fp_dim) !Cutoff for each AB interaction
real(8):: r_cut(fp_dim,fp_dim)
real(8), parameter :: pi=3.141592653589793238462643383279502884197d0
real(8):: min_bond(fp_dim,fp_dim),ylm_r,ylm_i
end subroutine get_fp_malypso
subroutine get_distance_malypso(fp1,fp2,fp_dim,nat,kinds,nl,dist)
integer:: fp_dim,nl,yll,i,j,mode,nat,iarr,i_dim,ii,jj,kinds(nat),k(nat),iii,jjj,nmat,imax,imin
real(8):: fp1(nl,fp_dim,nat),fp2(nl,fp_dim,nat),dist,a(nat,nat),summ,vec(nl),vec1(nl),vec2(nl),norm1,norm2
end subroutine get_distance_malypso
SUBROUTINE assndx(mode, a, n, m, k, sum)
INTEGER, INTENT(IN)   :: mode
REAL(8), INTENT(IN OUT)  :: a(:,:)
INTEGER, INTENT(IN)   :: n
INTEGER, INTENT(IN)   :: m
INTEGER, INTENT(OUT)  :: k(:)
REAL(8), INTENT(OUT)  :: sum
end subroutine assndx
! ./src/fingerprint.F90 :
subroutine init_fp(parini,fp_len,latvec)
use mod_parini, only: typ_parini
type(typ_parini), intent(in):: parini
integer:: fp_len,iat,nmax
real(8):: convert,latvec(3,3),vol
end subroutine init_fp
subroutine get_fp(parini,fp_len,pos_red,latvec,fp)
use mod_parini, only: typ_parini
use fingerprint, only: fp_15_fp_size, fp_method, fp_11_rcut, fp_11_sigma, fp_11_dbin
type(typ_parini), intent(in):: parini
integer:: fp_len,iat,natmol
real(8):: fp(fp_len),pos_red(3,parini%nat),latvec(3,3),rxyz(3,parini%nat),vol,rcov_arr(parini%nat),fp_coganov_atomic(3,fp_15_fp_size,parini%ntypat_global,parini%nat)
real(8):: rvan(parini%nat) !nat*molecules)
character(len=2):: finalchar(parini%nat) ! dimension(nat*molecules)
end subroutine get_fp
! ./src/fingerprint_GOM.F90 :
subroutine get_fp_gauss(nat, ntypat, natx_sphere, typat, lseg, width_cutoff, nex_cutoff, alat, rxyz, rcov, fp)
  integer, intent(in) :: nat, ntypat, natx_sphere, lseg
  integer, dimension(nat), intent(in) :: typat
  real(8), intent(in) :: width_cutoff, nex_cutoff
  real(8), dimension(3,3), intent(in) :: alat
  real(8), dimension(3,nat), intent(in) :: rxyz
  real(8), dimension(nat), intent(in) :: rcov
  real(8), dimension(lseg*(ntypat+1), nat), intent(out) :: fp
  integer, parameter :: nwork = 100
  integer, dimension(lseg*natx_sphere) :: ind_small
  real(8), dimension(natx_sphere) :: amplitude
  real(8), dimension(lseg*natx_sphere) :: fpp
  real(8), dimension(3, natx_sphere) :: rxyz_sphere
  real(8), dimension(natx_sphere) :: rcov_sphere
  real(8), dimension(lseg*(ntypat+1), lseg*(ntypat+1)) :: omsa, omsb, omsaa, omsbb
end subroutine get_fp_gauss
subroutine get_distance_gauss(fp1, fp2, lseg, nat, ntypat, typat, fpd)
  integer, intent(in)  :: lseg, nat, ntypat
  integer, dimension(nat) :: typat
  real(8), dimension(lseg*(ntypat+1), nat) :: fp1, fp2
  real(8), intent(out) :: fpd
  real(8), dimension(nat, nat) :: cost
  real(8), dimension(ntypat)   :: dist
end subroutine get_distance_gauss
subroutine mltampl_4(nat,amplitude,om)
end subroutine mltampl_4
subroutine mltampl_1(nat,amplitude,om)
end subroutine mltampl_1
subroutine create_om_1(nat,rxyz,rcov,om)
end subroutine create_om_1
subroutine create_om_4(nat,rxyz,rcov,om)
end subroutine create_om_4
! ./src/fingerprint_MOLGOM.F90 :
subroutine get_distance_molgom(fp1,fp2,dist,lseg,molecules,molecules_sphere,principleev)
integer:: lseg,molecules,molecules_sphere,principleev
integer:: iassign(molecules)
real(8):: tt,dtt,dist
real(8):: cost(molecules,molecules)
real(8):: fp1(lseg*molecules_sphere*principleev,molecules),fp2(lseg*molecules_sphere*principleev,molecules)
end subroutine get_distance_molgom
subroutine create_contracted_om_1(width_overlap,principleev,nat,molecules,rxyz,rvan,amplitude,fp_t,lseg,write_files)
  logical :: write_files
  integer :: principleev, xyz, alpha, beta, mu, nu, k, i, j, m, n
  integer :: lwork, nat, molecules, lseg, info
  real*8 :: width_overlap
  real*8, dimension(3, molecules*nat):: rxyz
  real*8, dimension(molecules) :: amplitude
  real*8, dimension(3,nat) :: rxyz_temp
  real*8, dimension(molecules*nat) :: rvan
  real*8, dimension(nat) :: rvan_temp
  real*8, dimension(nat,principleev,molecules) :: em
  real*8, dimension(nat,nat) :: om
  real*8, dimension(nat*molecules,nat*molecules) :: om_b
  real*8, dimension(molecules,principleev,molecules,principleev) :: om_t
  real*8, dimension(nat) :: fp
  real*8, dimension(principleev*molecules) :: fp_t
end subroutine create_contracted_om_1
subroutine create_molom_1(nat,rxyz,rvan,om,width_overlap)
   real*8 :: width_overlap
end subroutine create_molom_1
subroutine periodic_fingerprint(parini,rxyz,alat0,finalchar,rvan,fpsall,nat)
   use mod_parini, only: typ_parini
   type(typ_parini), intent(in):: parini
   integer:: nat
   integer, dimension (parini%fp_18_expaparameter+1) :: shifting
   real*8, dimension (3,(parini%fp_18_expaparameter+1)**3) :: possibilites
   real*8, dimension(3,nat*parini%fp_18_molecules):: rxyz
   real*8, dimension(nat*parini%fp_18_molecules) :: rvan
   real*8, dimension(parini%fp_18_lseg*parini%fp_18_molecules_sphere*parini%fp_18_principleev,parini%fp_18_molecules) :: fpsall
   real*8, dimension (3,3)::alat,alat0
   character(len=2), dimension(nat*parini%fp_18_molecules) :: finalchar
   logical, dimension(parini%fp_18_molecules,(parini%fp_18_expaparameter+1)**3) :: is_copied
end subroutine periodic_fingerprint
subroutine findmolecule(parini,rxyz,alat0,finalchar,xred,char_type,typat,ntypat,nat)
    use mod_parini, only: typ_parini
    type(typ_parini), intent(in):: parini
   integer, dimension(nat*parini%fp_18_molecules), intent(in) :: typat
   character(2), dimension(ntypat), intent(in) :: char_type
   integer, intent(in):: ntypat,nat
   real*8, dimension(3,3), intent(in) :: alat0
   real*8, dimension(3,nat*parini%fp_18_molecules), intent(out) :: rxyz
   real*8, dimension(3,nat*parini%fp_18_molecules), intent(in)  :: xred
   character(len=2), dimension(nat*parini%fp_18_molecules), intent(out) :: finalchar
   real*8, dimension(3, nat*parini%fp_18_molecules,27) :: rxyz_b
   real*8, dimension(nat*parini%fp_18_molecules) :: rcov
   logical, dimension (2,nat*parini%fp_18_molecules,27) :: is_true
   logical, dimension (2,nat*parini%fp_18_molecules) :: first_search_is_true
   integer, dimension (nat*parini%fp_18_molecules) :: first_search_molecule_number
   integer, dimension (nat*parini%fp_18_molecules,27) :: molecule_number
   integer, dimension (nat*parini%fp_18_molecules,2) :: atom_number
   integer, dimension(parini%fp_18_molecules,2) :: finalmolecules
   integer, dimension (nat*parini%fp_18_molecules,27) :: finalmolecule_number
   character(len=2), dimension(nat*parini%fp_18_molecules,27) :: is_char
end subroutine findmolecule
subroutine sym2rcov(sym,rcov)
  real(8)  :: rcov
  character(len=2) :: sym  ! chemical symbol 
end subroutine sym2rcov
subroutine sym2rvan(sym,rvan)
  real(8)  :: rvan
  character(len=2) :: sym  ! chemical symbol 
end subroutine sym2rvan
! ./src/fit_elecpot.F90 :
subroutine subtask_fit_elecpot(parini)
    use mod_parini, only: typ_parini
    type(typ_parini), intent(in):: parini
end subroutine subtask_fit_elecpot
subroutine fit_elecpot(parini)
    use mod_parini, only: typ_parini
    use mod_electrostatics, only: typ_poisson
    use mod_atoms, only: typ_atoms, typ_atoms_arr, get_rat, set_rat, update_ratp
    use mod_ann, only: typ_cent, typ_ann_arr
    type(typ_parini), intent(in):: parini
end subroutine fit_elecpot
subroutine put_pot_sym_rzx(rat,hgx,hgy,hgz,nat,qat,gw,ng,lcn,reset,weight,dft_pot,cent_pot,qpar,apar,rpar)
    logical :: reset
    integer , intent(in):: nat, ng(1:3), lcn
    real(8) , intent(in):: rat(1:3,1:nat), hgx, hgy, hgz, qat(1:lcn,1:nat),gw(1:lcn,1:nat),weight(1:ng(1),1:ng(2),1:ng(3))
    real(8) , intent(in):: dft_pot(1:ng(1),1:ng(2),1:ng(3))
    real(8) , intent(out):: cent_pot(1:ng(1),1:ng(2),1:ng(3)), apar(1:lcn,1:nat), qpar(1:lcn,1:nat), rpar(1:3,1:nat)
    real(8) :: cent_pot_a_par(1:ng(1),1:ng(2),1:ng(3)), cent_pot_q_par(1:ng(1),1:ng(2),1:ng(3))
    real(8) :: cent_pot_x_par(1:ng(1),1:ng(2),1:ng(3)), cent_pot_y_par(1:ng(1),1:ng(2),1:ng(3)), cent_pot_z_par(1:ng(1),1:ng(2),1:ng(3))
end subroutine put_pot_sym_rzx
subroutine stdval_rzx(f,f_len,mean,std,var)
    integer, intent(in) :: f_len
    real(8), intent(in) :: f(f_len)
    real(8), intent(out) :: mean, std, var
    real(8) :: g(f_len)
end subroutine stdval_rzx
! ./src/flame_as_potential_mod.F90 :
! ./src/flame.F90 :
! ./src/flame_init_fini.F90 :
subroutine alborz_init(parini,parres,file_ini)
    use mod_task, only: typ_file_ini, time_start
    use mod_parini, only: typ_parini
    type(typ_file_ini), intent(inout):: file_ini
    type(typ_parini), intent(inout):: parini
    type(typ_parini), intent(inout):: parres
    character(len=*), parameter:: filename='flame_log.yaml'
end subroutine alborz_init
subroutine alborz_initialize_timing_categories
    character(len=*), parameter :: pscpt1='Alborz init-fin'
    character(len=*), parameter :: pscpt2='potential'
end subroutine alborz_initialize_timing_categories
subroutine alborz_final(parini,file_ini)
    use mod_parini, only: typ_parini
    use mod_task, only: typ_file_ini, time_start, time_end
    type(typ_parini), intent(inout):: parini
    type(typ_file_ini), intent(inout):: file_ini
end subroutine alborz_final
subroutine init_random_seed(parini)
    use mod_parini, only: typ_parini
    type(typ_parini), intent(in):: parini
end subroutine init_random_seed
subroutine set_atomc_types_info(parini)
    use mod_parini, only: typ_parini
    type(typ_parini), intent(inout):: parini
end subroutine set_atomc_types_info
subroutine flm_print_logo(parini)
    use mod_parini, only: typ_parini
    type(typ_parini), intent(inout):: parini
end subroutine flm_print_logo
! ./src/forcefield.F90 :
subroutine forcefield_init(parini,atoms)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, set_typat, set_qat
    use mod_electrostatics, only: typ_poisson
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
end subroutine forcefield_init
subroutine calculate_forces_energy_ff(parini,atoms)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
end subroutine calculate_forces_energy_ff
subroutine forcefield_final(parini,atoms)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(in):: atoms
end subroutine forcefield_final
! ./src/fp_distance.F90 :
subroutine get_fp_distance(parini,fp_len,fp1,fp2,fp_dist)
use mod_parini, only: typ_parini
type(typ_parini), intent(in):: parini
integer:: fp_len
real(8):: fp(fp_len),pos_red(3,parini%nat),latvec(3,3),rxyz(3,parini%nat),fp1(fp_len),fp2(fp_len),fp_dist
end subroutine get_fp_distance
! ./src/fpos_flat.F90 :
subroutine fpos_flat(parini,pressure,fpos,flat,strten,fcart,latvec,md_type) 
use mod_parini, only: typ_parini
type(typ_parini), intent(in):: parini
integer:: iat,i,j,md_type
real(8),dimension(3,parini%nat):: fcart,fpos
real(8),dimension(3,3)  :: latvec,tmplat,pressure,a,velmat,sigma,flat,str_matrix
real(8):: amass(parini%nat),latmass,crossp(3),strten(6),vol,vpostmp(3),volvel,trace3,vol_1_3
end subroutine fpos_flat
! ./src/fragments.F90 :
subroutine fragments(parini,latvec,xred,nfrag,xcart,fragarr,fragsize)
use mod_parini, only: typ_parini
type(typ_parini), intent(in):: parini
real(8),dimension(3,parini%nat), INTENT(IN) :: xred
real(8):: latvec(3,3),rotmat(3,3),dproj(6)
integer :: nfrag, nfragold
real(8):: ekin,vcm1,vcm2,vcm3,ekin0,scale,xcart(3,parini%nat)
integer, dimension(parini%nat):: fragarr,fragsize
end subroutine fragments
subroutine get_fragsize(fragsize,lhead,llist,nat,nmol)
integer:: nat,nmol,iat,ifrag,llist(nat),lhead(nmol),fragsize(nmol)
end subroutine get_fragsize
subroutine refragment(fragarr,nat)
integer:: fragarr(nat),nat,iat,jat,cnt,find,fragarr_tmp(nat)
end subroutine refragment
subroutine make_linked_list(fragarr,fragsize,lhead,llist,nat,nmol)
integer:: fragarr(nat),nat,iat,nmol,ifrag
integer:: lhead(nmol),llist(nat),fragsize(nmol)
end subroutine make_linked_list
subroutine get_cmass(cmass,masstot,xcart,amass,lhead,llist,nat,nmol)
integer:: nat,nmol,iat,ifrag,lhead(nmol),llist(nat)
real(8):: xcart(3,nat),amass(nat),masstot(nmol),cmass(3,nmol)
end subroutine get_cmass
subroutine get_inertia_tensor(parini,intens,inprin,inaxis,cmass,xcart,amass,lhead,llist,nat,nmol)
use mod_parini, only: typ_parini
type(typ_parini), intent(in):: parini
integer:: nat,nmol,iat,ifrag,i,j,llist(nat),lhead(nmol),LWORK,info
real(8):: xcart(3,nat),amass(nat),cmass(3,nmol),intens(3,3,nmol),dist2,xtmp(3)
real(8):: inprin(3,nmol),inaxis(3,3,nmol),diag_inert(3,3),tmp_vec(3),tmp_val
end subroutine get_inertia_tensor
subroutine get_fcm_torque(fcm,torque,fcart,quat,xcart_mol,lhead,llist,nat,nmol)
integer:: nat,nmol,iat,ifrag,llist(nat),lhead(nmol)
real(8),intent(in):: fcart(3,nat),quat(4,nmol),xcart_mol(3,nat)
real(8):: fcm(3,nmol),torque(3,nmol),crossp(3),xtmp(3),rotmat(3,3)
end subroutine get_fcm_torque
subroutine init_cm_mol(parini,latvec,xred,xcart_shifted,xred_cm,quat,amass,masstot,intens,inprin,inaxis,lhead,llist,nat,nmol)
use mod_parini, only: typ_parini
type(typ_parini), intent(in):: parini
real(8),intent(in):: latvec(3,3),xred(3,nat)
integer:: nat,nmol,iat,llist(nat),lhead(nmol),fragsize(nmol),imol,jmol,kmol
real(8):: xcart_in(3,nat),xcart_shifted(3,nat),cmass(3,nmol),amass(nat)
real(8):: masstot(nmol),angbohr,quat(4,nmol),xred_cm(3,nmol),xcart_tmp(3,nat)
real(8):: circular(3,3),tol,rot_c(3,3),rot_all(3,3),inprin(3,nmol),intens(3,3,nmol)
real(8):: inaxis(3,3,nmol),ident(3,3),tmp(3,nmol),quat_tmp(4),tmp_real(4),tmp_mat(3,3)
logical:: symtop(nmol),tmp_logical
end subroutine init_cm_mol
! ./src/gaussdist.F90 :
      subroutine gausdist(nat,vxyz,amass)
      real(8):: t1,t2,tt,amass(nat)
      real(8),parameter:: eps=1.d-8
      real(8),dimension(3*nat)::  vxyz
      integer:: nat,i
end subroutine gausdist
      subroutine gausdist_cell(latvec,vlat)
      real(8),parameter:: eps=1.d-8
      real(8)::  vlat(9),latvec(9)
end subroutine gausdist_cell
! ./src/genconf_diatomic.F90 :
subroutine genconf_diatomic(parini,genconf)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms_all, typ_file_info, atom_all_allocate, atom_all_deallocate
    use mod_genconf, only: typ_genconf
    type(typ_parini), intent(in):: parini
    type(typ_genconf), intent(in):: genconf
end subroutine genconf_diatomic
! ./src/genconf_mod.F90 :
! ./src/genconf_random.F90 :
subroutine genrandom(parini,genconf)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, typ_file_info, atom_allocate_old, atom_deallocate_old
    use mod_genconf, only: typ_genconf
    type(typ_parini), intent(inout):: parini
    type(typ_genconf), intent(in):: genconf
    real(8),parameter::pi=4.d0*atan(1.d0)
end subroutine genrandom
! ./src/genconf_rangrow.F90 :
subroutine rangrow(parini,genconf)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, typ_file_info, atom_allocate_old, atom_deallocate_old
    use mod_genconf, only: typ_genconf
    type(typ_parini), intent(inout):: parini
    type(typ_genconf), intent(in):: genconf
end subroutine rangrow
! ./src/genconf_trimer.F90 :
subroutine genconf_trimer(parini,genconf)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms_all, typ_file_info, atom_all_allocate, atom_all_deallocate
    use mod_genconf, only: typ_genconf
    type(typ_parini), intent(in):: parini
    type(typ_genconf), intent(in):: genconf
end subroutine genconf_trimer
! ./src/grid_basic.F90 :
subroutine get_glimitsphere(hx,hy,hz,nbgpx,nbgpy,nbgpz,mboundg)
    real(8), intent(in):: hx, hy, hz
    integer, intent(in):: nbgpx, nbgpy, nbgpz
    integer, intent(out):: mboundg(1:2,-nbgpy:nbgpy,-nbgpz:nbgpz)
end subroutine get_glimitsphere
subroutine init_grid_param(nat,rxyz,cv,rgcut,ngx,ngy,ngz,ratred,vol,nlimsq,nagx,nagy,nagz,nbgx,nbgy,nbgz)
    integer, intent(in):: nat
    real(8), intent(in):: rxyz(3,nat)
    real(8), intent(in):: cv(3,3)
    real(8), intent(in):: rgcut
    integer, intent(in):: ngx, ngy, ngz
    real(8), intent(out):: ratred(3,nat), vol
    integer, intent(out):: nlimsq, nagx, nagy, nagz, nbgx, nbgy, nbgz
end subroutine init_grid_param
subroutine charge_back_to_cell(ngx,ngy,ngz,nagx,nagy,nagz,ibcx,wa,rho)
    integer, intent(in):: ngx, ngy, ngz, nagx, nagy, nagz, ibcx
    real(8), intent(in):: wa(1-nagx:ngx+nagx,1-nagy:ngy+nagy,1-nagz:ngz+nagz)
    real(8), intent(inout):: rho(ngx,ngy,ngz)
end subroutine charge_back_to_cell
subroutine potential_on_extended_grid(lda,ngx,ngy,ngz,nagx,nagy,nagz,ibcx,pot,wa)
    integer, intent(in):: lda, ngx, ngy, ngz, nagx, nagy, nagz, ibcx
    real(8), intent(in):: pot(lda,ngy,ngz)
    real(8), intent(out):: wa(1-nagx:ngx+nagx,1-nagy:ngy+nagy,1-nagz:ngz+nagz)
end subroutine potential_on_extended_grid
! ./src/grid_gto_sym.F90 :
subroutine put_gto_sym(parini,bc,reset,nat,rxyz,qat,gw,rgcut,ngx,ngy,ngz,hgrid,rho)
    use mod_parini, only: typ_parini
    type(typ_parini), intent(in):: parini
    character(*), intent(in):: bc
    logical, intent(in):: reset
    integer, intent(in):: nat
    real(8), intent(in):: rxyz(3,nat)
    real(8), intent(in):: qat(nat)
    real(8), intent(in):: gw(nat)
    real(8), intent(in):: rgcut
    integer, intent(in):: ngx, ngy, ngz
    real(8), intent(in):: hgrid(3,3)
    real(8), intent(inout):: rho(ngx,ngy,ngz)
end subroutine put_gto_sym
subroutine rqgrad_gto_sym(parini,bc,nat,rxyz,qat,gw,rgcut,lda,ngx,ngy,ngz,hgrid,pot,rgrad,qgrad)
    use mod_parini, only: typ_parini
    type(typ_parini), intent(in):: parini
    character(*), intent(in):: bc
    integer, intent(in):: nat
    real(8), intent(in):: rxyz(3,nat)
    real(8), intent(in):: qat(nat)
    real(8), intent(in):: gw(nat)
    real(8), intent(in):: rgcut
    integer, intent(in):: lda, ngx, ngy, ngz
    real(8), intent(in):: hgrid(3,3)
    real(8), intent(in):: pot(lda,ngy,ngz)
    real(8), intent(out):: rgrad(3,nat), qgrad(nat)
end subroutine rqgrad_gto_sym
subroutine force_gto_sym(parini,bc,nat,rxyz,qat,gw,rgcut,lda,ngx,ngy,ngz,hgrid,pot,fat)
    use mod_parini, only: typ_parini
    type(typ_parini), intent(in):: parini
    character(*), intent(in):: bc
    integer, intent(in):: nat
    real(8), intent(in):: rxyz(3,nat)
    real(8), intent(in):: qat(nat)
    real(8), intent(in):: gw(nat)
    real(8), intent(in):: rgcut
    integer, intent(in):: lda, ngx, ngy, ngz
    real(8), intent(in):: hgrid(3,3)
    real(8), intent(in):: pot(lda,ngy,ngz)
    real(8), intent(out):: fat(3,nat)
end subroutine force_gto_sym
subroutine gwrqgrad_gto_sym(parini,bc,nat,rxyz,qat,gw,rgcut,lda,ngx,ngy,ngz,hgrid,pot,rgrad,qgrad,agrad)
    use mod_parini, only: typ_parini
    type(typ_parini), intent(in):: parini
    character(*), intent(in):: bc
    integer, intent(in):: nat
    real(8), intent(in):: rxyz(3,nat)
    real(8), intent(in):: qat(nat) 
    real(8), intent(in):: gw(nat)
    real(8), intent(in):: rgcut
    integer, intent(in):: lda, ngx, ngy, ngz
    real(8), intent(in):: hgrid(3,3)
    real(8), intent(in):: pot(lda,ngy,ngz)
    real(8), intent(out):: rgrad(3,nat), qgrad(nat), agrad(nat)
end subroutine gwrqgrad_gto_sym
subroutine rhograd_gto_sym(parini,bc,reset,nat,rxyz,cv,qat,gw,rgcut,ngx,ngy,ngz,rho,rho_q_par,rho_a_par)
    use mod_parini, only: typ_parini
    type(typ_parini), intent(in):: parini
    character(*), intent(in):: bc
    logical, intent(in):: reset
    integer, intent(in):: nat
    real(8), intent(in):: rxyz(3,nat)
    real(8), intent(in):: cv(3,3)
    real(8), intent(in):: qat(nat)
    real(8), intent(in):: gw(nat)
    real(8), intent(in):: rgcut
    integer, intent(in):: ngx, ngy, ngz
    real(8), intent(inout):: rho(ngx,ngy,ngz),rho_a_par(ngx,ngy,ngz),rho_q_par(ngx,ngy,ngz)
end subroutine rhograd_gto_sym
! ./src/grid_gto_sym_ortho.F90 :
subroutine put_gto_sym_ortho(parini,bc,reset,nat,rxyz,qat,gw,rgcut,ngx,ngy,ngz,hgrid,rho)
    use mod_parini, only: typ_parini
    type(typ_parini), intent(in):: parini
    character(*), intent(in):: bc
    logical, intent(in):: reset
    integer, intent(in):: nat
    real(8), intent(in):: rxyz(3,nat)
    real(8), intent(in):: qat(nat)
    real(8), intent(in):: gw(nat)
    real(8), intent(in):: rgcut
    integer, intent(in):: ngx, ngy, ngz
    real(8), intent(in):: hgrid(3,3)
    real(8), intent(inout):: rho(ngx,ngy,ngz)
end subroutine put_gto_sym_ortho
subroutine qgrad_gto_sym_ortho(parini,bc,nat,rxyz,qat,gw,rgcut,lda,ngx,ngy,ngz,hgrid,pot,g)
    use mod_parini, only: typ_parini
    type(typ_parini), intent(in):: parini
    character(*), intent(in):: bc
    integer, intent(in):: nat
    real(8), intent(in):: rxyz(3,nat)
    real(8), intent(in):: qat(nat)
    real(8), intent(in):: gw(nat)
    real(8), intent(in):: rgcut
    integer, intent(in):: lda, ngx, ngy, ngz
    real(8), intent(in):: hgrid(3,3)
    real(8), intent(in):: pot(lda,ngy,ngz)
    real(8), intent(inout):: g(nat)
end subroutine qgrad_gto_sym_ortho
subroutine force_gto_sym_ortho(parini,bc,nat,rxyz,qat,gw,rgcut,lda,ngx,ngy,ngz,hgrid,pot,fat)
    use mod_parini, only: typ_parini
    type(typ_parini), intent(in):: parini
    character(*), intent(in):: bc
    integer, intent(in):: nat
    real(8), intent(in):: rxyz(3,nat)
    real(8), intent(in):: qat(nat)
    real(8), intent(in):: gw(nat)
    real(8), intent(in):: rgcut
    integer, intent(in):: lda, ngx, ngy, ngz
    real(8), intent(in):: hgrid(3,3)
    real(8), intent(inout):: pot(lda,ngy,ngz)
    real(8), intent(out):: fat(3,nat)
end subroutine force_gto_sym_ortho
! ./src/grid_rp4gto_sym.F90 :
subroutine put_rp4gto_sym(parini,bc,reset,nat,rxyz,cv,qat,gw,rgcut,ngx,ngy,ngz,rho,rho_q_par,rho_a_par)
    use mod_parini, only: typ_parini
    type(typ_parini), intent(in):: parini
    character(*), intent(in):: bc
    logical, intent(in):: reset
    integer, intent(in):: nat
    real(8), intent(in):: rxyz(3,nat)
    real(8), intent(in):: cv(3,3)
    real(8), intent(in):: qat(nat)
    real(8), intent(in):: gw(nat)
    real(8), intent(in):: rgcut
    integer, intent(in):: ngx, ngy, ngz
    real(8), intent(inout):: rho(ngx,ngy,ngz),rho_a_par(ngx,ngy,ngz),rho_q_par(ngx,ngy,ngz)
end subroutine put_rp4gto_sym
! ./src/hung.F90 :
subroutine hung(N,A,F,Z)
      integer:: n
      real(8)::  A(n,n),Z,U(n),V(n)
      integer F(N),FB(n), RC(n)
end subroutine hung
subroutine INCR_inalborz(n,F,J,FB,RC)
      integer:: n,I,J,JJ,  F(n),FB(n),RC(n)
end subroutine incr_inalborz
subroutine INIT_inalborz(N,A,F,M,U,V,FB,P)
      integer:: n,m, F(n),FB(n),P(n)
      real(8) A(n,n) , U(n),V(n)
      real(8), parameter :: INF = 1.d9
end subroutine init_inalborz
subroutine PATH_inalborz(N,A,II,F,JJ,U,V,FB,RC)
      integer:: N 
      real(8)::  A(n,n),U(n),V(N),PI(n), IA, MIN
      integer:: F(N),LR(n),UC(n)
      integer:: FB(n),RC(n)
      real(8), parameter :: INF = 1.d9
      integer::  i,j,k,L,ii,jj,NUC,NLR,R
end subroutine path_inalborz
! ./src/identical.F90 :
subroutine identical(parini,nlminx,nlmin,fp_method,fp_len,ent_wpos,fp_wpos,ent_arr,fp_arr,&
           &ent_delta,fp_delta,newmin,kid,fp_dist_min,k_e_wpos,n_unique,n_nonuni,lid,nid)
use mod_parini, only: typ_parini
type(typ_parini), intent(in):: parini
integer:: nlminx,nlmin,fp_len,kid,k_e_wpos,n_unique,n_nonuni
integer:: i,l,klow,k,khigh,fp_method,lid(nlminx),nid
real(8):: fp_arr(fp_len,nlminx),fp_wpos(fp_len),ent_arr(nlminx),ent_wpos,fp_delta,ent_delta,fp_dist_min,fp_dist
logical newmin,inrange
end subroutine identical
! ./src/inertia_tensor.F90 :
subroutine inertia_tensor(nat,xcart,cmass,amass,intens)
integer:: nat,iat,i,j
real(8):: xcart(3,nat),amass(nat),intens(3,3),cmass(3),xtmp(3),dist2
end subroutine inertia_tensor
! ./src/init_rotvels.F90 :
subroutine init_rotvels(parini,nat,xred,latvec,temp,amass,vel)
use mod_parini, only: typ_parini
type(typ_parini), intent(in):: parini
integer,intent(in):: nat
real(8),intent(in):: xred(3,nat),latvec(3,3),temp,amass(nat)
real(8),intent(out):: vel(3,nat)
real(8):: xcart(3,nat),ekin_rot,ekin_trans,rotmat(3,3),dproj(6),angbohr,erot_tmp,ekin_tot,v2gauss,vtest,tmp(3)
integer, dimension(nat):: fragarr,fragsize
end subroutine init_rotvels
subroutine assign_vel(nat,xcart,cmass,omega,vel)
integer:: nat,iat
real(8):: xcart(3,nat),cmass(3),xtmp(3),vel(3,nat),omega(3)
end subroutine assign_vel
subroutine rot_ener(omega,intens,erot)
real(8):: omega(3),intens(3,3),erot
end subroutine rot_ener
! ./src/init_vel.F90 :
subroutine init_vel(parini,parres,vel,vel_lat,vel_vol,latvec,pos_red,latmass,temp,nsoften,folder)
 use mod_parini, only: typ_parini
 type(typ_parini), intent(in):: parini
 type(typ_parini), intent(inout):: parres
 real(8) :: acell_in(3),xred_in(3,parini%nat),vel_in(3,parini%nat),fcart_in(3,parini%nat),etot_in,strten_in(6),vel_lat_in(3,3)
 real(8):: vel(3,parini%nat),temp,pos_red(3,parini%nat),vcm(3),vel_vol
 integer:: i,iat,idim,nsoften
 real(8):: amass(parini%nat),s1,s2,v2gauss,vtest,rescale_vel,vel_lat(3,3),latvec(3,3),latmass
 real(8), parameter :: temp_fac_lat=1.d-1 !This percentage of the temperature that should be given to the lattice 
 character(40):: folder
end subroutine init_vel
! ./src/insert.F90 :
subroutine insert(nlminx,nlmin,fp_len,nat,k_e_wpos,e_wpos,ent_wpos,fp_wpos,wpos_red,&
  &wpos_latvec,wpos_fcart,wpos_strten,spg_wpos,spgtol_wpos,fdos_wpos,&
  &e_arr,ent_arr,fp_arr,pl_arr,lat_arr,f_arr,str_arr,spg_arr,spgtol_arr,dos_arr,ct_arr)
  integer:: fp_len,ct_arr(nlminx),spg_arr(nlminx),nat,iat,spg_wpos
  integer:: k, nlmin, k_e_wpos, nlminx,i
  real(8):: e_wpos, ent_wpos, wpos_red(3,nat),wpos_latvec(3,3),spgtol_wpos,fdos_wpos,fp_wpos(fp_len)
  real(8):: e_arr(nlminx),ent_arr(nlminx),fp_arr(fp_len,nlminx),pl_arr(3,nat,nlminx),f_arr(3,nat,nlminx)
  real(8):: lat_arr(3,3,nlminx),spgtol_arr(nlminx),dos_arr(nlminx),fp1(fp_len),fp2(fp_len),str_arr(6,nlminx)
  real(8):: wpos_fcart(3,nat),wpos_strten(6)
end subroutine insert
! ./src/io_acf.F90 :
subroutine acf_write(file_info,atoms,atoms_all,strkey)
    use mod_atoms, only: typ_file_info, typ_atoms, typ_atoms_all, atom_copy_old
    type(typ_file_info), intent(inout):: file_info
    type(typ_atoms), optional, intent(in):: atoms
    type(typ_atoms_all), optional, intent(in):: atoms_all
    character(*), optional, intent(in):: strkey
end subroutine acf_write
subroutine acf_write_new(file_info,atoms_arr,strkey)
    use mod_atoms, only: typ_file_info, typ_atoms, typ_atoms_arr, atom_copy_old
    type(typ_file_info), intent(inout):: file_info
    type(typ_atoms_arr),intent(in):: atoms_arr
    character(*), optional, intent(in):: strkey
end subroutine acf_write_new
subroutine rotate4acf(nat,rat,cv,cvrot)
    integer, intent(in):: nat
    real(8), intent(inout):: rat(3,nat)
    real(8), intent(in):: cv(3,3), cvrot(3,3)
end subroutine rotate4acf
subroutine acf_force_write(file_info,atoms,atoms_all,strkey)
    use mod_atoms, only: typ_file_info, typ_atoms, typ_atoms_all, atom_copy_old
    type(typ_file_info), intent(inout):: file_info
    type(typ_atoms), optional, intent(in):: atoms
    type(typ_atoms_all), optional, intent(in):: atoms_all
    character(*), optional, intent(in):: strkey
end subroutine acf_force_write
subroutine acf_read(parini,filename,nconfmax,atoms,atoms_all)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, typ_atoms_all, atom_all_allocate, atom_copy_old
    type(typ_parini), intent(in):: parini
    character(*), intent(in):: filename
    integer, intent(in):: nconfmax
    type(typ_atoms), optional, intent(inout):: atoms
    type(typ_atoms_all), optional, intent(inout):: atoms_all
end subroutine acf_read
subroutine acf_read_new(parini,filename,nconfmax,atoms_arr)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, typ_atoms_arr, atom_allocate, atom_copy
    type(typ_parini), intent(in):: parini
    character(*), intent(in):: filename
    integer, intent(in):: nconfmax
    type(typ_atoms_arr), intent(inout):: atoms_arr
end subroutine acf_read_new
subroutine str_parse(str,line,atoms_t,c5)
    use mod_atoms, only: typ_atoms, typ_atoms_all
    character(*), intent(in):: str
    integer, intent(in):: line
    type(typ_atoms), intent(inout):: atoms_t
    character(*), intent(inout):: c5
end subroutine str_parse
subroutine str_motion2bemoved(str_motion,bemoved)
    character(*), intent(in):: str_motion
    logical, intent(inout):: bemoved(3)
end subroutine str_motion2bemoved
! ./src/io_ascii.F90 :
subroutine ascii_getsystem(parini,filename)
use mod_parini, only: typ_parini
type(typ_parini), intent(inout):: parini
character(40) :: filename
end subroutine ascii_getsystem
subroutine read_atomic_file_ascii(filename,nat,units,xred,latvec,fcart,strten,&
           &fixat,fixlat,readfix,fragarr,readfrag,printval1,printval2)
integer:: nat,natin,iat,ierror,io,n,k,fragarr(nat),fragarr_tmp,lhead(nat),llist(nat),nmol,m,l
logical:: fixat(nat),fixlat(7),readfix,reduced_tmp,readfrag
character(40):: filename,units
real(8):: pos(3,nat),xred(3,nat),latvec(3,3),dproj(6),strten(6),fcart(3,nat)
real(8):: angbohr,evhartree,enthalpy_at,printval1,printval2
end subroutine read_atomic_file_ascii
subroutine write_atomic_file_ascii(parini,filename,nat,units,xred,latvec0,fcart,strten,char_type,&
           &ntypat,typat,fixat,fixlat,energy,pressure,printval1,printval2)
use mod_parini, only: typ_parini
type(typ_parini), intent(in):: parini
integer:: nat,natin,iat,ntypat,typat(nat),j
character(40):: filename,units
real(8):: pos(3,nat),xred(3,nat),latvec(3,3),latvec0(3,3),dproj(6),rotmat(3,3),v(3,3),ucvol,fcart(3,nat),strten(6)
real(8):: energy, etotal, enthalpy, enthalpy_at,pressure,printval1,printval2,tmp(3)
character(2):: char_type(ntypat)
logical:: fixat(nat),fixlat(7)
end subroutine write_atomic_file_ascii
! ./src/io_bin.F90 :
subroutine read_bin_conf(parini,filename,atoms_arr)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms_arr, atom_allocate
    type(typ_parini), intent(in):: parini
    character(*), intent(in):: filename
    type(typ_atoms_arr), intent(inout):: atoms_arr
end subroutine read_bin_conf
subroutine read_bin_conf_v1(parini,filename,iunit,atoms_arr)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms_arr, iatom_to_sat, atom_allocate, update_rat
    type(typ_parini), intent(in):: parini
    character(*), intent(in):: filename
    integer, intent(in):: iunit
    type(typ_atoms_arr), intent(inout):: atoms_arr
end subroutine read_bin_conf_v1
subroutine write_bin_conf(file_info,atoms,strkey)
    use mod_atoms, only: typ_atoms, typ_file_info
    type(typ_file_info), intent(inout):: file_info
    type(typ_atoms), intent(in):: atoms
    character(*), optional, intent(in):: strkey
end subroutine write_bin_conf
subroutine write_bin_conf_v1(filename,file_position,iunit,atoms)
    use mod_atoms, only: typ_atoms, sat_to_iatom, get_rat
    character(*), intent(in):: filename
    character(*), intent(in):: file_position
    integer, intent(in):: iunit
    type(typ_atoms), intent(in):: atoms
end subroutine write_bin_conf_v1
! ./src/io_cube.F90 :
subroutine cube_read(filename,atoms,poisson)
    use mod_atoms, only: typ_atoms, iatom_to_sat, atom_allocate_old, update_rat
    use mod_electrostatics, only: typ_poisson
    character(*), intent(in):: filename
    type(typ_atoms), intent(out):: atoms
    type(typ_poisson), intent(out):: poisson
end subroutine cube_read
subroutine cube_write(filename,atoms,poisson,rho_or_pot)
    use mod_atoms, only: typ_atoms, sat_to_iatom, get_rat
    use mod_electrostatics, only: typ_poisson
    character(*), intent(in):: filename
    type(typ_atoms), intent(in):: atoms
    type(typ_poisson), intent(in):: poisson
    character(*), intent(in):: rho_or_pot
end subroutine cube_write
! ./src/io_utils.F90 :
subroutine read_list_files_yaml(fname,nfiles_max,fn_list,nfiles)
    character(len=*), intent(in):: fname
    integer, intent(in):: nfiles_max
    character(len=256), intent(out):: fn_list(nfiles_max)
    integer, intent(out):: nfiles
end subroutine read_list_files_yaml
! ./src/io_vasp.F90 :
subroutine write_poscar(filename,nat,rat,latvec,ntypat,natarr,comment,vasp5,comment2,atom_motion)
    integer:: nat, ntypat, natarr(128)
    character(*):: filename
    real(8):: rat(3,nat),latvec(3,3)
    character(*):: comment, comment2
    logical:: vasp5, atom_motion(3,nat)
end subroutine write_poscar
! ./src/io_vasp_minhocao.F90 :
subroutine poscar_getsystem(parini,filename)
use mod_parini, only: typ_parini
type(typ_parini), intent(inout):: parini
character(*) :: filename
end subroutine poscar_getsystem
subroutine write_atomic_file_poscar(parini,filename,nat,units,xred,latvec0,fcart,strten,char_type,&
           &ntypat,typat,fixat,fixlat,energy,pressure,printval1,printval2)
use mod_parini, only: typ_parini
type(typ_parini), intent(in):: parini
integer:: nat,natin,iat,ntypat,typat(nat)
character(*):: filename,units
real(8):: pos(3,nat),xred(3,nat),latvec(3,3),latvec0(3,3),dproj(6),rotmat(3,3),v(3,3),ucvol,fcart(3,nat),strten(6)
real(8):: energy, etotal, enthalpy, enthalpy_at,pressure,printval1,printval2
character(2):: char_type(ntypat)
logical:: fixat(nat),fixlat(7)
end subroutine write_atomic_file_poscar
subroutine read_atomic_file_poscar(filename,nat,units,xred,latvec,fcart,strten,&
           &fixat,fixlat,readfix,fragarr,readfrag,printval1,printval2)
integer:: i,ntypat_tmp,nat,natin,iat,ierror,io,n,k,fragarr(nat),fragarr_tmp,lhead(nat),llist(nat),nmol
logical:: fixat(nat),fixlat(7),readfix,reduced_tmp,readfrag
character(*):: filename,units
real(8):: pos(3,nat),xred(3,nat),latvec(3,3),dproj(6),strten(6),fcart(3,nat)
real(8):: angbohr,evhartree,enthalpy_at,printval1,printval2,scaling
end subroutine read_atomic_file_poscar
! ./src/io_xyz.F90 :
subroutine writexyz(filename,fn_position,nat,rat,bemoved,sat,cellvec,boundcond,comment)
    integer, intent(in):: nat
    real(8), intent(in):: rat(3,nat)
    logical, intent(in):: bemoved(3,nat)
    character(5), intent(in):: sat(nat)
    real(8), intent(in):: cellvec(3,3)
    character(*), intent(in):: filename, fn_position, boundcond, comment
end subroutine writexyz
subroutine readxyznat(filename,nat)
    integer:: nat
    character(*):: filename
end subroutine readxyznat
subroutine readxyz(filename,nat,rat,sat,comment1,comment2,atom_motion)
    integer:: nat
    real(8):: rat(3,nat)
    character(5):: sat(nat)
    character(*):: filename, comment1, comment2
    logical:: atom_motion(3,nat)
end subroutine readxyz
! ./src/io_yaml_conf.F90 :
subroutine read_yaml_conf(parini,filename,nconfmax,atoms_arr)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms_arr, atom_allocate, update_rat
    type(typ_parini), intent(in):: parini
    character(*), intent(in):: filename
    integer, intent(in):: nconfmax
    type(typ_atoms_arr), intent(inout):: atoms_arr
end subroutine read_yaml_conf
subroutine write_yaml_conf(file_info,atoms,strkey)
    use mod_atoms, only: typ_file_info, typ_atoms, update_ratp, get_rat
    type(typ_file_info), intent(inout):: file_info
    type(typ_atoms), intent(in):: atoms
    character(*), optional, intent(in):: strkey
end subroutine write_yaml_conf
! ./src/latticetools_minhocao.F90 :
 subroutine dist2line(point,ppoint1,ppoint2,dist)
 real(8), intent(in) :: point(3),ppoint1(3),ppoint2(3)
 real(8), intent(out):: dist
end subroutine dist2line
 subroutine dist2plane(point,nvec,ppoint,dist)
 real(8), intent(in) :: point(3),nvec(3),ppoint(3)
 real(8), intent(out):: dist
end subroutine dist2plane
subroutine dist_ang2latvec(dist_ang,latvec,pi)
real(8):: dist_ang(6),latvec(3,3),pi,convang
end subroutine dist_ang2latvec
subroutine dist_latvec2ang(dist_ang,latvec,pi)
real(8):: dist_ang(6),latvec(3,3),pi,convang
end subroutine dist_latvec2ang
 subroutine dproj2latvec(dproj,latvec)
 real*8:: dproj(6),latvec(3,3)
end subroutine dproj2latvec
 subroutine latvec2dproj(dproj,latvec,rotmat,rxyz,nat)
 integer,intent(in)  :: nat
 real*8,intent(inout):: dproj(6),latvec(3,3),rotmat(3,3),rxyz(3,nat)
end subroutine latvec2dproj
 subroutine latvec2acell_rprim(latvec,acell,rprim)
 real(8):: latvec(3,3), rprim(3,3), acell(3)
end subroutine latvec2acell_rprim
subroutine getvol(latvec,vol)
real(8):: latvec(3,3),v(3,3),vol
end subroutine getvol
 subroutine acell_rprim2latvec(latvec,acell,rprim)
 real(8):: latvec(3,3), rprim(3,3), acell(3)
end subroutine acell_rprim2latvec
 subroutine nveclatvec(latvec,nvec)
 real*8, intent(in) :: latvec(3,3)
 real*8, intent(out):: nvec(3,3)
end subroutine nveclatvec
 subroutine rxyz_cart2int(latvec,rxyzint,rxyzcart,nat)
 real(8):: rxyzint(3,nat), rxyzcart(3,nat),latvec(3,3),latvecinv(3,3)
 integer:: nat,iat
end subroutine rxyz_cart2int
 subroutine rxyz_int2cart(latvec,rxyzint,rxyzcart,nat)
 real(8):: rxyzint(3,nat), rxyzcart(3,nat),latvec(3,3)
 integer:: nat,iat
end subroutine rxyz_int2cart
 subroutine fxyz_cart2int(nat,fxyz_cart,fxyz_int,latvec)
 real(8):: fxyz_cart(3,nat),fxyz_int(3,nat),latvec(3,3),transmat(3,3)
 integer:: nat,iat
end subroutine fxyz_cart2int
subroutine k_expansion(parini,latvec,xred,ka,kb,kc,k_latvec,k_xcart)
use mod_parini, only: typ_parini
type(typ_parini), intent(in):: parini
real(8):: latvec(3,3),k_latvec(3,3),k_xcart(3,parini%nat,ka,kb,kc),xred(3,parini%nat) 
integer:: iat,k,l,m,ka,kb,kc
end subroutine k_expansion
 subroutine strten2flat(strten,flat,latvec,press)
 real(8):: strten(6),flat(3,3),latvec(3,3),press,pressmat(3,3),str_matrix(3,3),latvect(3,3),latvectinv(3,3),vol
end subroutine strten2flat
subroutine rotate_stresstensor(strten,rotmat)
real(8):: strten(6),rotmat(3,3),stress(3,3)
end subroutine rotate_stresstensor
 subroutine find_kpt(k1, k2, k3, lat, gridden)
   integer, intent(out) :: k1,k2,k3
   real(8), intent(in)  :: lat(3,3), gridden
end subroutine find_kpt
   subroutine track_kpt(gridden, glen, kpt)
     real(8), intent(in) :: gridden, glen
     integer :: kpt,j
end subroutine track_kpt
subroutine rotmat_fcart_stress(latvec_init,latvec_trans,rotmat)
real(8):: latvec_init(3,3),latvec_trans(3,3),latvec_trans_inv(3,3),rotmat(3,3)
end subroutine rotmat_fcart_stress
 subroutine updaterxyz(latvecold,latvecnew,rxyz,nat)
 real(8), intent(in)   :: latvecold(3,3), latvecnew(3,3)
 real(8), intent(inout):: rxyz(3,nat)
 integer, intent(in)   :: nat
end subroutine updaterxyz
! ./src/lenosky_tightbinding.F90 :
subroutine lenoskytb_alborz(parini,atoms,natsi,count_md)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_potl, only: potl_typ
    use mod_tightbinding, only: typ_partb, lenosky
    use mod_linked_lists, only: typ_pia_arr, typ_linked_lists
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
    integer, intent(in):: natsi
    real(8), intent(inout):: count_md
end subroutine lenoskytb_alborz
subroutine lenoskytb_init(partb,atoms,natsi,linked_lists)
    use mod_tightbinding, only: typ_partb, lenosky
    use mod_atoms, only: typ_atoms
    use mod_linked_lists, only: typ_linked_lists
    type(typ_partb), intent(inout):: partb
    type(typ_atoms), intent(in):: atoms
    integer, intent(in):: natsi
    type(typ_linked_lists), intent(in):: linked_lists
end subroutine lenoskytb_init
subroutine totalenergy(pia_arr,linked_lists,parini,partb,atoms,natsi,pplocal)
    use mod_linked_lists, only: typ_pia_arr, typ_linked_lists
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_potl, only: potl_typ
    use mod_tightbinding, only: typ_partb, lenosky
    type(typ_linked_lists), intent(in):: linked_lists
    type(typ_pia_arr), intent(in):: pia_arr
    type(typ_parini), intent(in):: parini
    type(typ_partb), intent(inout):: partb
    type(typ_atoms), intent(inout):: atoms
    integer, intent(in):: natsi
    type(potl_typ), intent(inout):: pplocal
end subroutine totalenergy
subroutine pairenergy(parini,partb,atoms,pplocal,natsi)
    use mod_parini, only: typ_parini
    use mod_tightbinding, only: typ_partb
    use mod_atoms, only: typ_atoms
    use mod_potl, only: potl_typ
    use mod_linked_lists, only: typ_pia_arr, typ_linked_lists
    type(typ_parini), intent(in):: parini
    type(typ_partb), intent(inout):: partb
    type(typ_atoms), intent(inout):: atoms
    integer, intent(in):: natsi
    type(potl_typ), intent(in):: pplocal
end subroutine pairenergy
subroutine lenoskytb_final(partb)
    use mod_tightbinding, only: typ_partb
    type(typ_partb), intent(inout):: partb
end subroutine lenoskytb_final
subroutine radelmgeneralsp(r,radar,dradar,atomtypei,atomtypej,pplocal)
    use mod_potl, only: potl_typ
    type(potl_typ), intent(in):: pplocal
    real(8), intent(in):: r
    real(8), intent(out):: radar(0:3), dradar(0:3)
    integer, intent(in):: atomtypei, atomtypej
end subroutine radelmgeneralsp
subroutine clssplint(str_action,s,xt,yt,derivt,extype)
    use mod_splinetb, only: NSPMAX, spline_typ
    character(*), intent(in) :: str_action
    type(spline_typ), intent(in) :: s
    integer, intent(in)::  extype 
    real(8), intent(in) :: xt
    real(8), intent(out) :: yt, derivt
end subroutine clssplint
subroutine eselfgeneral(eself)
    real(8), intent(inout):: eself(0:3)
end subroutine eselfgeneral
subroutine prmst38c(partb,pplocal)
    use mod_tightbinding, only: typ_partb
    use mod_potl, only: potl_typ
    type(typ_partb):: partb
    type(potl_typ):: pplocal
end subroutine prmst38c
subroutine clsfread_spline(unit,s)
    use mod_splinetb, only: NSPMAX, spline_typ
    integer:: unit
    type(spline_typ), intent(out):: s
end subroutine clsfread_spline
subroutine clsspline(s)
    use mod_splinetb, only: NSPMAX, spline_typ
    type(spline_typ), intent(inout):: s
end subroutine clsspline
! ./src/linked_lists_mod.F90 :
! ./src/lmder_modified.F90 :
subroutine init_lmder_modified(parlm,m,ldfjac)
    use mod_parlm, only: typ_parlm
    type(typ_parlm), intent(inout):: parlm
    integer, intent(in):: m, ldfjac
end subroutine init_lmder_modified
subroutine final_lmder_modified(parlm)
    use mod_parlm, only: typ_parlm
    type(typ_parlm), intent(inout):: parlm
end subroutine final_lmder_modified
subroutine lmder_modified(parlm,m,ldfjac)
    use mod_parlm, only: typ_parlm
      type(typ_parlm), intent(inout):: parlm
      integer:: m,ldfjac
end subroutine lmder_modified
! ./src/logo_minhocao.F90 :
subroutine print_logo()
end subroutine print_logo
! ./src/md.F90 :
subroutine dynamics(parini)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, typ_file_info, set_ndof, atom_deallocate_old
    type(typ_parini), intent(inout):: parini
end subroutine dynamics
subroutine md_nve(parini,atoms)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, typ_file_info, atom_copy_old, set_atomic_mass
    type(typ_parini), intent(in):: parini
    type(typ_atoms):: atoms, atoms_old
end subroutine md_nve
subroutine md_nph(parini,atoms)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, typ_file_info, set_rat, update_ratp, update_rat
    type(typ_parini), intent(in):: parini
    type(typ_atoms):: atoms
end subroutine md_nph
subroutine set_velocities(atoms, ekin_arg)
    use mod_atoms, only: typ_atoms, update_ratp
    type(typ_atoms), intent(inout):: atoms
    real(8), optional::ekin_arg
end subroutine set_velocities
subroutine ekin_temprature(atoms,temp,vcm,rcm,totmass) 
    use mod_atoms, only: typ_atoms, update_ratp
    type(typ_atoms):: atoms
    real(8) :: vcm(3), rcm(3), t1, tmp
    real(8) :: aboltzmann, totmass, temp 
end subroutine ekin_temprature
! ./src/md_minhocao_andersen.F90 :
subroutine MD_ANDERSEN_MHM     (parini,parres,latvec_in,xred_in,fcart_in,strten_in,vel_in,vel_lat_in,vvol_in,etot_in,iprec,counter,folder)
 use mod_parini, only: typ_parini
 type(typ_parini), intent(in):: parini
 type(typ_parini), intent(inout):: parres
 real(8) :: latvec_in(3,3),xred_in(3,parini%nat),vel_in(3,parini%nat),fcart_in(3,parini%nat),etot_in,strten_in(6),vel_lat_in,vvol_in
 real(8):: amass(parini%nat)
 real(8),dimension(3,parini%nat):: xcart
 real(8),dimension(3,parini%nat):: fposcur
 real(8),dimension(3,parini%nat):: accposcur
 real(8),dimension(3,parini%nat):: accpospred
 real(8),dimension(3,parini%nat):: accposprev
 real(8),dimension(3,parini%nat):: fpospred
 real(8),dimension(3,parini%nat):: vpospred
 real(8),dimension(3,parini%nat):: poscur
 real(8),dimension(3,parini%nat):: vxyz
 real(8),dimension(3,parini%nat):: vposcur
 real(8),dimension(3,parini%nat):: pospred
 real(8),dimension(3,parini%nat):: dxred
 real(8):: counter
 integer:: iprec 
 character(40)::filename,folder
end subroutine md_andersen_mhm
subroutine ekin_at_lat_andersen(amass,latmass,latvec,vpos,vlat,vvol,ekinat,ekinlat,f0,md_type,nat)
integer:: iat,i,md_type,nat
real(8):: latvec(3,3),vpos(3,nat),vlat(3,3),ekinat,ekinlat,rkin,vposcurtmp(3),crossp(3),f0(3,3),vol,vvol
real(8):: latmass,amass(nat),lattrans(3,3),latdottrans(3,3),ekintrace(3,3),sigma(3,3),sigmatrans(3,3),vol_1_3
end subroutine ekin_at_lat_andersen
! ./src/md_minhocao.F90 :
subroutine MD_MHM   (parini,parres,latvec_in,xred_in,fcart_in,strten_in,vel_in,vel_lat_in,vvol_in,etot_in,iprec,counter,folder)
 use mod_parini, only: typ_parini
 type(typ_parini), intent(in):: parini
 type(typ_parini), intent(inout):: parres
 real(8) :: latvec_in(3,3),xred_in(3,parini%nat),vel_in(3,parini%nat),fcart_in(3,parini%nat),etot_in,strten_in(6),vel_lat_in(3,3),vvol_in
 real(8):: amass(parini%nat)
 real(8),dimension(3,parini%nat):: xcart
 real(8),dimension(3,parini%nat):: fposcur
 real(8),dimension(3,parini%nat):: accposcur
 real(8),dimension(3,parini%nat):: accpospred
 real(8),dimension(3,parini%nat):: accposprev
 real(8),dimension(3,parini%nat):: fpospred
 real(8),dimension(3,parini%nat):: vpospred
 real(8),dimension(3,parini%nat):: poscur
 real(8),dimension(3,parini%nat):: vxyz
 real(8),dimension(3,parini%nat):: vposcur
 real(8),dimension(3,parini%nat):: pospred
 real(8),dimension(3,parini%nat):: dxred
 real(8),dimension(parres%nmd_dynamics):: ensave   !zl
 real(8),dimension(parres%nmd_dynamics):: ensmoth  !zl
 real(8):: counter
 integer:: iprec 
 character(40)::filename,folder
end subroutine md_mhm
subroutine ekin_at_lat(amass,latmass,latvec,vpos,vlat,ekinat,ekinlat,f0,md_type,nat)
integer:: iat,i,md_type,nat
real(8):: latvec(3,3),vpos(3,nat),vlat(3,3),ekinat,ekinlat,rkin,vposcurtmp(3),crossp(3),f0(3,3),vol
real(8):: latmass,amass(nat),lattrans(3,3),latdottrans(3,3),ekintrace(3,3),sigma(3,3),sigmatrans(3,3)
end subroutine ekin_at_lat
subroutine stress_velocity(vpos,latvec,amass,nat,vpressure)
real(8):: velmat(3,3),vpostmp(3),latvec(3,3),vpos(3,nat),vpressure,a(3,3),vol,amass(nat)
integer:: iat,nat,i,j
end subroutine stress_velocity
! ./src/md_minhocao_rbmd.F90 :
subroutine MD_MHM_ROT(parini,parres,latvec_in,xred_in,xred_cm_in,xcart_mol,quat_in,fcart_in,strten_in,&
                      &vel_in,vel_cm_in,vel_lat_in,l_in,vvol_in,etot_in,&
                      &masstot,intens,inprin,inaxis,lhead,llist,nmol,iprec,counter,folder)
 use mod_parini, only: typ_parini
 type(typ_parini), intent(in):: parini
 type(typ_parini), intent(inout):: parres
 real(8) ::latvec_in(3,3),xred_in(3,parini%nat),vel_in(3,parini%nat),fcart_in(3,parini%nat),etot_in,strten_in(6),vel_lat_in(3,3),vvol_in
 real(8):: amass(nmol)
 real(8),dimension(3,nmol):: xred_cm_in
 real(8),dimension(3,nmol):: fcart_cm
 real(8),dimension(3,nmol):: torque
 real(8),dimension(3,parini%nat) :: xcart_mol
 real(8),dimension(3,nmol):: l_in
 real(8),dimension(3,nmol):: vel_cm_in
 real(8),dimension(4,nmol):: quat_in
 real(8),dimension(4,nmol):: quatcur
 real(8),dimension(4,nmol):: quatpred
 real(8),dimension(3,3,nmol):: intens
 real(8),dimension(3,nmol):: inprin
 real(8),dimension(nmol):: masstot
 real(8),dimension(3,3,nmol):: inaxis
 integer,dimension(nmol):: lhead
 integer,dimension(parini%nat):: llist
 integer:: nmol
 real(8),dimension(3,nmol):: xcart
 real(8),dimension(3,nmol):: fposcur
 real(8),dimension(3,nmol):: accposcur
 real(8),dimension(3,nmol):: accpospred
 real(8),dimension(3,nmol):: accposprev
 real(8),dimension(3,nmol):: fpospred
 real(8),dimension(3,nmol):: vpospred
 real(8),dimension(3,nmol):: poscur
 real(8),dimension(3,nmol):: vxyz
 real(8),dimension(3,nmol):: vposcur
 real(8),dimension(3,parini%nat):: pospred
 real(8),dimension(3,parini%nat):: dxred
 real(8):: counter
 integer:: iprec 
 character(40)::filename,folder
end subroutine md_mhm_rot
subroutine rbmd_symasym_s1(T_t,L_t,dt,L_til_t)
real(8),intent(in) :: T_t(3), L_t(3),dt
real(8),intent(out):: L_til_t(3)
end subroutine rbmd_symasym_s1
subroutine rbmd_sym_s23(Inprin,L_til_t,dt,L_til_t5,L_til_t10)
real(8),intent(in) :: Inprin(3), L_til_t(3),dt
real(8),intent(out):: L_til_t10(3),L_til_t5(3)
end subroutine rbmd_sym_s23
subroutine rbmd_symasym_s4(Inprin,L_til_t5,quat_t,dt,quat_t10)
real(8),intent(in) :: Inprin(3),L_til_t5(3),dt,quat_t(4)
real(8),intent(out):: quat_t10(4)
end subroutine rbmd_symasym_s4
function A_omega(omega)
real(8),dimension(4,4):: A_omega
real(8)::omega(3)
end function a_omega
subroutine rbmd_symasym_s5(T_t10,L_til_t10,dt,L_t10)
real(8),intent(in) :: T_t10(3), L_til_t10(3),dt
real(8),intent(out):: L_t10(3)
end subroutine rbmd_symasym_s5
subroutine rbmd_asym_s23(Inprin,L_til_t,dt,L_til_t5,L_til_t10)
real(8),intent(in) :: Inprin(3), L_til_t(3),dt
real(8),intent(out):: L_til_t10(3),L_til_t5(3)
end subroutine rbmd_asym_s23
subroutine rbmd_driver(quat_t,T_t,L_t,quat_t10,T_t10,L_t10,dt,inprin,&
           &fragsize,symtop,nmol)
real(8),intent(in) :: inprin(3,nmol),L_t(3,nmol),T_t(3,nmol),quat_t(4,nmol),dt,T_t10(3,nmol)
integer,intent(in) :: nmol,fragsize(nmol)
logical,intent(in) :: symtop(nmol)
real(8),intent(out):: L_t10(3,nmol),quat_t10(4,nmol)
real(8) :: L_til_t(3,nmol),L_til_t5(3,nmol),L_til_t10(3,nmol)
end subroutine rbmd_driver
subroutine expand_rigid(latvec,xred_cm,quat,xcart_mol,lhead,llist,nat,nmol,xred_in)
real(8),intent(in):: latvec(3,3),xred_cm(3,nmol),quat(4,nmol),xcart_mol(3,nat)
real(8),intent(out):: xred_in(3,nat)
real(8):: rotmat(3,3),xcart_tmp(3,nat),cmass(3,nmol)
integer:: nat,nmol,iat,imol,llist(nat),lhead(nmol)
end subroutine expand_rigid
! ./src/md_NVT.F90 :
subroutine md_nvt_langevin(parini,atoms)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, typ_file_info
    type(typ_parini), intent(inout):: parini
    type(typ_atoms):: atoms
    real(8):: eta(3,atoms%nat)
    real(8):: langev(atoms%nat), forces_langevin(3,atoms%nat)
    real(8):: rat_next(3,atoms%nat), vat_old(3,atoms%nat)
    real(8):: rat_init(3,atoms%nat)
end subroutine md_nvt_langevin
subroutine md_nvt_nose_hoover_cp(parini,atoms)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, typ_file_info, set_rat, get_rat, update_ratp
    type(typ_parini), intent(inout):: parini
    type(typ_atoms):: atoms
    real(8):: forces_nosehoover(3,atoms%nat)
    real(8):: rat_next(3,atoms%nat), rat_prev(3,atoms%nat),vat_old(3,atoms%nat) 
    real(8):: rat_init(3,atoms%nat)
end subroutine md_nvt_nose_hoover_cp
subroutine md_nvt_nose_hoover_chain(parini,atoms)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, typ_file_info, get_rat, update_ratp, update_rat
    type(typ_parini), intent(inout):: parini
    type(typ_atoms):: atoms
    real(8):: rat_init(3,atoms%nat)
    integer:: jj(3,atoms%nat), vfile
end subroutine md_nvt_nose_hoover_chain
subroutine set_langevin_randforce(eta,nat)
    integer :: nat
    real(8) ::eta(3,nat), sum1, sum2, sum3
end subroutine set_langevin_randforce
subroutine back_to_cell(atoms)
    use mod_atoms, only: typ_atoms, update_ratp, update_rat
    type(typ_atoms):: atoms
end subroutine back_to_cell
subroutine plane_repulsion(atoms)
    use mod_atoms, only: typ_atoms, update_ratp
    type(typ_atoms):: atoms
end subroutine plane_repulsion
subroutine thermostat_evolution(atoms,zeta_next,zeta,zeta_prev,dzeta,mass_q,kt,ntherm,imd)
    use mod_atoms, only: typ_atoms, typ_file_info
    type(typ_atoms):: atoms
    integer :: ntherm, imd 
    real(8) :: kt, t1
    real(8) :: zeta_next(3,atoms%nat,ntherm), zeta(3,atoms%nat,ntherm),zeta_prev(3,atoms%nat,ntherm)
    real(8) :: dzeta(3,atoms%nat,ntherm), mass_q(ntherm)
    real(8) :: force_therm(3,atoms%nat,ntherm)
end subroutine thermostat_evolution
subroutine get_atomic_mass(atoms,totmass)
    use mod_atoms, only: typ_atoms
    type(typ_atoms):: atoms
    real(8):: totmass,mass_conv = 1822.888484264545
end subroutine get_atomic_mass
subroutine write_trajectory_velocity(parini,atoms,file_info,rat_init,imd,ntherm,zeta,dzeta)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, typ_file_info, update_ratp
    type(typ_parini), intent(inout):: parini
    type(typ_atoms):: atoms
    type(typ_file_info):: file_info
    integer:: imd, ntherm, ith
    real(8):: zeta(ntherm), dzeta(ntherm)
    real(8):: rat_init(3,atoms%nat)
end subroutine write_trajectory_velocity
! ./src/minhocao_enthalpyrelax.F90 :
subroutine enthalpyrelax(parini,parres,latvec,xred,tolmin,tolmax,ntol,findsym)
 use mod_parini, only: typ_parini
       type(typ_parini), intent(in):: parini
       type(typ_parini), intent(inout):: parres
       real(8):: latvec(3,3),xred(3,parini%nat),pinit,pfinal,psteps,pcur,latvec0(3,3),xred0(3,parini%nat),vel_vol_in
       real(8):: vel_in(3,parini%nat),vel_lat_in(3,3),fcart(3,parini%nat),strten(6)
       real(8):: tolmin,tolmax,spgtol
       integer:: ntol,spgint
       logical:: findsym
end subroutine enthalpyrelax
! ./src/minhocao_pathintegral.F90 :
subroutine pathintegral(parini,parres,latvec,xred)
 use mod_parini, only: typ_parini
       type(typ_parini), intent(in):: parini
       type(typ_parini), intent(inout):: parres
       real(8):: rxyz0(3,parini%nat),fxyz(3,parini%nat),displ(3,parini%nat)
       real(8):: evals(3),s2(3,3),dmat(3,3),dproj(6),rotmat(3,3),xred(3,parini%nat)
       real(8):: stepsize_at,stepsize_lat,t1,t2,t3,path,xred_in(3,parini%nat),latvec_in(3,3),strten_in(6)
       real(8):: str_matrix(3,3),transformed(3,3),transformed_inv(3,3),fcart_in(3,parini%nat)
       real(8):: dlat(6),latvec(3,3),latvecinv(3,3),stress(3,3),displat(3,3),tstress(3,3)
end subroutine pathintegral
! ./src/minhocao_plot_fp_grid.F90 :
subroutine plot_fp_grid(parini,nlminx,nlmin,nat,fp_len,fp_arr,lat_arr,pl_arr)
use mod_parini, only: typ_parini
type(typ_parini), intent(in):: parini
integer:: nlminx,nlmin,fp_len,i,kk,nat
real(8):: fp_arr(fp_len,nlminx),fp_dist
real(8):: tmp_acell(3),tmp_real,tmp_rprim(3,3),lat_arr(3,3,nlminx),pl_arr(3,nat,nlminx),randpos(3)
end subroutine plot_fp_grid
! ./src/minhocao_poslowrelax.F90 :
subroutine poslowrelax(parini,parres,latvec,xred,tolmin,tolmax,ntol)
 use mod_parini, only: typ_parini
       type(typ_parini), intent(in):: parini
       type(typ_parini), intent(inout):: parres
       integer::  iprec,nstruct,i,ntol,spgint
       real(8):: latvec(3,3),xred(3,parini%nat),pinit,pfinal,psteps,latvec0(3,3),xred0(3,parini%nat),vel_vol_in
       real(8):: vel_in(3,parini%nat),vel_lat_in(3,3),fcart(3,parini%nat),strten(6)
       real(8):: counter,count_geopt,enthalpy,energy,vol,ext_press,tolmin,tolmax,spgtol_pos,spg_pos
end subroutine poslowrelax
! ./src/minhocao_rotate_like_crazy.F90 :
subroutine rotate_like_crazy(parini,parres,latvec,xred,tolmin,tolmax,ntol)
 use mod_parini, only: typ_parini
       type(typ_parini), intent(in):: parini
       type(typ_parini), intent(inout):: parres
       integer::  iprec,nstruct,i,ntol,spgint
       real(8):: latvec(3,3),xred(3,parini%nat),pinit,pfinal,psteps,latvec0(3,3),xred0(3,parini%nat),vel_vol_in
       real(8):: vel_in(3,parini%nat),vel_lat_in(3,3),fcart(3,parini%nat),strten(6),axis(3),rotmat(3,3),angle
       real(8):: counter,count_geopt,enthalpy,energy,vol,ext_press,tolmin,tolmax,spgtol_pos,spg_pos
end subroutine rotate_like_crazy
! ./src/minhocao_varvol.F90 :
subroutine varvol(parini,parres,latvec,xred,tolmin,tolmax,ntol,findsym)
 use mod_parini, only: typ_parini
       type(typ_parini), intent(in):: parini
       type(typ_parini), intent(inout):: parres
       real(8):: latvec(3,3),xred(3,parini%nat),latvec0(3,3),xred0(3,parini%nat)
       real(8):: vel_in(3,parini%nat),vel_lat_in(3,3),fcart(3,parini%nat),strten(6)
       real(8):: tolmin,tolmax,spgtol
       integer:: ntol,spgint,itime
       logical:: findsym,is_percentage
end subroutine varvol
! ./src/minhopp_allocation.F90 :
subroutine allocate_minhopp_arrays1(nproc)
    integer, intent(in):: nproc
end subroutine allocate_minhopp_arrays1
subroutine allocate_minhopp_arrays2(nat,nproc)
    integer, intent(in):: nat, nproc
end subroutine allocate_minhopp_arrays2
subroutine deallocate_minhopp_arrays
end subroutine deallocate_minhopp_arrays
! ./src/minhopp.F90 :
subroutine minimahopping(parini)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, typ_atoms_arr, atom_copy
    use mod_opt, only: typ_paropt
    type(typ_parini), intent(in):: parini
end subroutine minimahopping
subroutine init_minimahopping(parini,atoms_curr,atoms_hopp,atoms_allproc,atoms_locmin,paropt,paropt_prec)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, typ_atoms_arr, atom_copy, set_ndof
    use mod_opt, only: typ_paropt
    type(typ_parini), intent(in):: parini
    type(typ_paropt), intent(inout):: paropt_prec, paropt
    type(typ_atoms), intent(inout):: atoms_curr, atoms_hopp
    type(typ_atoms_arr), intent(inout):: atoms_allproc
    type(typ_atoms_arr), intent(inout):: atoms_locmin
end subroutine init_minimahopping
subroutine final_minimahopping(parini,atoms_curr,atoms_hopp,atoms_allproc,atoms_locmin,paropt,paropt_prec)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, typ_atoms_arr, atom_deallocate
    use mod_opt, only: typ_paropt
    type(typ_parini), intent(in):: parini
    type(typ_paropt), intent(inout):: paropt_prec, paropt
    type(typ_atoms), intent(inout):: atoms_curr, atoms_hopp
    type(typ_atoms_arr), intent(inout):: atoms_allproc
    type(typ_atoms_arr), intent(inout):: atoms_locmin
end subroutine final_minimahopping
subroutine set_amass(atoms_hopp)
    use mod_atoms, only: typ_atoms
    type(typ_atoms), intent(inout):: atoms_hopp
end subroutine set_amass
subroutine relax_minhopp(parini,atoms,paropt_prec,paropt)
    use mod_parini, only: typ_parini
    use mod_opt, only: typ_paropt
    use mod_atoms, only: typ_atoms, get_rat, set_rat
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
    type(typ_paropt), intent(inout):: paropt_prec, paropt
end subroutine relax_minhopp
subroutine print_minhopp_parameters
end subroutine print_minhopp_parameters
subroutine read_earr
end subroutine read_earr
subroutine readnat(atoms_curr)
    use mod_atoms, only: typ_atoms
    type(typ_atoms):: atoms_curr
end subroutine readnat
subroutine read_poscur_alborz(atoms_curr,atoms_allproc)
    use mod_atoms, only: typ_atoms, typ_atoms_arr, atom_copy
    type(typ_atoms):: atoms_curr
    type(typ_atoms_arr):: atoms_allproc
end subroutine read_poscur_alborz
subroutine read_minhopp_parameters 
end subroutine read_minhopp_parameters
subroutine minhopp_newrun_initialization(atoms_curr,atoms_locmin)
    use mod_atoms, only: typ_atoms, typ_atoms_arr, atom_copy
    type(typ_atoms), intent(in):: atoms_curr
    type(typ_atoms_arr), intent(inout):: atoms_locmin
end subroutine minhopp_newrun_initialization
subroutine read_poslow(atoms_locmin)
    use mod_atoms, only: typ_atoms_arr
    type(typ_atoms_arr), intent(inout):: atoms_locmin
end subroutine read_poslow
subroutine send_minimum_to_all(atoms_curr)
    use mod_atoms, only: typ_atoms, get_rat
    type(typ_atoms), intent(in):: atoms_curr
end subroutine send_minimum_to_all
subroutine send_minhopp_parameters_to_all(atoms_curr)
    use mod_atoms, only: typ_atoms, get_rat
    type(typ_atoms), intent(in):: atoms_curr
end subroutine send_minhopp_parameters_to_all
subroutine mdescape(parini,atoms_hopp)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, typ_file_info, update_rat, update_ratp
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms_hopp
end subroutine mdescape
subroutine collect_data_from_all_processors(ntry,atoms_curr,atoms_allproc,atoms_locmin)
    use mod_atoms, only: typ_atoms, typ_atoms_arr, atom_copy
    integer, intent(in):: ntry
    type(typ_atoms), intent(in):: atoms_curr
    type(typ_atoms_arr), intent(inout):: atoms_allproc, atoms_locmin
end subroutine collect_data_from_all_processors
subroutine request_receive(atoms_allproc)
    use mod_atoms, only: typ_atoms_arr
    type(typ_atoms_arr), intent(in):: atoms_allproc
end subroutine request_receive
subroutine test_receive(atoms_allproc,atoms_locmin)
    use mod_atoms, only: typ_atoms_arr, set_rat
    type(typ_atoms_arr), intent(inout):: atoms_allproc, atoms_locmin
end subroutine test_receive
subroutine cancel_excessive_irecv
end subroutine cancel_excessive_irecv
subroutine insert_alborz(kepos,epos)
    integer, intent(in):: kepos
    real(8), intent(in):: epos
end subroutine insert_alborz
subroutine save_low_conf_alborz(atoms,atoms_locmin)
    use mod_atoms, only: typ_atoms, typ_atoms_arr, atom_copy
    type(typ_atoms), intent(in):: atoms
    type(typ_atoms_arr), intent(inout):: atoms_locmin
end subroutine save_low_conf_alborz
subroutine velopt(parini,atoms)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, update_ratp
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
end subroutine velopt
subroutine soften(parini,nstep,atoms0,count_soften,count_soften_tot)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, typ_file_info, atom_deallocate, atom_copy
    type(typ_parini), intent(in):: parini
    integer, intent(in):: nstep
    type(typ_atoms), intent(inout):: atoms0
    real(8), intent(inout):: count_soften, count_soften_tot
end subroutine soften
subroutine write_minhopp(atoms_allproc,atoms_locmin)
    use mod_atoms, only: typ_atoms_arr, typ_file_info
    type(typ_atoms_arr), intent(inout):: atoms_allproc, atoms_locmin
end subroutine write_minhopp
subroutine write_minhopp_parameters
end subroutine write_minhopp_parameters
subroutine write_earr
end subroutine write_earr
subroutine escape_failed(parini,erat,erathopp)
    use mod_parini, only: typ_parini
    type(typ_parini), intent(in):: parini
    real(8), intent(in):: erat, erathopp
end subroutine escape_failed
subroutine local_minimum_accepted(atoms_hopp,atoms_curr,atoms_locmin)
    use mod_atoms, only: typ_atoms, typ_atoms_arr, atom_copy
    type(typ_atoms), intent(in):: atoms_hopp
    type(typ_atoms), intent(inout):: atoms_curr
    type(typ_atoms_arr), intent(inout):: atoms_locmin
end subroutine local_minimum_accepted
subroutine local_minimum_rejected(atoms_hopp)
    use mod_atoms, only: typ_atoms
    type(typ_atoms), intent(in):: atoms_hopp
end subroutine local_minimum_rejected
subroutine report_minhopp_iteration_info(atoms_curr)
    use mod_atoms, only: typ_atoms, typ_atoms_arr
    type(typ_atoms), intent(in):: atoms_curr
end subroutine report_minhopp_iteration_info
subroutine already_visited_minimum(parini)
    use mod_parini, only: typ_parini
    type(typ_parini), intent(in):: parini
end subroutine already_visited_minimum
subroutine new_minimum(atoms_hopp)
    use mod_atoms, only: typ_atoms, typ_file_info
    type(typ_atoms), intent(in):: atoms_hopp
end subroutine new_minimum
subroutine print_final_statistics
end subroutine print_final_statistics
subroutine MPI_atom_arr_copy(nat,atoms_arr)
    use mod_atoms, only: typ_atoms, typ_atoms_arr, get_rat, set_rat
    integer,intent(in)::nat
    type(typ_atoms_arr),intent(inout):: atoms_arr
end subroutine mpi_atom_arr_copy
! ./src/minhopp_mod.F90 :
! ./src/minhopp_pot.F90 :
subroutine setpot_init(parini,atoms_curr,paropt,paropt_prec)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_opt, only: typ_paropt
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms_curr
    type(typ_paropt), intent(inout):: paropt, paropt_prec
end subroutine setpot_init
subroutine setpot_final(parini,atoms_curr)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(in):: atoms_curr
end subroutine setpot_final
subroutine setpot_mdescape
end subroutine setpot_mdescape
subroutine setpot_soften
end subroutine setpot_soften
subroutine setpot_geopt_prec
end subroutine setpot_geopt_prec
subroutine setpot_geopt
end subroutine setpot_geopt
! ./src/mpi_utilities.F90 :
subroutine cal_matvec_mpi(n,p,g,v1)
    integer, intent(in):: n
    real(8), intent(in):: p(n,n), g(n)
    real(8), intent(out):: v1(n)
end subroutine cal_matvec_mpi
! ./src/optimizer_bfgs.F90 :
subroutine mybfgs(iproc,nr,x,epot,f,nwork,work,paropt)
    use mod_opt, only: typ_paropt, frmt_base
    integer, intent(in):: iproc, nr, nwork
    real(8):: x(nr), f(nr), epot, work(nwork)
    type(typ_paropt), intent(inout):: paropt
    character(51), parameter:: frmt='('//frmt_base//',2es12.4,i3,a)'
end subroutine mybfgs
subroutine init_mybfgs(paropt,epot,fmax)
    use mod_opt, only: typ_paropt
    type(typ_paropt), intent(inout):: paropt
    real(8), intent(in):: epot, fmax
end subroutine init_mybfgs
! ./src/optimizer_bfgs_minhocao.F90 :
subroutine geopt_init()
end subroutine geopt_init
subroutine GEOPT_RBFGS_MHM(parini,parres,latvec_in,xred_in,fcart_in,strten_in,etot_in,iprec,counter,folder)
 use mod_parini, only: typ_parini
type(typ_parini), intent(in):: parini
type(typ_parini), intent(inout):: parres
integer:: iprec,iat,i,istr
real(8):: latvec_in(3,3),xred_in(3,parini%nat),fcart_in(3,parini%nat),etot_in,strten_in(6),counter,flat(9)
character(40):: folder
end subroutine geopt_rbfgs_mhm
subroutine bfgs_driver_atoms(parini,parres,latvec_in,xred_in,fcart_in,strten_in,etot_in,iprec,counter,fmax_tol,folder)
    use mod_parini, only: typ_parini
    type(typ_parini), intent(in):: parini
    type(typ_parini), intent(inout):: parres
    real(8) :: latvec_in(3,3),xred_in(3,parini%nat),fcart_in(3,parini%nat),etot_in,strten_in(6),enthalpy
    real(8), intent(inout) :: counter
    real(8), dimension(3*parini%nat) :: rxyz
    real(8), dimension(3*parini%nat) :: fxyz
    real(8) :: fmax,fmax_at,fmax_lat,fmax_tol,en0000,betax
    character(len=40) :: comment,filename,coord,folder
    integer ::  nwork,iprec
end subroutine bfgs_driver_atoms
subroutine bfgs_driver_lattice(parini,parres,latvec_in,xred_in,fcart_in,strten_in,etot_in,iprec,counter,fmax_tol,folder)
    use mod_parini, only: typ_parini
    type(typ_parini), intent(in):: parini
    type(typ_parini), intent(inout):: parres
    real(8) :: latvec_in(3,3),xred_in(3,parini%nat),fcart_in(3,parini%nat),etot_in,strten_in(6),enthalpy,latvec(9)
    real(8), intent(inout) :: counter
    real(8) :: fmax,fmax_at,fmax_lat,fmax_tol,en0000
    character(len=40) :: comment, filename,coord,folder
    integer ::  nwork,iprec
end subroutine bfgs_driver_lattice
subroutine init_parameters(r0,fc)
    real(kind=8) :: r0(4,4),fc(4,4)
end subroutine init_parameters
subroutine pseudohess(nat,rat,nbond,indbond1,indbond2,sprcons,xl0,hess)
    integer :: nat,nbond,indbond1(nbond),indbond2(nbond)
    real(kind=8) :: rat(3,nat),sprcons(nbond),xl0(nbond),hess(3*nat,3*nat)
end subroutine pseudohess
subroutine bfgs_reza(nat,nr,x,epot,f,nwork,work,alphax_at,alphax_lat,fmax,fmax_at,fmax_lat,counter,coord)
    integer :: nat,nr,nwork,mf,my,ms,nrsqtwo,iw1,iw2,iw3,iw4,info,i,j,l,mx
    integer :: counter
    real(kind=8) :: x(nr),f(nr),epot,work(nwork),alphax_at,alphax_lat,alphax
    real(kind=8) :: DDOT,tt1,tt2,de,fnrm,fmax,beta,beta_lat,fmax_at,fmax_lat
    character(40):: coord
end subroutine bfgs_reza
subroutine lbfgs_driver_lattice(parini,parres,latvec_in,xred_in,fcart_in,strten_in,etot_in,iprec,counter,fail,fmax_tol,folder)
 use mod_parini, only: typ_parini
  type(typ_parini), intent(in):: parini
  type(typ_parini), intent(inout):: parres
  real(8) :: latvec_in(3*3),xred_in(3,parini%nat),fcart_in(3,parini%nat),etot_in,strten_in(6),enthalpy,enthalpyprev
  real(8), intent(inout) :: counter
  logical, intent(out) :: fail
  real(8) :: fmax,fmax_lat,fmax_at,fmax_tol,latvec_write(3*3),pressure,strtarget(6),dstr(6),de,str_matrix(3,3),vol
  real(8) :: strten_write(6),fcart_write(3,parini%nat),etot_write
  integer :: check,istr,iexit,iprec
  character(40):: filename,folder
end subroutine lbfgs_driver_lattice
subroutine atomic_copymoving_forward(nat,n,x,nr,xa)
    integer :: n,nr,i,iat,ixyz,ir,nat
    real(kind=8) :: x(n),xa(nr)
end subroutine atomic_copymoving_forward
subroutine atomic_copymoving_backward(nat,nr,xa,n,x)
    integer :: n,nr,i,iat,ixyz,ir,nat
    real(kind=8) :: x(n),xa(nr)
end subroutine atomic_copymoving_backward
subroutine get_BFGS_forces_PR(parini,parres,pos_all,force_all,enthalpy,getwfk,iprec,latvec_in,xred_in,etot_in,fcart_in,strten_in)
use mod_parini, only: typ_parini
type(typ_parini), intent(in):: parini
type(typ_parini), intent(inout):: parres
integer:: iprec,iat
real(8):: pos_all(3*parini%nat+9)
real(8):: force_all(3*parini%nat+9)
real(8):: enthalpy,pressure
real(8):: xred_in(3,parini%nat),fcart_in(3,parini%nat),strten_in(6),etot_in,latvec_in(3,3)
logical:: getwfk
end subroutine get_bfgs_forces_pr
subroutine getvol_strain(strain,latvec0,vol)
real(8), dimension(3,3):: latvec0,latvec,strain,unitmat
real(8):: vol
end subroutine getvol_strain
subroutine  get_BFGS_forces_strainlatt(parini,parres,pos_all,force_all,enthalpy,getwfk,iprec,latvec0,&
           &lattdeg,latvec_in,xred_in,etot_in,fcart_in,strten_in)
use mod_parini, only: typ_parini
type(typ_parini), intent(in):: parini
type(typ_parini), intent(inout):: parres
integer:: iprec,iat,lattdeg
real(8):: pos_all(3*parini%nat+9)
real(8):: force_all(3*parini%nat+9)
real(8):: enthalpy,pressure,vol,unitmat(3,3)
real(8):: xred_in(3,parini%nat),fcart_in(3,parini%nat),strten_in(6),etot_in,latvec_in(3,3),transformed(3,3),transformed_inv(3,3)
real(8):: str_matrix(3,3),flat(3,3),pressure_mat(3,3),tmplat(3,3),sigma(3,3),crossp(3),stressvol(3,3),latvec0(3,3)
logical:: getwfk
end subroutine get_bfgs_forces_strainlatt
subroutine correct_hessin(hess,hessin,latvec,ndim,hessupdate,lattdeg)
integer:: ndim,LWORK,info,i,j,hessupdate,lattdeg
real(8):: hessin(ndim,ndim),hess(ndim,ndim),hess_tmp(ndim,ndim),dmat(ndim,ndim),latvec(3,3)
end subroutine correct_hessin
SUBROUTINE unit_matrix(mat,ndim)
real(8),DIMENSION(ndim,ndim), INTENT(INOUT) :: mat
integer:: ndim
end subroutine unit_matrix
subroutine get_BFGS_forces_max(parini,parres,pos_all,force_all,enthalpy,getwfk,iprec,latvec_in,xred_in,etot_in,fcart_in,strten_in)
use mod_parini, only: typ_parini
type(typ_parini), intent(in):: parini
type(typ_parini), intent(inout):: parres
integer:: iprec,iat
real(8):: pos_all(3*parini%nat+9)
real(8):: force_all(3*parini%nat+9)
real(8):: enthalpy,pressure,vol
real(8):: xred_in(3,parini%nat),fcart_in(3,parini%nat),strten_in(6),etot_in,latvec_in(3,3),transformed(3,3),transformed_inv(3,3)
logical:: getwfk
end subroutine get_bfgs_forces_max
subroutine get_BFGS_forces_atom(parini,parres,pos,force,latvec,enthalpy,getwfk,iprec,latvec_in,xred_in,etot_in,fcart_in,strten_in)
use mod_parini, only: typ_parini
type(typ_parini), intent(in):: parini
type(typ_parini), intent(inout):: parres
integer:: iprec,iat
real(8):: pos(3*parini%nat),latvec(3,3)
real(8):: force(3*parini%nat)
real(8):: enthalpy,pressure,vol
real(8):: xred_in(3,parini%nat),fcart_in(3*parini%nat),strten_in(6),etot_in,latvec_in(3,3),transformed(3,3),transformed_inv(3,3)
logical:: getwfk
end subroutine get_bfgs_forces_atom
subroutine get_BFGS_forces_lattice(parini,parres,pos,force,latvec,enthalpy,getwfk,iprec,latvec_in,xred_in,etot_in,fcart_in,strten_in)
 use mod_parini, only: typ_parini
type(typ_parini), intent(in):: parini
type(typ_parini), intent(inout):: parres
integer:: iprec,iat
real(8):: pos(3,parini%nat),latvec(3,3)
real(8):: force(9)
real(8):: enthalpy,pressure,vol
real(8):: xred_in(3,parini%nat),fcart_in(3*parini%nat),strten_in(6),etot_in,latvec_in(3,3),transformed(3,3),transformed_inv(3,3)
logical:: getwfk
end subroutine get_bfgs_forces_lattice
subroutine sd_minhocao(nat,nr,x,epot,f,betax,betax_lat,fmax,iter)
integer:: nr,i,iter,nat
real(8):: x(nr),f(nr),epot,betax,betax_lat,fmax
end subroutine sd_minhocao
        subroutine stress_volume(latvec,vol,pressure,stressvol)
        real(8):: latvec(3,3),vol,stressvol(3,3),inv_latvec(3,3),pressure
end subroutine stress_volume
subroutine get_fmax(parini,fcart_in,strten_in,fmax,fmax_at,fmax_lat)
use mod_parini, only: typ_parini
type(typ_parini), intent(in):: parini
real(8):: fcart_in(3,parini%nat),strten_in(6),fmax,fmax_at,fmax_lat
end subroutine get_fmax
subroutine init_hessinv(parini,hessin,latvec,omega,b0,lattdeg) 
use mod_parini, only: typ_parini
type(typ_parini), intent(in):: parini
integer:: itype,iat,i,j,k,lattdeg
real(8):: omega,b0,hessin(3*parini%nat+9,3*parini%nat+9),diagat,avmass,diaglat
real(8):: amass(parini%nat),rcov,amass_u(parini%ntypat_global),vol
real(8),dimension(3,3):: diagat_lat,diagat_lat_inv,latvec,latvectrans
end subroutine init_hessinv
subroutine GEOPT_MBFGS_MHM(parini,parres,latvec_in,xred_in,fcart_in,strten_in,etot_in,iprec,counter,folder)
use mod_parini, only: typ_parini
type(typ_parini), intent(in):: parini
type(typ_parini), intent(inout):: parres
REAL(8) :: fret, counter
REAL(8), INTENT(INOUT) :: xred_in(3*parini%nat),latvec_in(9),fcart_in(3*parini%nat),strten_in(6),etot_in
INTEGER, PARAMETER :: ITMAX=4000
REAL(8), PARAMETER :: STPMX=1.0d0,EPS=epsilon(xred_in),TOLX=4.0d0*EPS
REAL(8):: dg(3*parini%nat+9),g(3*parini%nat+9),hdg(3*parini%nat+9),pnew(3*parini%nat+9),xi(3*parini%nat+9),p(3*parini%nat+9)
REAL(8):: tp(3*parini%nat+9),tg(3*parini%nat+9),dvin(3*parini%nat+9),vout(3*parini%nat+9),vout_prev(3*parini%nat+9)
REAL(8):: vin_min(3*parini%nat+9),vin(3*parini%nat+9),vel_in(3*parini%nat),vel_lat_in(9)
REAL(8):: vout_min(3*parini%nat+9),dedv_min(3*parini%nat+9)
REAL(8), DIMENSION(3*parini%nat+9,3*parini%nat+9) :: hessin,hessin0,hess_tmp,hess,hessin_dsyev
REAL(8) :: alpha_pl,dlatvec(9),dxred(3*parini%nat)
INTEGER :: choice,status,sumstatus,iprec,iexit,lattdeg,hessupdate
REAL(8) :: rxyz0(3*parini%nat),eval(3*parini%nat+9),fmax,fmax_at,fmax_lat,pressure
character(40)::filename,folder
end subroutine geopt_mbfgs_mhm
subroutine GEOPT_MBFGS_MHM_OLD(parini,parres,latvec_in,xred_in,fcart_in,strten_in,etot_in,iprec,counter,folder)
use mod_parini, only: typ_parini
type(typ_parini), intent(in):: parini
type(typ_parini), intent(inout):: parres
REAL(8) :: fret, counter
REAL(8), INTENT(INOUT) :: xred_in(3*parini%nat),latvec_in(9),fcart_in(3*parini%nat),strten_in(6),etot_in
INTEGER, PARAMETER :: ITMAX=4000
REAL(8), PARAMETER :: STPMX=1.0d0,EPS=epsilon(xred_in),TOLX=4.0d0*EPS
REAL(8):: dg(3*parini%nat+9),g(3*parini%nat+9),hdg(3*parini%nat+9),pnew(3*parini%nat+9),xi(3*parini%nat+9),p(3*parini%nat+9)
REAL(8):: tp(3*parini%nat+9),tg(3*parini%nat+9),dvin(3*parini%nat+9),vout(3*parini%nat+9),vout_prev(3*parini%nat+9)
REAL(8):: vin_min(3*parini%nat+9),vin(3*parini%nat+9)
REAL(8):: vout_min(3*parini%nat+9),dedv_min(3*parini%nat+9)
REAL(8), DIMENSION(3*parini%nat+9,3*parini%nat+9) :: hessin,hessin0
INTEGER :: choice,status,sumstatus,iprec,iexit
REAL(8) :: latvec0(9),rxyz0(3*parini%nat),eval(3*parini%nat+9),fmax,fmax_at,fmax_lat,pressure
character(40)::filename,folder
end subroutine geopt_mbfgs_mhm_old
FUNCTION vabs(v) result(res)
real(8),dimension(:):: v
real(8):: res
end function vabs
function outerprod(a,b)
real(8),dimension(:),intent(in)::a,b
real(8),dimension(size(a),size(b))::outerprod
end function outerprod
FUNCTION assert_eq(n1,n2,n3,n4,string) result(res)
CHARACTER(LEN=*), INTENT(IN) :: string
INTEGER, INTENT(IN) :: n1,n2,n3,n4
INTEGER :: res
end function assert_eq
subroutine findmin(choice,dedv_1,dedv_2,dedv_predict,&
& d2edv2_1,d2edv2_2,d2edv2_predict,&
& etotal_1,etotal_2,etotal_predict,&
& lambda_1,lambda_2,lambda_predict,status)
 integer,intent(in) :: choice
 integer,intent(out) :: status
 real(8),intent(in) :: dedv_1,dedv_2,etotal_1,etotal_2,lambda_1,lambda_2
 real(8),intent(out) :: d2edv2_1,d2edv2_2,d2edv2_predict,dedv_predict
 real(8),intent(out) :: etotal_predict,lambda_predict
end subroutine findmin
! ./src/optimizer_bfgs_qe.F90 :
     FUNCTION dot_product_( vec1, vec2 ) result(dotprod)
       REAL(8), INTENT(IN) :: vec1(:), vec2(:)
       REAL(8)             :: dotprod
end function dot_product_
     FUNCTION external_product_( vec1, vec2 )
       REAL(8), INTENT(IN) :: vec1(:), vec2(:)
       REAL(8)             :: external_product_(SIZE( vec1 ))
end function external_product_
     FUNCTION norm( vec ) result(dnrm)
       REAL(8), INTENT(IN) :: vec(:)
       REAL(8)             :: dnrm
end function norm
     FUNCTION matrix_times_vector( mat, vec )
       REAL(8), INTENT(IN) :: vec(:)
       REAL(8), INTENT(IN) :: mat(:,:)
       REAL(8)             :: matrix_times_vector(SIZE( vec ))
       REAL(8)             :: aux(SIZE( vec ))
end function matrix_times_vector
     FUNCTION vector_times_matrix( vec, mat )
       REAL(8), INTENT(IN) :: vec(:)
       REAL(8), INTENT(IN) :: mat(:,:)
       REAL(8)             :: vector_times_matrix(SIZE( vec ))
       REAL(8)             :: aux(SIZE( vec ))
end function vector_times_matrix
     FUNCTION matrix( vec1, vec2 )
       REAL(8), INTENT(IN) :: vec1(:), vec2(:)
       REAL(8)             :: matrix(SIZE( vec1 ),SIZE( vec2 ))
       REAL(8)             :: aux(SIZE( vec1 ),SIZE( vec2 ))
end function matrix
     FUNCTION identity( dim ) result(iden)
       INTEGER, INTENT(IN) :: dim
       REAL(8)            :: iden(dim,dim)
end function identity
   SUBROUTINE bfgs( pos_in, h, energy, grad_in, fcell, fixion, scratch, stdout,&
                 energy_thr, grad_thr, cell_thr, energy_error, grad_error,     &
                 cell_error, istep, nstep, step_accepted, stop_bfgs, lmovecell)
      REAL(8),         INTENT(INOUT) :: pos_in(:)
      REAL(8),         INTENT(INOUT) :: h(3,3)
      REAL(8),         INTENT(INOUT) :: energy
      REAL(8),         INTENT(INOUT) :: grad_in(:)
      REAL(8),         INTENT(INOUT) :: fcell(3,3)
      INTEGER,          INTENT(IN)    :: fixion(:)
      CHARACTER(LEN=*), INTENT(IN)    :: scratch
      INTEGER,          INTENT(IN)    :: stdout
      REAL(8),         INTENT(IN)    :: energy_thr, grad_thr, cell_thr
      INTEGER,          INTENT(OUT)   :: istep
      INTEGER,          INTENT(IN)    :: nstep
      REAL(8),         INTENT(OUT)   :: energy_error, grad_error, cell_error
      LOGICAL,          INTENT(OUT)   :: step_accepted, stop_bfgs
      LOGICAL,          INTENT(IN)    :: lmovecell
end subroutine bfgs
SUBROUTINE gdiis_step()
 REAL(8), ALLOCATABLE :: res(:,:), overlap(:,:), work(:)
 INTEGER,  ALLOCATABLE :: iwork(:)
 INTEGER               :: k, k_m, info
 REAL(8)              :: gamma0
end subroutine gdiis_step
SUBROUTINE reset_bfgs( n )
INTEGER, INTENT(IN) :: n
end subroutine reset_bfgs
SUBROUTINE read_bfgs_file( pos, grad, fixion, energy, scratch, n, stdout )
REAL(8),         INTENT(INOUT) :: pos(:)
REAL(8),         INTENT(INOUT) :: grad(:)
INTEGER,          INTENT(IN)    :: fixion(:)
CHARACTER(LEN=*), INTENT(IN)    :: scratch
INTEGER,          INTENT(IN)    :: n
INTEGER,          INTENT(IN)    :: stdout
REAL(8),         INTENT(INOUT) :: energy
end subroutine read_bfgs_file
SUBROUTINE write_bfgs_file( pos, energy, grad, scratch, n)
      INTEGER,         INTENT(IN) :: n
      REAL(8),         INTENT(IN) :: pos(:)
      REAL(8),         INTENT(IN) :: energy
      REAL(8),         INTENT(IN) :: grad(:)
      CHARACTER(LEN=*), INTENT(IN) :: scratch
end subroutine write_bfgs_file
   SUBROUTINE update_inverse_hessian( pos, grad, n, stdout )
      REAL(8), INTENT(IN)  :: pos(:)
      REAL(8), INTENT(IN)  :: grad(:)
      INTEGER,  INTENT(IN)  :: n
      INTEGER,  INTENT(IN)  :: stdout
end subroutine update_inverse_hessian
   SUBROUTINE check_wolfe_conditions( lwolfe, energy, grad )
      REAL(8), INTENT(IN)  :: energy
      REAL(8), INTENT(IN)  :: grad(:)
      LOGICAL,  INTENT(OUT) :: lwolfe
end subroutine check_wolfe_conditions
   FUNCTION energy_wolfe_condition ( energy ) result(res)
      REAL(8), INTENT(IN)  :: energy
      LOGICAL:: res
end function energy_wolfe_condition
   FUNCTION gradient_wolfe_condition ( grad ) result(res)
      REAL(8), INTENT(IN)  :: grad(:)
      LOGICAL:: res
end function gradient_wolfe_condition
   SUBROUTINE compute_trust_radius( lwolfe, energy, grad, n, stdout )
      LOGICAL,  INTENT(IN)  :: lwolfe
      REAL(8), INTENT(IN)  :: energy
      REAL(8), INTENT(IN)  :: grad(:)
      INTEGER,  INTENT(IN)  :: n
      INTEGER,  INTENT(IN)  :: stdout
end subroutine compute_trust_radius
   FUNCTION scnorm1( vect ) result(res)
      REAL(8), INTENT(IN) :: vect(:)
      REAL(8):: res
end function scnorm1
   FUNCTION scnorm( vect ) result(res)
      REAL(8), INTENT(IN) :: vect(:)
      REAL(8):: res
end function scnorm
   SUBROUTINE terminate_bfgs( energy, energy_thr, grad_thr, cell_thr, &
                              lmovecell, stdout, scratch )
      REAL(8),         INTENT(IN) :: energy, energy_thr, grad_thr, cell_thr
      LOGICAL,          INTENT(IN) :: lmovecell
      INTEGER,          INTENT(IN) :: stdout
      CHARACTER(LEN=*), INTENT(IN) :: scratch
end subroutine terminate_bfgs
subroutine GEOPT_qbfgs(parini,parres,latvec_in,xred_in,fcart_in,strten_in,etot_in,iprec,counter,folder)
  use mod_parini, only: typ_parini
  type(typ_parini), intent(in):: parini
  type(typ_parini), intent(inout):: parres
  real(8):: latvec_in(3,3),xred_in(3,parini%nat),fcart_in(3,parini%nat),strten_in(6),etot_in,counter,xred(3,parini%nat),fcart(3,parini%nat),latvec(3,3)
  integer:: iprec
  character(40):: folder
end subroutine geopt_qbfgs
subroutine recips (a1, a2, a3, b1, b2, b3)
  real(8) :: a1 (3), a2 (3), a3 (3), b1 (3), b2 (3), b3 (3)
end subroutine recips
subroutine cryst_to_cart (nvec, vec, trmat, iflag)
  integer, intent(in) :: nvec, iflag
  real(8), intent(in) :: trmat (3, 3)
  real(8), intent(inout) :: vec (3, nvec)
end subroutine cryst_to_cart
  subroutine cell_force( fcell, ainv, stress, omega, press)!, wmassIN )
    REAL(8), intent(out) :: fcell(3,3)
    REAL(8), intent(in) :: stress(3,3), ainv(3,3)
    REAL(8), intent(in) :: omega, press
end subroutine cell_force
subroutine invmat (n, a, a_inv, da)
  integer :: n
  real(8), DIMENSION (n,n) :: a, a_inv
  real(8) :: da
  integer :: info, lda, lwork, ipiv (n)
  real(8) :: work (n) 
end subroutine invmat
subroutine qe_volume (alat, a1, a2, a3, omega)
  real(8) :: alat, a1 (3), a2 (3), a3 (3), omega
end subroutine qe_volume
! ./src/optimizer_cg.F90 :
subroutine cgminimum(iproc,n,nr,x,f,epot,paropt,nwork,work)
    use mod_opt, only: typ_paropt, frmt_base
    type(typ_paropt):: paropt
    integer, intent(in):: iproc, n, nr, nwork
    real(8):: x(n), f(n), epot, work(nwork)
    character(46), parameter:: frmt='('//frmt_base//',e12.4,a)'
end subroutine cgminimum
subroutine init_cgminimum(paropt,n,nr,f,nwork,work,epot,fnrm)
    use mod_opt, only: typ_paropt
    type(typ_paropt), intent(inout):: paropt
    integer, intent(in):: n, nr, nwork
    real(8), intent(in):: f(n), epot, fnrm
    real(8), intent(inout):: work(nwork)
end subroutine init_cgminimum
! ./src/optimizer_dfp.F90 :
subroutine mydfp(nr,x,epot,f,nwork,work,paropt)
    use mod_opt, only: typ_paropt
    integer::nr,nwork,mf,my,ms,nrsq,iw1,iw2,iw3,info,i,j,l,mx
    real(8)::x(nr),f(nr),epot,work(nwork)
    type(typ_paropt)::paropt
end subroutine mydfp
! ./src/optimizer_diis.F90 :
subroutine diisminimum(n,nr,x,epot,f,paropt,nwork,work)
    use mod_opt, only: typ_paropt
    integer:: n, nr, nwork, info, id, jd
    real(8):: x(n), f(n), epot, work(nwork), fnrm, dnrm2, ddot, fmax
    type(typ_paropt):: paropt
end subroutine diisminimum
! ./src/optimizer_drivers.F90 :
subroutine minimize(parini,iproc,atoms,paropt)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, typ_file_info, update_ratp, update_rat
    use mod_opt, only: typ_paropt
    type(typ_parini), intent(in):: parini
    integer, intent(in):: iproc
    type(typ_atoms), intent(inout):: atoms
    type(typ_paropt), intent(inout):: paropt
end subroutine minimize
subroutine test_convergence(n,f,paropt)
    use mod_opt, only: typ_paropt
    integer, intent(in):: n
    real(8), intent(in):: f(n)
    type(typ_paropt), intent(inout):: paropt
end subroutine test_convergence
subroutine x_to_xr(n,x,f,bemoved,nr,xr,fr)
    integer, intent(in):: n, nr
    real(8), intent(in):: x(n), f(n)
    logical, intent(in):: bemoved(3,n/3)
    real(8), intent(inout):: xr(nr), fr(nr)
end subroutine x_to_xr
subroutine xr_to_x(nr,xr,n,bemoved,x)
    integer, intent(in):: n, nr
    real(8), intent(in):: xr(nr)
    logical, intent(in):: bemoved(3,n/3)
    real(8), intent(inout):: x(n)
end subroutine xr_to_x
subroutine report_param(paropt)
    use mod_opt, only: typ_paropt
    type(typ_paropt), intent(inout):: paropt
end subroutine report_param
subroutine initminimize(paropt)
    use mod_opt, only: typ_paropt
    type(typ_paropt), intent(inout):: paropt
end subroutine initminimize
subroutine finalminimize(paropt)
    use mod_opt, only: typ_paropt
    type(typ_paropt), intent(inout):: paropt
end subroutine finalminimize
! ./src/optimizer_drivers_vc.F90 :
subroutine vc_minimize(parini,iproc,atoms,paropt)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, typ_file_info
    use mod_opt, only: typ_paropt
    type(typ_parini), intent(in):: parini
    integer, intent(in):: iproc !, n, nr
    type(typ_atoms), intent(inout):: atoms
    type(typ_paropt), intent(inout):: paropt
end subroutine vc_minimize
subroutine vc_test_convergence(n,f,paropt)
    use mod_opt, only: typ_paropt
    integer, intent(in):: n
    real(8), intent(in):: f(n)
    type(typ_paropt), intent(inout):: paropt
end subroutine vc_test_convergence
subroutine vc_x_to_xr(atoms,nr,xr,fr)
    use mod_atoms, only: typ_atoms, get_rat
    type(typ_atoms), intent(in):: atoms
    integer, intent(in):: nr
    real(8), intent(inout):: xr(nr), fr(nr)
end subroutine vc_x_to_xr
subroutine vc_xr_to_x(nr,xr,atoms)
    use mod_atoms, only: typ_atoms, update_rat
    integer, intent(in):: nr
    real(8), intent(in):: xr(nr)
    type(typ_atoms), intent(inout):: atoms
end subroutine vc_xr_to_x
subroutine vc_report_param(paropt)
    use mod_opt, only: typ_paropt
    type(typ_paropt), intent(inout):: paropt
end subroutine vc_report_param
! ./src/optimizer_fire.F90 :
subroutine fire(parini,iproc,n,x,epot,f,work,paropt)
    use mod_parini, only: typ_parini
    use mod_opt, only: typ_paropt, frmt_base
    type(typ_parini), intent(in):: parini
    integer, intent(in):: n, iproc
    real(8), intent(inout):: x(n), epot, f(n)
    real(8), intent(inout):: work(3*n) !1:n velocities, n+1:2*n previous force
    type(typ_paropt), intent(inout):: paropt
    character(59), parameter:: frmt='('//frmt_base//',3es12.4,i4,1es12.4,a)'
end subroutine fire
subroutine init_fire(n,f,epot,work,paropt)
    use mod_opt, only: typ_paropt
    integer, intent(in):: n
    real(8), intent(in):: epot, f(n)
    real(8), intent(inout):: work(3*n)
    type(typ_paropt), intent(inout):: paropt
end subroutine init_fire
! ./src/optimizer_fire_minhocao.F90 :
subroutine GEOPT_FIRE_MHM(parini,parres,latvec_in,xred_in,fcart_in,strten_in,vel_in,vel_lat_in,vvol_in,etot_in,iprec,counter,folder)
 use mod_parini, only: typ_parini
 type(typ_parini), intent(in):: parini
 type(typ_parini), intent(inout):: parres
 real(8) :: latvec_in(3,3),xred_in(3,parini%nat),vel_in(3,parini%nat),fcart_in(3,parini%nat),etot_in,strten_in(6),vel_lat_in(3,3),vvol_in
 real(8):: amass(parini%nat)
 real(8),dimension(3,parini%nat):: xcart
 real(8),dimension(3,parini%nat):: fposcur
 real(8),dimension(3,parini%nat):: fpospred
 real(8),dimension(3,parini%nat):: vpospred
 real(8),dimension(3,parini%nat):: poscur
 real(8),dimension(3,parini%nat):: vxyz
 real(8),dimension(3,parini%nat):: vposcur
 real(8),dimension(3,parini%nat):: pospred
 real(8),dimension(3,parini%nat):: accposcur
 real(8),dimension(3,parini%nat):: accpospred
 real(8),dimension(3,parini%nat):: accposprev
 real(8),dimension(3,parini%nat):: dxred
 real(8):: counter 
 integer:: iprec
 character(40)::filename,folder
 real(8):: alpha,P,P_at,P_lat,fmax,fmax_at,fmax_lat,fall(3,parini%nat+3),fallnorm,vall(3,parini%nat+3),vallnorm
end subroutine geopt_fire_mhm
subroutine acceleration_fire(pressure,accpos,acclat,accvol,vpos,vlat,vvol,strten,fcart,latvec,amass,latmass,f0inv,md_type,nat) 
integer:: iat,i,j,md_type,nat
real(8),dimension(3,nat):: accpos,vpos,fcart,fpos
real(8),dimension(3,3)  :: acclat,vlat,latvec,tmplat,pressure,a,velmat,sigma,lattrans,latdottrans,gdot,g,ginv,gtot,str_matrix
real(8),dimension(3,3)  :: term1,term2,term3,term4,term5,term5_1,term5_2,sigmatrans,f0inv
real(8):: amass(nat),latmass,crossp(3),strten(6),vol,vpostmp(3),volvel,trace3
real(8):: accvol,vvol,vol_1_3
end subroutine acceleration_fire
! ./src/optimizer_gmdfire.F90 :
subroutine gmdfire(nr,x,epot,f,work,paropt)
    use mod_opt, only: typ_paropt
    integer:: nr
    real(8):: x(nr), epot, f(nr), de, DDOT, fnrm, fmax, vnrm, dt, p
    real(8):: work(5*nr) !1:nr velocities, nr+1:2*nr previous force
    type(typ_paropt)::paropt
end subroutine gmdfire
! ./src/optimizer_nlbfgs.F90 :
      SUBROUTINE NLBFGS(N,M,X,F,G,DIAG,W,paropt)
      use mod_opt, only: typ_paropt
      INTEGER N,M
      real(8):: X(N),G(N),DIAG(N),W(N*(2*M+1)+2*M)
      real(8):: F,fmax
      type(typ_paropt):: paropt
end subroutine nlbfgs
      SUBROUTINE NLB1(ITER,NFUN,GNORM,N,M,X,F,G,STP,FINISH,paropt)
      use mod_opt, only: typ_paropt
      type(typ_paropt):: paropt
      INTEGER ITER,NFUN,N,M
      LOGICAL FINISH
end subroutine nlb1
      SUBROUTINE NMCSRCH(N,X,F,G,S,STP,FTOL,MAXFEV,INFO,NFEV,WA,paropt)
      use mod_opt, only: typ_paropt
      type(typ_paropt):: paropt
      INTEGER N,MAXFEV,INFO,NFEV
end subroutine nmcsrch
      SUBROUTINE NMCSTEP(STX,FX,DX,STY,FY,DY,STP,FP,DP,BRACKT,STPMIN,STPMAX,INFO)
      INTEGER INFO
      LOGICAL BRACKT,BOUND
end subroutine nmcstep
! ./src/optimizer_sd.F90 :
subroutine sdminimum(parini,iproc,nr,x,f,epot,paropt,nwork,work)
    use mod_parini, only: typ_parini
    use mod_opt, only: typ_paropt, frmt_base
    type(typ_parini), intent(in):: parini
    integer, intent(in):: iproc, nr, nwork
    real(8), intent(inout):: x(nr), f(nr), work(nwork)
    real(8), intent(in):: epot
    type(typ_paropt), intent(inout):: paropt
    character(52), parameter:: frmt='('//frmt_base//',e12.4,i5,l2,a)'
end subroutine sdminimum
subroutine init_sdminimum(paropt,nr,x,nwork,work)
    use mod_opt, only: typ_paropt
    type(typ_paropt), intent(inout):: paropt
    integer, intent(in):: nr, nwork
    real(8), intent(in):: x(nr)
    real(8), intent(inout):: work(nwork)
end subroutine init_sdminimum
subroutine what_is_condition_of_feedback(paropt,de1,df1,feedbackcondition)
    use mod_opt, only: typ_paropt
    type(typ_paropt), intent(in):: paropt
    real(8), intent(in):: de1, df1
    logical, intent(out):: feedbackcondition
end subroutine what_is_condition_of_feedback
subroutine test_saturation(paropt,de1,de2,df2,fnrm)
    use mod_opt, only: typ_paropt
    type(typ_paropt), intent(inout):: paropt
    real(8), intent(in):: de1, de2, df2, fnrm !, fnrmitm1
end subroutine test_saturation
subroutine final_sdminimum(paropt)
    use mod_opt, only: typ_paropt
    type(typ_paropt), intent(inout):: paropt
end subroutine final_sdminimum
! ./src/optimizer_sd_minhocao.F90 :
subroutine GEOPT_SD(parini,parres,latvec_in,xred_in,fcart_in,strten_in,etot_in,iprec,counter,folder)
   use mod_parini, only: typ_parini
   type(typ_parini), intent(in):: parini
   type(typ_parini), intent(inout):: parres
   character(len=*), parameter :: subname='sqnm'
   integer :: it,i,iat,l,j,idim,jdim,ihist,icheck !<counter variables
   integer:: ncount_cluster_x,iexit,iprec,lattdeg
   real(8):: latvec_in(3,3),xred_in(3,parini%nat),fcart_in(3,parini%nat),strten_in(6),etot_in,counter,pressure,latvec0(3,3),enthalpy_old
   character(40)::filename,folder
end subroutine geopt_sd
! ./src/optimizer_simplex.F90 :
subroutine simplex(vertices,fval,step,ndim,ftol,functn,iter)
    integer, intent(in):: ndim
    real(8), intent(in):: ftol, step
    real(8), intent(inout):: vertices(ndim,ndim+1), fval(ndim+1)
    integer, intent(out):: iter
    real(8), parameter:: alpha=1.d0 !reflection coefficient, a positive value
    real(8), parameter:: beta=0.5d0 !contraction coefficient, greater than one
    real(8), parameter:: gama=2.d0 !expansion coefficient, lies between zero and one
    integer, parameter:: nmax=20, itmax=50000
end subroutine simplex
    subroutine functn(n,p,func)
        integer, intent(in)  :: n
        real(8), intent(in)  :: p(n)
        real(8), intent(out) :: func
end subroutine functn
! ./src/optimizer_sqnm.F90 :
subroutine sqnm(parini,atoms,paropt,count_sqnm,fail)
   use mod_parini, only: typ_parini
   use mod_atoms, only: typ_atoms, set_rat, get_rat
   use mod_opt, only: typ_paropt, frmt_base
   type(typ_parini), intent(in):: parini
   type(typ_atoms), intent(inout):: atoms
   type(typ_paropt), intent(inout):: paropt
   real(8), intent(inout):: count_sqnm
   logical, intent(out):: fail
   character(len=*), parameter :: subname='sqnm'
   character(41), parameter:: frmt='('//frmt_base//',i5)'
   integer :: nat    !< number of atoms
   real(8) :: trustr !< a single atoms is not allowed to be dsiplaced more than by trustr
end subroutine sqnm
subroutine minenergyandforces_alborz(parini,iproc,nproc,eeval,imode,atoms,nat,rat,fat,fstretch,fxyzraw,epot,iconnect,nbond_,wold,alpha_stretch0,alpha_stretch,infocode)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, set_rat
    type(typ_parini), intent(in):: parini
    integer, intent(in)           :: iproc,nproc,imode
    type(typ_atoms), intent(inout):: atoms
    integer, intent(in)           :: nat
    integer, intent(in)           :: nbond_
    integer, intent(in)           :: iconnect(2,nbond_)
    real(8),intent(inout)        :: rat(3,nat)
    real(8),intent(out)          :: fxyzraw(3,nat)
    real(8),intent(inout)        :: fat(3,nat)
    real(8),intent(out)          :: fstretch(3,nat)
    real(8), intent(inout)       :: wold(nbond_)
    real(8), intent(in)          :: alpha_stretch0
    real(8), intent(inout)       :: alpha_stretch
    real(8), intent(inout)       :: epot
    logical, intent(in)           :: eeval
    integer,intent(out) :: infocode
end subroutine minenergyandforces_alborz
subroutine give_rcov_sqnm(iproc,atoms,rcov)
    use mod_atoms, only: typ_atoms
    integer, intent(in) :: iproc
    type(typ_atoms), intent(in) :: atoms
    real(8), intent(out) :: rcov(atoms%nat)
end subroutine give_rcov_sqnm
subroutine getSubSpaceEvecEval(label,iproc,verbosity,nat,nhist,nhistx,ndim,cutoffratio,lwork,work,idx,rxyz,fxyz,aa,rr,ff,rrr,fff,eval,res,success)
    integer, intent(in) :: iproc,verbosity,nat,nhist,nhistx,lwork
    character(len=*), intent(in) :: label
    integer, intent(out) :: ndim
    integer, intent(in) :: idx(0:nhistx)
    real(8), intent(in) :: rxyz(3,nat,0:nhistx),fxyz(3,nat,0:nhistx)
    real(8), intent(out) :: aa(nhistx,nhistx),eval(nhistx)
    real(8), intent(out) :: work(lwork)
    real(8), intent(out) :: rr(3,nat,0:nhistx), ff(3,nat,0:nhistx)
    real(8), intent(out) :: rrr(3,nat,0:nhistx), fff(3,nat,0:nhistx)
    real(8), intent(out) :: res(nhistx)
    real(8), intent(in) :: cutoffratio
    logical, intent(out) :: success
    real(8) :: rnorm(nhistx)
end subroutine getsubspaceeveceval
subroutine modify_gradient(nat,ndim,rrr,eval,res,fxyz,alpha,dd)
    integer, intent(in) :: nat
    integer, intent(in) :: ndim
    real(8), intent(out) :: dd(3,nat)
    real(8), intent(in) :: fxyz(3,nat)
    real(8), intent(in) :: rrr(3,nat,ndim)
    real(8), intent(in) :: eval(ndim)
    real(8), intent(in) :: res(ndim)
    real(8), intent(in) :: alpha
    real(8) :: scpr(ndim)
end subroutine modify_gradient
subroutine findbonds(label,iproc,verbosity,atoms,rcov,nbond,iconnect)
    use mod_atoms, only: typ_atoms, get_rat
    integer, intent(in) :: iproc,verbosity
    character(len=*), intent(in) :: label
    type(typ_atoms), intent(in) :: atoms
    real(8), intent(in) :: rcov(atoms%nat)
    integer, intent(out) :: nbond
    integer, intent(out) :: iconnect(2,1000)
end subroutine findbonds
subroutine projectbond(nat,nbond,rat,fat,fstretch,iconnect,wold,alpha_stretch0,alpha_stretch)
    integer, intent(in) :: nat
    integer, intent(in) :: nbond
    real(8), intent(in) :: rat(3,nat)
    real(8), intent(inout) :: fat(3,nat)
    real(8), intent(inout) :: fstretch(3,nat)
    integer, intent(in) :: iconnect(2,nbond)
    real(8), intent(inout) :: wold(nbond)
    real(8), intent(in) :: alpha_stretch0
    real(8), intent(inout) :: alpha_stretch
    real(8) :: ss(nbond,nbond),w(nbond),vv(3,nat,nbond)
end subroutine projectbond
! ./src/optimizer_sqnm_minhocao.F90 :
subroutine GEOPT_sqnm(parini,parres,latvec_in,xred_in,fcart_in,strten_in,etot_in,iprec,counter,folder)
   use mod_parini, only: typ_parini
   type(typ_parini), intent(in):: parini
   type(typ_parini), intent(inout):: parres
   character(len=*), parameter :: subname='sqnm'
   integer :: it,i,iat,l,j,idim,jdim,ihist,icheck !<counter variables
   integer:: ncount_cluster_x,iexit,iprec
   real(8):: latvec_in(3,3),xred_in(3,parini%nat),fcart_in(3,parini%nat),strten_in(6),etot_in,counter,pressure
   character(40)::filename,folder
   real(8):: pos_tmp(3,parini%nat),latvec_old(3,3)
end subroutine geopt_sqnm
subroutine minenergyandforces(parini,parres,eeval,imode,nat,rat,rxyzraw,fat,fstretch,&
           fxyzraw,epot,alpha_stretch0,alpha_stretch,&
           latvec_in,xred_in,etot_in,fcart_in,strten_in,iprec)
    use mod_parini, only: typ_parini
    type(typ_parini), intent(in):: parini
    type(typ_parini), intent(inout):: parres
    integer, intent(in)           :: imode
    integer, intent(in)           :: nat
    real(8),intent(inout)        :: rat(3,nat+3)
    real(8),intent(out)          :: rxyzraw(3,nat+3)
    real(8),intent(out)          :: fxyzraw(3,nat+3)
    real(8),intent(inout)        :: fat(3,nat+3)
    real(8),intent(out)          :: fstretch(3,nat+3)
    real(8), intent(in)          :: alpha_stretch0
    real(8), intent(inout)       :: alpha_stretch
    real(8), intent(inout)       :: epot
    logical, intent(in)          :: eeval
    real(8):: rxyz(3,nat+3)
    real(8):: fxyz(3,nat+3)
    real(8):: force_all(3,nat+3)
    real(8):: latvec0(3,3),latvec_in(3,3),xred_in(3,nat),etot_in,fcart_in(3,nat),strten_in(6) 
    integer:: lattdeg=1,iprec
end subroutine minenergyandforces
subroutine sqnm_invhess(nat,h,metric,hessinv)
integer:: nat,info,i,j,k
real(8):: metric(3*(nat+3),3*(nat+3)),hessinv(3*(nat+3),3*(nat+3))
real(8):: h(3,3),hinv(3,3),g(3,3),ginv(3,3)
real(8):: hessinv_at(3*nat,3*nat),hessinv_lat(9,9)
real(8):: eval(3*(nat+3)),eval_at(3*nat),eval_lat(9)
end subroutine sqnm_invhess
! ./src/opt_mod.F90 :
! ./src/parini_mod.F90 :
! ./src/parser_all.F90 :
subroutine get_main_parameters(file_ini,parini)
    use mod_task, only: typ_file_ini
    use mod_parini, only: typ_parini
    type(typ_file_ini), intent(inout):: file_ini
    type(typ_parini), intent(inout):: parini
end subroutine get_main_parameters
subroutine get_minhopp_parameters(file_ini,parini)
    use mod_parini, only: typ_parini
    use mod_task, only: typ_file_ini
    type(typ_file_ini), intent(inout):: file_ini
    type(typ_parini), intent(inout):: parini
end subroutine get_minhopp_parameters
subroutine get_opt_param(file_ini,paropt)
    use mod_task, only: typ_file_ini
    use mod_opt, only: typ_paropt
    type(typ_file_ini), intent(inout):: file_ini
    type(typ_paropt), intent(inout):: paropt
end subroutine get_opt_param
subroutine get_geopt_parameters(file_ini,parini)
    use mod_parini, only: typ_parini
    use mod_task, only: typ_file_ini
    type(typ_file_ini), intent(inout):: file_ini
    type(typ_parini), intent(inout):: parini
end subroutine get_geopt_parameters
subroutine get_geopt_prec_parameters(file_ini,parini)
    use mod_parini, only: typ_parini
    use mod_task, only: typ_file_ini
    type(typ_file_ini), intent(inout):: file_ini
    type(typ_parini), intent(inout):: parini
end subroutine get_geopt_prec_parameters
subroutine get_saddle_1s_opt_parameters(file_ini,parini)
    use mod_parini, only: typ_parini
    use mod_task, only: typ_file_ini
    type(typ_file_ini), intent(inout):: file_ini
    type(typ_parini), intent(inout):: parini
end subroutine get_saddle_1s_opt_parameters
subroutine get_saddle_1s_parameters(file_ini,parini)
    use mod_parini, only: typ_parini
    use mod_task, only: typ_file_ini
    type(typ_file_ini), intent(inout):: file_ini
    type(typ_parini), intent(inout):: parini
end subroutine get_saddle_1s_parameters
subroutine get_potential_parameters(file_ini,parini)
    use mod_parini, only: typ_parini
    use mod_task, only: typ_file_ini
    type(typ_file_ini), intent(inout):: file_ini
    type(typ_parini), intent(inout):: parini
end subroutine get_potential_parameters
subroutine get_ann_parameters(file_ini,parini)
    use mod_parini, only: typ_parini
    use mod_task, only: typ_file_ini
    type(typ_file_ini), intent(inout):: file_ini
    type(typ_parini), intent(inout):: parini
end subroutine get_ann_parameters
subroutine get_dynamics_parameters(file_ini,parini)
    use mod_parini, only: typ_parini
    use mod_task, only: typ_file_ini
    type(typ_file_ini), intent(inout):: file_ini
    type(typ_parini), intent(inout):: parini
end subroutine get_dynamics_parameters
subroutine get_bader_parameters(file_ini,parini)
    use mod_parini, only: typ_parini
    use mod_task, only: typ_file_ini
    type(typ_file_ini), intent(inout):: file_ini
    type(typ_parini), intent(inout):: parini
end subroutine get_bader_parameters
subroutine get_genconf_parameters(file_ini,parini)
    use mod_parini, only: typ_parini
    use mod_task, only: typ_file_ini
    type(typ_file_ini), intent(inout):: file_ini
    type(typ_parini), intent(inout):: parini
end subroutine get_genconf_parameters
subroutine get_conf_comp_parameters(file_ini,parini)
    use mod_parini, only: typ_parini
    use mod_task, only: typ_file_ini
    type(typ_file_ini), intent(inout):: file_ini
    type(typ_parini), intent(inout):: parini
end subroutine get_conf_comp_parameters
subroutine get_testforces_parameters(file_ini,parini)
    use mod_parini, only: typ_parini
    use mod_task, only: typ_file_ini
    type(typ_file_ini), intent(inout):: file_ini
    type(typ_parini), intent(inout):: parini
end subroutine get_testforces_parameters
subroutine get_single_point_parameters(file_ini,parini)
    use mod_parini, only: typ_parini
    use mod_task, only: typ_file_ini
    type(typ_file_ini), intent(inout):: file_ini
    type(typ_parini), intent(inout):: parini
end subroutine get_single_point_parameters
subroutine get_ewald_parameters(file_ini,parini)
    use mod_parini, only: typ_parini
    use mod_task, only: typ_file_ini
    type(typ_file_ini), intent(inout):: file_ini
    type(typ_parini), intent(inout):: parini
end subroutine get_ewald_parameters
subroutine get_misc_parameters(file_ini,parini)
    use mod_parini, only: typ_parini
    use mod_task, only: typ_file_ini
    type(typ_file_ini), intent(inout):: file_ini
    type(typ_parini), intent(inout):: parini
end subroutine get_misc_parameters
! ./src/parser_core.F90 :
subroutine read_file_input(file_ini)
    use mod_task, only: typ_file_ini
    type(typ_file_ini), intent(inout):: file_ini
end subroutine read_file_input
subroutine get_header_location(file_ini,str_header)
    use mod_task, only: typ_file_ini
    type(typ_file_ini), intent(inout):: file_ini
    character(*), intent(in):: str_header
end subroutine get_header_location
subroutine split_line(file_ini)
    use mod_task, only: typ_file_ini
    type(typ_file_ini), intent(inout):: file_ini
end subroutine split_line
subroutine get_one_param(file_ini,var_name,int_var,real_var,char_var,char_line_var,log_var)
    use mod_task, only: typ_file_ini
    type(typ_file_ini), intent(inout):: file_ini
    character(*), intent(in):: var_name
    integer, optional, intent(out):: int_var
    real(8), optional, intent(out):: real_var
    character(*), optional, intent(out):: char_var
    character(*), optional, intent(out):: char_line_var
    logical, optional, intent(out):: log_var
end subroutine get_one_param
! ./src/parser_minhocao.F90 :
subroutine params_read(parini)
use mod_parini, only: typ_parini
type(typ_parini), intent(inout):: parini
end subroutine params_read
subroutine params_read_for_yaml(parini)
    use mod_parini, only: typ_parini
    type(typ_parini), intent(inout):: parini
end subroutine params_read_for_yaml
subroutine params_defaults(parini,mdmin_in,dtion_md_in,alpha_lat_in,alpha_at_in,read_poscur)
use mod_parini, only: typ_parini
type(typ_parini), intent(inout):: parini
integer:: mdmin_in,itype,i,j
real(8):: dtion_md_in,alpha_lat_in,alpha_at_in
logical:: read_poscur
end subroutine params_defaults
subroutine params_check(parini)
use mod_parini, only: typ_parini
type(typ_parini), intent(in):: parini
end subroutine params_check
subroutine params_echo(parini)
use mod_parini, only: typ_parini
type(typ_parini), intent(in):: parini
end subroutine params_echo
subroutine fp_assign(parini)
use mod_parini, only: typ_parini
type(typ_parini), intent(inout):: parini
end subroutine fp_assign
! ./src/parser_yaml.F90 :
subroutine yaml_get_parameters(parini)
    use mod_parini, only: typ_parini
    type(typ_parini), intent(inout):: parini
end subroutine yaml_get_parameters
subroutine yaml_get_main_parameters(parini)
    use mod_parini, only: typ_parini
    type(typ_parini), intent(inout):: parini
end subroutine yaml_get_main_parameters
subroutine yaml_get_minhopp_parameters(parini)
    use mod_parini, only: typ_parini
    type(typ_parini), intent(inout):: parini
end subroutine yaml_get_minhopp_parameters
subroutine yaml_get_opt_parameters(parini,paropt)
    use mod_parini, only: typ_parini
    use mod_opt, only: typ_paropt
    type(typ_parini), intent(in):: parini
    type(typ_paropt), intent(inout):: paropt
end subroutine yaml_get_opt_parameters
subroutine yaml_get_geopt_parameters(parini)
    use mod_parini, only: typ_parini
    type(typ_parini), intent(inout):: parini
end subroutine yaml_get_geopt_parameters
subroutine yaml_get_geopt_prec_parameters(parini)
    use mod_parini, only: typ_parini
    type(typ_parini), intent(inout):: parini
end subroutine yaml_get_geopt_prec_parameters
subroutine yaml_get_saddle_1s_opt_parameters(parini)
    use mod_parini, only: typ_parini
    type(typ_parini), intent(inout):: parini
end subroutine yaml_get_saddle_1s_opt_parameters
subroutine yaml_get_saddle_1s_parameters(parini)
    use mod_parini, only: typ_parini
    type(typ_parini), intent(inout):: parini
end subroutine yaml_get_saddle_1s_parameters
subroutine yaml_get_potential_parameters(parini)
    use mod_parini, only: typ_parini
    type(typ_parini), intent(inout):: parini
end subroutine yaml_get_potential_parameters
subroutine yaml_get_ann_parameters(parini)
    use mod_parini, only: typ_parini
    type(typ_parini), intent(inout):: parini
end subroutine yaml_get_ann_parameters
subroutine yaml_get_dynamics_parameters(parini)
    use mod_parini, only: typ_parini
    type(typ_parini), intent(inout):: parini
end subroutine yaml_get_dynamics_parameters
subroutine yaml_get_bader_parameters(parini)
    use mod_parini, only: typ_parini
    type(typ_parini), intent(inout):: parini
end subroutine yaml_get_bader_parameters
subroutine yaml_get_genconf_parameters(parini)
    use mod_parini, only: typ_parini
    type(typ_parini), intent(inout):: parini
end subroutine yaml_get_genconf_parameters
subroutine yaml_get_conf_comp_parameters(parini)
    use mod_parini, only: typ_parini
    type(typ_parini), intent(inout):: parini
end subroutine yaml_get_conf_comp_parameters
subroutine yaml_get_testforces_parameters(parini)
    use mod_parini, only: typ_parini
    type(typ_parini), intent(inout):: parini
end subroutine yaml_get_testforces_parameters
subroutine yaml_get_single_point_parameters(parini)
    use mod_parini, only: typ_parini
    type(typ_parini), intent(inout):: parini
end subroutine yaml_get_single_point_parameters
subroutine yaml_get_ewald_parameters(parini)
    use mod_parini, only: typ_parini
    type(typ_parini), intent(inout):: parini
end subroutine yaml_get_ewald_parameters
subroutine yaml_get_misc_parameters(parini)
    use mod_parini, only: typ_parini
    type(typ_parini), intent(inout):: parini
end subroutine yaml_get_misc_parameters
subroutine yaml_get_confinement_parameters(parini)
    use mod_parini, only: typ_parini
    type(typ_parini), intent(inout):: parini
end subroutine yaml_get_confinement_parameters
subroutine yaml_get_fingerprint_parameters(parini)
    use mod_parini, only: typ_parini
    type(typ_parini), intent(inout):: parini
end subroutine yaml_get_fingerprint_parameters
subroutine yaml_get_fit_elecpot_parameters(parini)
    use mod_parini, only: typ_parini
    type(typ_parini), intent(inout):: parini
end subroutine yaml_get_fit_elecpot_parameters
subroutine set_dict_parini_default(parini)
    use mod_parini, only: typ_parini
    type(typ_parini), intent(inout):: parini
end subroutine set_dict_parini_default
subroutine set_dict_parini_user(parini)
    use mod_parini, only: typ_parini
    type(typ_parini), intent(inout):: parini
    character(len=*), parameter:: fname="flame_in.yaml"
end subroutine set_dict_parini_user
subroutine check_nonoptional_parameters(parini)
    use mod_parini, only: typ_parini
    type(typ_parini), intent(in):: parini
end subroutine check_nonoptional_parameters
! ./src/pbc_distance.F90 :
subroutine pbc_distance0(latvec,xred_1,xred_2,distance2,dxyz)
real(8):: xred_1(3),xred_2(3),diff(3),distance2,latvec(3,3),dxyz(3)
end subroutine pbc_distance0
subroutine pbc_distance1(latvec,xred_1,xred_2,distance2)
real(8):: xred_1(3),xred_2(3),diff(3),distance2,latvec(3,3)
end subroutine pbc_distance1
subroutine pbc_distance2(latvec,xred_1,xcart_1,xred_2,xcart_2,distance2)
real(8):: xred_1(3),xred_2(3),diff(3),distance2,latvec(3,3),xcart_1(3),xcart_2(3),xcart_tmp(3),xcart_20(3),xcart_10(3)
end subroutine pbc_distance2
! ./src/phonon.F90 :
subroutine cal_hessian_4p(parini)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, typ_atoms_arr, typ_file_info, atom_copy_old
    type(typ_parini), intent(in):: parini
end subroutine cal_hessian_4p
subroutine projectout_rotation(atoms,hess,rlarge)
    use mod_atoms, only: typ_atoms, update_ratp
    type(typ_atoms), intent(inout):: atoms
    real(8), intent(inout):: hess(3*atoms%nat,3*atoms%nat)
    real(8), intent(in):: rlarge
end subroutine projectout_rotation
! ./src/plain_ewald.F90 :
subroutine plain_ewald(atoms,en)
    use mod_atoms, only: typ_atoms
    type(typ_atoms), intent(inout):: atoms
    real(8):: pi, q(1:atoms%nat)
    real(8):: k, k2, alpha, sum_en, en
end subroutine plain_ewald
subroutine structur_factor_comput(atoms,kx,ky,kz,s1,s2,sfactor_norm)
    use mod_atoms, only: typ_atoms
    type(typ_atoms), intent(inout):: atoms
    real(8), allocatable, intent(inout)::sfactor_norm(:,:,:) 
    integer, intent(inout):: kx, ky, kz
    real(8), allocatable:: s1(:,:,:), s2(:,:,:)
    real(8):: pi, q(1:atoms%nat)
end subroutine structur_factor_comput
! ./src/potential_ANN.F90 :
subroutine init_potential_ann(parini,atoms)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
end subroutine init_potential_ann
subroutine cal_potential_ann(parini,atoms)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, atom_deallocate_old, get_rat
    use mod_symfunc, only: typ_symfunc
    use mod_opt_ann, only: typ_opt_ann
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
end subroutine cal_potential_ann
subroutine final_potential_ann
end subroutine final_potential_ann
subroutine add_repulsive_potential(parini,atoms)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, set_rcov, update_ratp
    use mod_linked_lists, only: typ_linked_lists
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
end subroutine add_repulsive_potential
! ./src/potential_BigDFT.F90 :
subroutine init_potential_forces_bigdft(atoms_t)
    use mod_atoms, only: typ_atoms
    type(typ_atoms), intent(in):: atoms_t
end subroutine init_potential_forces_bigdft
subroutine final_potential_forces_bigdft
end subroutine final_potential_forces_bigdft
subroutine cal_potential_forces_bigdft(atoms)
    use mod_atoms, only: typ_atoms, update_ratp
    type(typ_atoms), intent(inout):: atoms
end subroutine cal_potential_forces_bigdft
subroutine writexyz_bigdft(filename,nat,rat,comment)
    integer:: nat,iat
    real(8)::rat(3,nat),x,y,z,cellx,celly,cellz
    character(*)::filename,comment
end subroutine writexyz_bigdft
subroutine get_output_bigdft(iproc,filename,nat,fat,epot,success)
    integer, intent(in):: iproc, nat
    character(*):: filename
    real(8):: fat(3,nat),epot
    logical:: success
end subroutine get_output_bigdft
! ./src/potential_BLJ_vc.F90 :
subroutine init_lennardjones_vc(nat,sat)
    integer, intent(in):: nat
    character(5), intent(in):: sat(nat)
    real(8), parameter:: alphalj=2.5d0
end subroutine init_lennardjones_vc
subroutine lennardjones_vc(iproc,nat,xred0,latvec,pressure,fxyz,celldv,stress,etot,enth)
    integer, intent(in):: iproc, nat
    real(8):: xred(3,nat),fxyz(3,nat),xred0(3,nat),dxyz(3),r1red(3),r2red(3),rcut2(2,2)
    real(8):: etot,latvec_x(3,3),rec_nec(3),enth,stressvol(3,3),vol
    real(8):: latvec(3,3),latvecinv(3,3),celldv(3,3),pressure,gradvol(3,3),stress(3,3),tmplat(3,3)
end subroutine lennardjones_vc
subroutine stress_volume_alborz(latvec,vol,pressure,stressvol)
    real(8):: latvec(3,3),vol,stressvol(3,3),inv_latvec(3,3),pressure
end subroutine stress_volume_alborz
subroutine cell_vol(nat,latvec,vol)
    integer:: nat
    real(8):: latvec(3,3),vol,a(3,3)
end subroutine cell_vol
! ./src/potential_confinement.F90 :
subroutine init_confinement_parser(parini)
use mod_parini, only: typ_parini
type(typ_parini), intent(inout):: parini
end subroutine init_confinement_parser
subroutine confinement_energy_forces(parini,nat,xred,latvec,energy,forces,strten)
use mod_parini, only: typ_parini
type(typ_parini), intent(in):: parini
integer:: nat,iconf,iat
real(8):: xred(3,nat),latvec(3,3),energy,forces(3,nat),dist,dist_av,nvec(3,3),point0(3),point(3)
real(8):: xcart(3,nat),tt,flat(3,3),xred_ppoint(3),str(3,3),strten(6),vol,fcart_all(3),ft(3)
end subroutine confinement_energy_forces
subroutine conf_latforce(latvec,conf_dim,xred_point,xred_ppoint,str)
real(8):: latvec(3,3),xred_point(3),xred_ppoint(3),str(3,3)
integer:: conf_dim
end subroutine conf_latforce
! ./src/potential_DFTB.F90 :
subroutine init_potential_forces_dftb(atoms)
    use mod_atoms, only: typ_atoms
    type(typ_atoms), intent(in):: atoms
end subroutine init_potential_forces_dftb
subroutine final_potential_forces_dftb
end subroutine final_potential_forces_dftb
subroutine cal_potential_forces_dftb(atoms)
    use mod_atoms, only: typ_atoms, set_typat
    type(typ_atoms), intent(inout):: atoms
end subroutine cal_potential_forces_dftb
subroutine writexyz_dftb(filename,atoms)
    use mod_atoms, only: typ_atoms, get_rat_iat
    character(*):: filename !,comment
    type(typ_atoms), intent(in):: atoms
end subroutine writexyz_dftb
subroutine get_output_dftb_alborz(filename,atoms,success)
    use mod_atoms, only: typ_atoms
    character(*), intent(in):: filename
    type(typ_atoms), intent(inout):: atoms
    logical, intent(out):: success
end subroutine get_output_dftb_alborz
! ./src/potential_FF.F90 :
subroutine init_potential_forces_ff(parini,atoms)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
end subroutine init_potential_forces_ff
subroutine cal_potential_forces_ff(parini,atoms)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
end subroutine cal_potential_forces_ff
subroutine final_potential_forces_ff(parini,atoms)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(in):: atoms
end subroutine final_potential_forces_ff
! ./src/potential_flame.F90 :
subroutine call_to_alborz_init(parini,nat)
    use mod_parini, only: typ_parini
    type(typ_parini), intent(in):: parini
    integer, intent(in):: nat
end subroutine call_to_alborz_init
subroutine call_to_alborz_get(boundcond,nat,latvec,xred,fcart,energy,strten)
    character(*), intent(in):: boundcond
    integer, intent(in):: nat
    real(8), intent(in):: latvec(3,3), xred(3,nat)
    real(8), intent(inout):: fcart(3,nat), energy, strten(6)
end subroutine call_to_alborz_get
subroutine call_to_alborz_final
end subroutine call_to_alborz_final
! ./src/potential_LJ.F90 :
subroutine init_lennardjones
end subroutine init_lennardjones
subroutine lennardjones(atoms)
    use mod_atoms, only: typ_atoms, update_ratp
    type(typ_atoms), intent(inout):: atoms
end subroutine lennardjones
! ./src/potential_LTB.F90 :
subroutine init_lenosky_tb(atoms_t)
    use mod_atoms, only: typ_atoms
    type(typ_atoms), intent(in):: atoms_t
end subroutine init_lenosky_tb
subroutine lenosky_tb(parini,atoms)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, update_ratp, update_rat
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
end subroutine lenosky_tb
! ./src/potential_main.F90 :
subroutine init_potential_forces(parini,atoms)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
end subroutine init_potential_forces
subroutine cal_potential_forces(parini,atoms)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, get_rat, update_ratp, set_rat
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
end subroutine cal_potential_forces
subroutine final_potential_forces(parini,atoms)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(in):: atoms
end subroutine final_potential_forces
subroutine remove_drift(atoms)
    use mod_atoms, only: typ_atoms
    type(typ_atoms), intent(inout):: atoms
end subroutine remove_drift
! ./src/potential_main_vc.F90 :
subroutine vc_init_potential_forces(atoms)
    use mod_atoms, only: typ_atoms
    type(typ_atoms), intent(inout):: atoms
end subroutine vc_init_potential_forces
subroutine cal_potential_forces_vc(iproc,nat,rat,cellvec,pressure,fat,celldv,stress,epot,enth)
    integer, intent(in):: iproc, nat
    real(8), intent(in):: rat(3,nat), cellvec(3,3), pressure
    real(8), intent(inout):: fat(3,nat), celldv(3,3), stress(3,3), epot, enth
end subroutine cal_potential_forces_vc
! ./src/potential_mod.F90 :
! ./src/potential_MPMD.F90 :
subroutine init_potential_forces_mpmd(atoms_t)
    use mod_atoms, only: typ_atoms
    type(typ_atoms), intent(in):: atoms_t
end subroutine init_potential_forces_mpmd
subroutine cal_potential_forces_mpmd(atoms)
    use mod_atoms, only: typ_atoms, update_ratp, update_rat
    type(typ_atoms), intent(inout):: atoms
end subroutine cal_potential_forces_mpmd
subroutine mpmd_init
end subroutine mpmd_init
! ./src/potential_NetSock.F90 :
  subroutine cal_potential_forces_netsock(atoms)
  use mod_potential, only: sock_socket, sock_inet, sock_port,sock_host,MSGLEN,sock_extra_string,reset
  use mod_atoms, only: typ_atoms, update_ratp
  type(typ_atoms), intent(inout):: atoms
  real(8):: xred(3,atoms%nat)
  real(8):: fcart(3,atoms%nat),energy,strten(6)
end subroutine cal_potential_forces_netsock
  subroutine init_netsock(parini)
  use mod_parini, only: typ_parini
  type(typ_parini), intent(in):: parini
end subroutine init_netsock
  subroutine send_data(pos,latvec,nat,repid,msg,nmsg,latvec_rot)
  use mod_potential, only: sock_socket, sock_inet, sock_port,sock_host,MSGLEN,sock_extra_string
  integer,intent(in):: nat, nmsg, repid
  real(8),intent(in):: latvec(3,3),pos(3,nat)
  real(8),intent(out)::latvec_rot(3,3)
  real(8):: latvec_inv(3,3),pos_back(3,nat),pos_cart(3,nat),dist_ang(6)
  real(8),parameter :: pi = 3.141592653589793239d0
  character*1024:: msg
end subroutine send_data
  subroutine get_data(etot,fcart,strten,latvec,latvec_rot,nat)
  use mod_potential, only: sock_socket, sock_inet, sock_port,sock_host,MSGLEN,sock_extra_string
  integer,intent(in) :: nat
  real(8),intent(in) :: latvec(3,3),latvec_rot(3,3)
  real(8),intent(out):: fcart(3,nat),etot,strten(6)
end subroutine get_data
       subroutine rotmat_fcart_stress_other(latvec_init,latvec_trans,rotmat)
       real(8):: latvec_init(3,3),latvec_trans(3,3),latvec_trans_inv(3,3),rotmat(3,3)
end subroutine rotmat_fcart_stress_other
       subroutine rotate_stresstensor_other(strten,rotmat)
       real(8):: strten(6),rotmat(3,3),stress(3,3)
end subroutine rotate_stresstensor_other
  subroutine final_netsock()
  character*1024:: host
end subroutine final_netsock
! ./src/potential_QSC.F90 :
subroutine init_potential_forces_qsc(atoms_t)
    use mod_atoms, only: typ_atoms
    type(typ_atoms), intent(in):: atoms_t
end subroutine init_potential_forces_qsc
subroutine cal_potential_forces_qsc(atoms)
    use mod_atoms, only: typ_atoms, update_ratp, update_rat
    type(typ_atoms), intent(inout):: atoms
end subroutine cal_potential_forces_qsc
subroutine final_potential_forces_qsc
end subroutine final_potential_forces_qsc
! ./src/potential_sec_main.F90 :
subroutine init_potential_forces_sec(parini,atoms)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
end subroutine init_potential_forces_sec
subroutine cal_potential_forces_sec(parini,atoms)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
end subroutine cal_potential_forces_sec
subroutine final_potential_forces_sec(atoms)
    use mod_atoms, only: typ_atoms
    type(typ_atoms), intent(in):: atoms
end subroutine final_potential_forces_sec
! ./src/potential_SIESTA.F90 :
subroutine init_cal_potential_forces_siesta(atoms_t)
    use mod_atoms, only: typ_atoms
    type(typ_atoms), intent(in):: atoms_t
end subroutine init_cal_potential_forces_siesta
subroutine cal_potential_forces_siesta(atoms)
    use mod_atoms, only: typ_atoms
    type(typ_atoms), intent(inout):: atoms
end subroutine cal_potential_forces_siesta
subroutine final_potential_forces_siesta
end subroutine final_potential_forces_siesta
! ./src/potential_VASP.F90 :
subroutine init_potential_forces_vasp(atoms_t)
    use mod_atoms, only: typ_atoms, set_typat
    type(typ_atoms), intent(inout):: atoms_t
end subroutine init_potential_forces_vasp
subroutine cal_potential_forces_vasp(atoms)
    use mod_atoms, only: typ_atoms, update_ratp
    type(typ_atoms), intent(inout):: atoms
end subroutine cal_potential_forces_vasp
subroutine final_potential_forces_vasp
end subroutine final_potential_forces_vasp
subroutine add_repulsive_wall(iproc,nat,rat,cellvec,fat,epot)
    integer, intent(in):: iproc, nat
    real(8), intent(in):: rat(3,nat), cellvec(3,3)
    real(8), intent(inout):: fat(3,nat), epot
end subroutine add_repulsive_wall
subroutine get_output_vasp_geopt_alborz(filename1,filename2,nat,latvec,xred,fcart,energy,strten,success)
    character(*):: filename1
    character(*):: filename2
    integer:: nat
    real(8):: fcart(3,nat),energy,strten(6),value,latvec(3,3),xred(3,nat),str_matrix(3,3),vol,a(3,3),scaling
    logical:: success
end subroutine get_output_vasp_geopt_alborz
! ./src/processors.F90 :
subroutine initprocessors
end subroutine initprocessors
subroutine finalizeprocessors
end subroutine finalizeprocessors
! ./src/processors_mod.F90 :
! ./src/propagate.F90 :
subroutine propagate(parini,nat,xred,latvec0,dxred,dlatvec,xredout,latvecout)
use mod_parini, only: typ_parini
type(typ_parini), intent(in):: parini
integer::nat,i,iat,j
real(8):: xred(3,nat),latvec(3,3),dxred(3,nat),dlatvec(3,3),xredout(3,nat),latvecout(3,3),len1,len2
real(8):: orig_angle(3),new_angle(3),axis(3),rotmat(3,3),center(3),latvec0(3,3)
end subroutine propagate
! ./src/replace.F90 :
subroutine replace(nlminx,nlmin,fp_len,nat,kid,e_wpos,ent_wpos,fp_wpos,wpos_red,&
  &wpos_latvec,spg_wpos,spgtol_wpos,fdos_wpos,&
  &e_arr,ent_arr,fp_arr,pl_arr,lat_arr,spg_arr,spgtol_arr,dos_arr,ct_arr,findsym)
  integer:: fp_len,ct_arr(nlminx),spg_arr(nlminx),nat,iat,spg_wpos
  integer:: k, nlmin,nlminx,i,kid
  real(8):: e_wpos, ent_wpos, wpos_red(3,nat),wpos_latvec(3,3),spgtol_wpos,fdos_wpos,fp_wpos(fp_len)
  real(8):: e_arr(nlminx),ent_arr(nlminx),fp_arr(fp_len,nlminx),pl_arr(3,nat,nlminx)
  real(8):: lat_arr(3,3,nlminx),spgtol_arr(nlminx),dos_arr(nlminx)
  logical:: findsym
end subroutine replace
! ./src/saddle_1s_dimer.F90 :
subroutine dimmethimproved(parini,iproc,atoms_s,nat,ndof,rat,epot,fat,curv,uvn,paropt)
    use mod_parini, only: typ_parini
    use mod_opt, only: typ_paropt
    use mod_atoms, only: typ_atoms, typ_file_info, atom_copy_old, atom_normalizevector
    type(typ_parini), intent(in):: parini
    integer, intent(in):: iproc, nat, ndof !number coordinates
    type(typ_atoms), intent(inout):: atoms_s
    real(8), intent(inout):: rat(3,nat) !positions
    real(8), intent(out):: fat(3,nat) !forces
    real(8), intent(out):: uvn(3,nat) !unit vector along the dimer, \hat{n}.
    real(8), intent(out):: epot !potential energy 
    real(8), intent(out):: curv !curvature
    type(typ_paropt), intent(inout):: paropt
end subroutine dimmethimproved
subroutine lowestcurvature(parini,iproc,atoms_s,nat,ndof,rat,uvn,fat,angletol,maxitlc,curv0,curv,nw)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, atom_ddot, atom_normalizevector
    type(typ_parini), intent(in):: parini
    integer, intent(in):: iproc, nat, ndof, maxitlc, nw
    type(typ_atoms), intent(inout):: atoms_s
    real(8), intent(in):: rat(3,nat), fat(3,nat), angletol
    real(8), intent(inout):: uvn(3,nat), curv0, curv
end subroutine lowestcurvature
subroutine rotatedimer(parini,iproc,atoms_s,nat,ndof,rat,uvn,fat,curv0,curv,fnrm)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, atom_ddot, atom_copy_old, atom_normalizevector
    type(typ_parini), intent(in):: parini
    integer, intent(in):: iproc, nat, ndof
    real(8), intent(in):: rat(3,nat), fat(3,nat)
    real(8), intent(inout):: uvn(3,nat), curv0, curv, fnrm
    type(typ_atoms), intent(inout):: atoms_s
    real(8), allocatable:: uvp(:,:) !unit vector normal to uvn ie. \hat{\phi} 
end subroutine rotatedimer
! ./src/saddle_1s.F90 :
subroutine surface_walking(parini)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, typ_atoms_arr, typ_file_info, atom_deallocate_old
    use mod_opt, only: typ_paropt
    type(typ_parini), intent(in):: parini
end subroutine surface_walking
subroutine read_input(atoms_s) !,paropt)
    use mod_atoms, only: typ_atoms
    type(typ_atoms), intent(inout):: atoms_s
end subroutine read_input
subroutine random_move_atoms(nat,atom_motion,cellvec,rat)
    integer, intent(in):: nat
    logical, intent(in):: atom_motion(3,nat)
    real(8), intent(in):: cellvec(3,3)
    real(8), intent(inout):: rat(3,nat)
end subroutine random_move_atoms
subroutine find_minima(parini,iproc,atoms_s,paropt_m,paropt_m_prec,uvn,curv,epot_m0)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, typ_file_info, atom_deallocate_old, atom_allocate_old
    use mod_opt, only: typ_paropt
    type(typ_parini), intent(in):: parini
    integer, intent(in):: iproc !, nat, ndof
    type(typ_atoms), intent(in):: atoms_s
    type(typ_paropt), intent(inout):: paropt_m, paropt_m_prec
    real(8), intent(in):: uvn(3,atoms_s%nat), curv, epot_m0
end subroutine find_minima
subroutine alongnegcurvature(iproc,atoms,uvn,c)
    use mod_atoms, only: typ_atoms, atom_deallocate_old
    integer, intent(in):: iproc
    type(typ_atoms), intent(inout):: atoms
    real(8), intent(in):: uvn(3,atoms%nat) !unit vector along the dimer, \hat{n}.
    real(8), intent(inout):: c !curvature
end subroutine alongnegcurvature
! ./src/saddle_1s_optimizer.F90 :
subroutine optimizer_saddle(parini,iproc,atoms_s,n,nr,x,f,epot,paropt,uvn)
    use mod_parini, only: typ_parini
    use mod_opt, only: typ_paropt
    use mod_atoms, only: typ_atoms
    type(typ_parini), intent(in):: parini
    integer, intent(in):: iproc, n, nr
    type(typ_atoms), intent(inout):: atoms_s
    real(8), intent(inout):: x(n), f(n), uvn(n), epot
    type(typ_paropt), intent(inout):: paropt
end subroutine optimizer_saddle
subroutine cal_potential_forces_modified(parini,iproc,atoms_s,n,x,f,epot,nr,uvn,feff,curv0,curv,fold,paropt)
    use mod_parini, only: typ_parini
    use mod_opt, only: typ_paropt
    use mod_atoms, only: typ_atoms, atom_ddot, atom_copy_old, atom_calnorm, atom_deallocate_old
    type(typ_parini), intent(in):: parini
    integer, intent(in):: iproc, n, nr
    type(typ_atoms), intent(inout):: atoms_s
    real(8), intent(inout):: x(3,atoms_s%nat)
    real(8), intent(inout):: f(3,atoms_s%nat), epot, uvn(3,atoms_s%nat), feff(3,atoms_s%nat), curv0, curv, fold(n)
    type(typ_paropt), intent(inout):: paropt
end subroutine cal_potential_forces_modified
subroutine test_convergence_saddle(n,f,curv,paropt)
    use mod_opt, only: typ_paropt
    integer, intent(in):: n
    real(8), intent(in):: f(n), curv
    type(typ_paropt), intent(inout):: paropt
end subroutine test_convergence_saddle
! ./src/saddle_1s_pot.F90 :
subroutine pot_initialize(parini,atoms,paropt,paropt_m)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_opt, only: typ_paropt
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
    type(typ_paropt), intent(inout):: paropt, paropt_m
end subroutine pot_initialize
! ./src/saddle_mod.F90 :
! ./src/save_low_conf.F90 :
subroutine save_low_conf(nat,npmin,npminx,ent_wpos,e_wpos,pos,latvec,spg,spgtol,fdos,elocmin,poslocmin,latlocmin)
  integer:: iat,nat, npmin, npminx, kmax, k 
  real(8):: e_wpos, ent_wpos, emax,spg,spgtol,fdos
  real(8):: elocmin(npminx,5)
  real(8):: pos(3,nat),latvec(3,3),poslocmin(3,nat,npminx),latlocmin(3,3,npminx)
end subroutine save_low_conf
! ./src/shortrange.F90 :
subroutine shortrange_init(atoms,shortrange,linked_lists,spline)
    use mod_atoms, only: typ_atoms
    use mod_shortrange, only: typ_shortrange
    use mod_linked_lists, only: typ_linked_lists
    use mod_spline, only: typ_spline
    type(typ_atoms), intent(in):: atoms
    type(typ_shortrange), intent(inout):: shortrange
    type(typ_linked_lists), intent(in):: linked_lists
    type(typ_spline), intent(inout):: spline
end subroutine shortrange_init
subroutine shortrange_final(linked_lists,spline)
    use mod_linked_lists, only: typ_linked_lists
    use mod_spline, only: typ_spline
    type(typ_linked_lists), intent(in):: linked_lists
    type(typ_spline), intent(inout):: spline
end subroutine shortrange_final
subroutine set_interaction(atoms,shortrange)
    use mod_atoms, only: typ_atoms
    use mod_shortrange, only: typ_shortrange
    type(typ_atoms), intent(inout):: atoms
    type(typ_shortrange), intent(inout):: shortrange
end subroutine set_interaction
subroutine cal_shortenergy(parini,shortrange,atoms,linked_lists,spline,alpha,cell,epot_short)
    use mod_parini, only: typ_parini
    use mod_shortrange, only: typ_shortrange
    use mod_linked_lists, only: typ_linked_lists
    use mod_spline, only: typ_spline
    use mod_atoms, only: typ_atoms, update_ratp
    type(typ_parini), intent(in):: parini
    type(typ_shortrange), intent(in):: shortrange
    type(typ_atoms), intent(inout):: atoms
    type(typ_linked_lists), intent(inout):: linked_lists
    type(typ_spline), intent(in):: spline
    real(8), intent(in):: alpha
    real(8), intent(out):: cell(3)
    real(8), intent(out):: epot_short !short range electrostatic energy
end subroutine cal_shortenergy
! ./src/shortrange_mod.F90 :
! ./src/slab_stress.F90 :
subroutine slab_stress(flat,fix_z)
real(8):: flat(3,3),ekin1,ekin2
logical:: fix_z
end subroutine slab_stress
! ./src/soften.F90 :
 subroutine soften_pos(parini,parres,latvec,pos_red0,ddcart,curv0,curv,res,pressure,count_soft,amass,nsoft,folder)
 use mod_parini, only: typ_parini
 type(typ_parini), intent(in):: parini
 type(typ_parini), intent(inout):: parres
 real(8) :: acell_in(3),xred_in(3,parini%nat),vel_in(3,parini%nat),etot_in,strten_in(6),vel_lat_in(3,3)
        integer:: nsoft,i,it,nit,iprec,iat
        real(8):: curv0,curv,res,pressure,count_soft,alpha
        real(8):: rxyz(3*parini%nat)
        real(8):: latvec(9),latvec_in(9)
        real(8):: ddcart(3*parini%nat)
        real(8):: rxyzcart(3*parini%nat)
        real(8):: pos_red0(3*parini%nat)
        real(8):: pos_red_in(3*parini%nat)
        real(8):: amass(parini%nat)
        real(8):: wlat(9),wlatold(9),fxyzcart(3*parini%nat)
        character(40):: filename,folder
        real(8):: pos_prev(3*parini%nat),dir_prev(3*parini%nat),dir(3*parini%nat),angle,norm
end subroutine soften_pos
 subroutine soften_lat(parini,parres,latvec,pos_red0,ddlat,curv0,curv,res,pressure,count_soft,amass,nsoft,folder)
 use mod_parini, only: typ_parini
 type(typ_parini), intent(in):: parini
 type(typ_parini), intent(inout):: parres
 real(8) :: acell_in(3),xred_in(3,parini%nat),vel_in(3,parini%nat),fcart_in(3,parini%nat),etot_in,strten_in(6),vel_lat_in(3,3)
        integer:: nsoft,i,it,nit,iprec,iat
        real(8):: curv0,curv,res,pressure,count_soft,alpha,alphalat
        real(8):: rxyz(3*parini%nat)
        real(8):: latvec(9),latvec_in(9)
        real(8):: dd(3*parini%nat)
        real(8):: ddcart(3*parini%nat)
        real(8):: rxyzcart(3*parini%nat)
        real(8):: ddlat(9)
        real(8):: ddall(3*parini%nat+9)
        real(8):: pos_red0(3*parini%nat)
        real(8):: amass(parini%nat)
        real(8):: wlat(9),wlatold(9),fxyzcart(3*parini%nat)
        character(40):: filename,folder
end subroutine soften_lat
! ./src/solve_poisson_cube.F90 :
subroutine solve_poisson(parini)
    use mod_parini, only: typ_parini
    use mod_electrostatics, only: typ_poisson
    use mod_atoms, only: typ_atoms, get_rat, update_ratp
    type(typ_parini), intent(in):: parini
end subroutine solve_poisson
! ./src/spglib_int.F90 :
subroutine get_spg(num_atom,positions,lattice,atom_types,symprec,spg)
integer:: nat, typat(num_atom), spg
  integer, intent(in) :: num_atom!, max_num_sym, is_time_reversal
  real(8), intent(in) :: symprec
  integer, intent(in), dimension(num_atom) :: atom_types
  real(8), intent(in), dimension(3, 3) :: lattice
  real(8), intent(in), dimension(3, num_atom) :: positions
end subroutine get_spg
subroutine spg_cell_refine(nat_in,nat_out,nat_max,positions,lattice,atom_types,symprec,spg)
integer:: nat, spg
  integer:: nat_in,nat_max!, max_num_sym, is_time_reversal
  integer:: nat_out
  real(8), intent(in) :: symprec
  integer, intent(inout), dimension(nat_max) :: atom_types
  real(8), intent(inout), dimension(3, 3) :: lattice
  real(8), intent(inout), dimension(3, nat_max) :: positions
end subroutine spg_cell_refine
subroutine spg_cell_primitive(nat_in,nat_out,nat_max,positions,lattice,atom_types,symprec,spg)
integer:: nat, spg
  integer:: nat_in,nat_max!, max_num_sym, is_time_reversal
  integer:: nat_out
  real(8), intent(in) :: symprec
  integer, intent(inout), dimension(nat_max) :: atom_types
  real(8), intent(inout), dimension(3, 3) :: lattice
  real(8), intent(inout), dimension(3, nat_max) :: positions
end subroutine spg_cell_primitive
! ./src/spher_harm_mathematica.F90 :
subroutine ylm_mathematica(l,m,theta,phi,ylm_r,ylm_i)
integer:: l,m
real(8):: theta,phi,ylm_r,ylm_i
end subroutine ylm_mathematica
! ./src/spline_mod.F90 :
! ./src/symfunc_mod.F90 :
! ./src/task_ann.F90 :
subroutine task_ann(parini)
    use mod_parini, only: typ_parini
    type(typ_parini), intent(in):: parini
end subroutine task_ann
! ./src/task_bader.F90 :
subroutine task_bader(parini)
    use mod_parini, only: typ_parini
    type(typ_parini), intent(in):: parini
end subroutine task_bader
! ./src/task_confcomp.F90 :
subroutine conf_comp(parini)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms_all, atom_all_allocate, atom_all_deallocate, set_rcov
    type(typ_parini), intent(in):: parini
end subroutine conf_comp
subroutine set_fpall_ann(atoms_all)
    use mod_atoms, only: typ_atoms_all, set_rat
    type(typ_atoms_all), intent(inout):: atoms_all
end subroutine set_fpall_ann
subroutine set_fpall_angle(atoms_all)
    use mod_atoms, only: typ_atoms_all, set_rat
    type(typ_atoms_all), intent(inout):: atoms_all
end subroutine set_fpall_angle
subroutine set_fpall_distance(atoms_all)
    use mod_atoms, only: typ_atoms_all, set_rat
    type(typ_atoms_all), intent(inout):: atoms_all
end subroutine set_fpall_distance
subroutine build_images(atoms,natpmax,natp,ratp)
    use mod_atoms, only: typ_atoms, get_rat
    type(typ_atoms), intent(in):: atoms
    integer, intent(in):: natpmax
    integer, intent(inout):: natp
    real(8), intent(inout):: ratp(3,natpmax)
end subroutine build_images
! ./src/task_genconf.F90 :
subroutine task_genconf(parini)
    use mod_genconf, only: typ_genconf
    use mod_parini, only: typ_parini
    type(typ_parini), intent(inout):: parini
end subroutine task_genconf
! ./src/task_geopt.F90 :
subroutine geopt(parini)
    use mod_parini, only: typ_parini
    use mod_opt, only: typ_paropt
    use mod_atoms, only: typ_atoms_arr, typ_file_info, set_ndof, atom_deallocate
    type(typ_parini), intent(inout):: parini !poscar_getsystem must be called from parser
end subroutine geopt
subroutine init_geopt(parini,paropt,paropt_prec)
    use mod_parini, only: typ_parini
    use mod_opt, only: typ_paropt
    type(typ_parini), intent(in):: parini
    type(typ_paropt), intent(inout):: paropt, paropt_prec
end subroutine init_geopt
! ./src/task_lammps.F90 :
subroutine lammps_task(parini)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms_arr, atom_copy, atom_deallocate, set_typat
    use, intrinsic :: ISO_C_binding, only : C_double, C_ptr, C_int, C_FUNPTR
    type(typ_parini), intent(in):: parini
end subroutine lammps_task
subroutine lammps_write(parini,atoms)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, update_ratp
    type(typ_parini), intent(in):: parini
    type(typ_atoms):: atoms
end subroutine lammps_write
! ./src/task_linkedlist.F90 :
subroutine  linkedlist_test(parini)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, type_pairs, typ_file_info, atom_deallocate_old
    use mod_linked_lists, only: typ_linked_lists, typ_pia_arr
    type(typ_parini), intent(in):: parini
end subroutine linkedlist_test
subroutine callinkedlist(parini,atoms,rcut,posat1st,nim,conf)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, type_pairs, update_ratp
    use mod_linked_lists, only: typ_linked_lists
    type(typ_parini):: parini
    type(typ_atoms):: atoms 
    type(type_pairs):: posat1st(atoms%nat) 
    integer:: istat, nim(atoms%nat)
    integer:: iat, jat, niat,kat,kkz,conf
    real(8):: sclinv,cell(3) ,rcut
end subroutine callinkedlist
subroutine sort_alborz(i ,j ,k,conf)
    integer :: i, j, k
    integer ::conf
end subroutine sort_alborz
subroutine sort2_alborz(i ,j ,k,conf,num)
    integer :: i, j, k
    integer ::conf,num
end subroutine sort2_alborz
subroutine genrandomconf(atoms,numb,conf)
    use mod_atoms, only: typ_atoms
    use mod_atoms, only: typ_atoms, typ_file_info, update_rat
    integer ::mat,conf
    character(2):: numb
    type(typ_atoms):: atoms 
end subroutine genrandomconf
! ./src/task_minhocao.F90 :
subroutine task_minhocao(parini,parres)
 use mod_parini, only: typ_parini
 type(typ_parini), intent(inout):: parini
 type(typ_parini), intent(inout):: parres
  real(8), parameter :: beta1=1.10d0,beta2=1.10d0,beta3=1.d0/1.10d0
  real(8), parameter :: alpha1=1.d0/1.10d0,alpha2=1.10d0
end subroutine task_minhocao
! ./src/task_miscellaneous.F90 :
subroutine miscellaneous_task(parini)
    use mod_parini, only: typ_parini
    type(typ_parini), intent(in):: parini
end subroutine miscellaneous_task
! ./src/task_mod.F90 :
! ./src/task_potential.F90 :
subroutine alborz_as_potential_init(nat,sat)
    integer, intent(in):: nat
    character(5), intent(in):: sat(nat)
end subroutine alborz_as_potential_init
subroutine alborz_as_potential_get(boundcond,nat,cellvec,rat,sat,fat,epot,stress)
    character(*), intent(in):: boundcond
    integer, intent(in):: nat
    real(8), intent(in):: rat(3,nat), cellvec(3,3)
    character(5), intent(in):: sat(nat)
    real(8), intent(out):: fat(3,nat), epot, stress(3,3)
end subroutine alborz_as_potential_get
subroutine alborz_as_potential_final
end subroutine alborz_as_potential_final
! ./src/task_single_point.F90 :
subroutine single_point_task(parini)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms_arr, typ_file_info, set_ndof, atom_deallocate
    type(typ_parini), intent(inout):: parini !poscar_getsystem must be called from parser
end subroutine single_point_task
subroutine read_poscar_for_single_point(parini,atoms)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, atom_allocate_old, update_rat
    type(typ_parini), intent(inout):: parini !poscar_getsystem must be called from parser
    type(typ_atoms):: atoms
end subroutine read_poscar_for_single_point
! ./src/task_testforces.F90 :
subroutine task_testforces(parini)
    use mod_parini, only: typ_parini
    type(typ_parini), intent(in):: parini
end subroutine task_testforces
subroutine testforces_fd(parini)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, typ_atoms_arr, atom_copy_old, atom_deallocate
    type(typ_parini), intent(in):: parini
end subroutine testforces_fd
subroutine teststress_fd(parini)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, update_ratp
    type(typ_parini), intent(in):: parini
    integer, parameter:: m=3
    real(8), parameter:: h=5.d-5
end subroutine teststress_fd
subroutine teststress_fd_cellvec(parini)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, update_ratp, update_rat
    type(typ_parini), intent(in):: parini
    integer, parameter:: m=3
    real(8), parameter:: h=5.d-5
end subroutine teststress_fd_cellvec
! ./src/test_free_bps.F90 :
subroutine test_free_bps(parini)
    use mod_parini, only: typ_parini
    use mod_electrostatics, only: typ_poisson
    use mod_atoms, only: typ_atoms, update_ratp, get_rat
    type(typ_parini), intent(in):: parini
end subroutine test_free_bps
! ./src/tightbinding.F90 :
subroutine set_indorb(partb,atoms)
    use mod_atoms, only: typ_atoms
    use mod_tightbinding, only: typ_partb
    type(typ_partb), intent(inout):: partb
    type(typ_atoms), intent(inout):: atoms
end subroutine set_indorb
subroutine gammaenergy(pia_arr,linked_lists,partb,atoms,natsi,pplocal)
    use mod_linked_lists, only: typ_pia_arr, typ_linked_lists
    use mod_atoms, only: typ_atoms, set_typat
    use mod_potl, only: potl_typ
    use mod_tightbinding, only: typ_partb, lenosky
    type(typ_linked_lists), intent(in):: linked_lists
    type(typ_pia_arr), intent(in):: pia_arr
    type(typ_partb), intent(inout):: partb
    type(typ_atoms), intent(inout):: atoms
    integer, intent(in):: natsi
    type(potl_typ), intent(inout):: pplocal
end subroutine gammaenergy
subroutine gammamat(pia_arr,linked_lists,partb,atoms,natsi,flag2,pplocal)
    use mod_linked_lists, only: typ_pia_arr, typ_linked_lists
    use mod_tightbinding, only: typ_partb
    use mod_atoms, only: typ_atoms
    use mod_potl, only: potl_typ
    type(typ_linked_lists), intent(in):: linked_lists
    type(typ_pia_arr), intent(in):: pia_arr
    type(typ_partb), intent(inout):: partb
    type(typ_atoms), intent(inout):: atoms
    integer, intent(in):: natsi, flag2
    type(potl_typ), intent(in):: pplocal
end subroutine gammamat
subroutine forcediagonalizeg(partb)
    use mod_tightbinding, only: typ_partb
    type(typ_partb), intent(inout):: partb
end subroutine forcediagonalizeg
subroutine gammacoupling(pia,ib,partb,atoms,flag2,atomtypei,atomtypej,pplocal,rem)
    use mod_linked_lists, only: typ_pia !, typ_linked_lists
    use mod_tightbinding, only: typ_partb, lenosky
    use mod_atoms, only: typ_atoms
    use mod_potl, only: potl_typ
    type(typ_pia), intent(in):: pia
    type(typ_partb), intent(inout):: partb
    type(typ_atoms), intent(inout):: atoms
    integer, intent(in):: ib, flag2, atomtypei, atomtypej
    type(potl_typ), intent(in):: pplocal
    real(8), intent(out):: rem(partb%nstride,partb%nstride)
end subroutine gammacoupling
subroutine slatercoupling(u,r,hgen,dhgen,flag2,mat)
    real(8), intent(in):: r, u(3)
    real(8), intent(in):: hgen(4), dhgen(4)
    integer, intent(in):: flag2
    real(8), intent(out):: mat(4,4)
end subroutine slatercoupling
subroutine yfdocclocal(partb)
    use mod_tightbinding, only: typ_partb, lenosky
    type(typ_partb), intent(inout):: partb
    integer, parameter:: fdmaxit= 20000 !Maximum number of NR iterations
    real(8), parameter:: fdepsocc=1.d-9 !Allowed error in number of electrons
end subroutine yfdocclocal
subroutine Hamiltonian_der(u,flag2,mat)
    real(8), intent(in):: u(3)
    integer, intent(in):: flag2
    real(8), intent(out):: mat(4,4)
end subroutine hamiltonian_der
! ./src/tightbinding_mod.F90 :
! ./src/timing_mod.F90 :
! ./src/torque_cell.F90 :
subroutine torque_cell(latvec0,vlat,torquenrm)
real(8), intent(in)    :: latvec0(3,3)
real(8), intent(inout) :: vlat(3,3),torquenrm
end subroutine torque_cell
! ./src/tosifumi.F90 :
subroutine set_tosifumi(atoms,tosifumi)
    use mod_atoms, only: typ_atoms
    use mod_shortrange, only: typ_tosifumi
    type(typ_atoms), intent(inout):: atoms
    type(typ_tosifumi), intent(inout):: tosifumi
end subroutine set_tosifumi
subroutine coulomb_free_direct(atoms)
    use mod_atoms, only: typ_atoms, update_ratp
    type(typ_atoms), intent(inout):: atoms
end subroutine coulomb_free_direct
subroutine calenergyforces(atoms,tosifumi)
    use mod_atoms, only: typ_atoms, update_ratp
    use mod_shortrange, only: typ_tosifumi
    type(typ_atoms), intent(inout):: atoms
    type(typ_tosifumi), intent(in):: tosifumi
end subroutine calenergyforces
subroutine tosifumi_parameters(s,p)
    character(6), intent(out):: s(10)
    real(8), intent(out):: p(5,10)
end subroutine tosifumi_parameters
! ./src/unitsconversion_mod.F90 :
! ./src/write_restart.F90 :
subroutine winter(parini,nat,units,ent_pos,e_pos,pos_red,pos_latvec,pos_fcart,pos_strten,nlminx,nlmin,npminx,& 
   &ent_arr,e_arr,ct_arr,spg_arr,spgtol_arr,dos_arr,pl_arr,lat_arr,f_arr,str_arr,fp_arr,fp_len,ent_delta,fp_delta,& 
   &eref,ediff,ekinetic,ekinetic_max,dt,nsoften,char_type,ntypat,typat,fixat,fixlat,pressure)
  use mod_parini, only: typ_parini
  type(typ_parini), intent(in):: parini
  integer, intent(in) :: nlminx,nlmin,nsoften,nat,npminx,fp_len
  real(8), intent(in) :: eref,ediff,ekinetic,dt,e_pos,ent_pos,ekinetic_max,ent_delta,fp_delta
  real(8), intent(in) :: pos_latvec(3,3) 
  real(8), intent(in) :: pos_strten(6) 
  real(8), dimension(nlminx),      intent(in) :: ent_arr,e_arr,spgtol_arr,dos_arr
  real(8), dimension(3,3,nlminx),  intent(in) :: lat_arr
  real(8), dimension(6,nlminx),    intent(in) :: str_arr
  real(8), dimension(3,nat,nlminx),intent(in) :: pl_arr,f_arr
  real(8), dimension(fp_len,nlminx),intent(in):: fp_arr
  integer, dimension(nlminx),      intent(in) :: ct_arr,spg_arr
  character(2), intent(in):: char_type(ntypat) 
  integer, intent(in):: ntypat 
  integer, intent(in):: typat(nat) 
  real(8), intent(in):: pressure 
  real(8), intent(in):: pos_red(3,nat) 
  real(8), intent(in):: pos_fcart(3,nat) 
  logical :: fixat(nat),fixlat(7)
  character(len=40) :: units 
end subroutine winter
subroutine wtioput(ediff,ekinetic,ekinetic_max,nsoften)
  integer:: nsoften
  real(8):: ediff, ekinetic,ekinetic_max
end subroutine wtioput
! ./t2.f90 :
! ./t3.f90 :
end interface
end module mod_interface
!*****************************************************************************************
