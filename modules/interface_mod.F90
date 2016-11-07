!***************************************************************************************************
module mod_interface
    implicit none
interface
! ./modules/alborz_as_potential_mod.F90 :
! ./modules/ann_mod.F90 :
! ./modules/atoms_mod.F90 :
! ./modules/constants_mod.F90 :
! ./modules/dynamics_mod.F90 :
! ./modules/electrostatics_mod.F90 :
! ./modules/fsockets.F90 :
! ./modules/genconf_mod.F90 :
! ./modules/linked_lists_mod.F90 :
! ./modules/minhopp_mod.F90 :
! ./modules/opt_mod.F90 :
! ./modules/parini_mod.F90 :
! ./modules/potential_mod.F90 :
! ./modules/processors_mod.F90 :
! ./modules/saddle_mod.F90 :
! ./modules/shortrange_mod.F90 :
! ./modules/spline_mod.F90 :
! ./modules/task_mod.F90 :
! ./modules/tightbinding_mod.F90 :
! ./modules/timing_mod.F90 :
! ./modules/unitsconversion_mod.F90 :
! ./src/alborz.F90 :
! ./src/alborz_init_final.F90 :
subroutine alborz_init(parini,file_ini)
    use mod_task, only: typ_file_ini, time_start
    use mod_parini, only: typ_parini
    implicit none
    type(typ_file_ini), intent(inout):: file_ini
    type(typ_parini), intent(inout):: parini
end subroutine alborz_init
subroutine alborz_initialize_timing_categories
    implicit none
end subroutine alborz_initialize_timing_categories
subroutine alborz_final(parini,file_ini)
    use mod_parini, only: typ_parini
    use mod_task, only: typ_file_ini, time_start, time_end
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_file_ini), intent(inout):: file_ini
end subroutine alborz_final
subroutine init_random_seed(parini)
    use mod_parini, only: typ_parini
    implicit none
    type(typ_parini), intent(in):: parini
end subroutine init_random_seed
subroutine set_atomc_types_info(parini)
    use mod_parini, only: typ_parini
    implicit none
    type(typ_parini), intent(inout):: parini
end subroutine set_atomc_types_info
! ./src/ann_basic.F90 :
subroutine ann_allocate(ekf,ann_arr)
    use mod_ann, only: typ_ann_arr, typ_ekf
    implicit none
    type(typ_ekf), intent(in):: ekf
    type(typ_ann_arr), intent(inout):: ann_arr
end subroutine ann_allocate
subroutine ann_deallocate(ann_arr)
    use mod_ann, only: typ_ann_arr
    implicit none
    type(typ_ann_arr), intent(inout):: ann_arr
end subroutine ann_deallocate
! ./src/ann_best_symfunc.F90 :
subroutine ann_best_symfunc(parini)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann, typ_symfunc_arr
    use mod_atoms, only: typ_atoms_arr
    implicit none
    type(typ_parini), intent(in):: parini
end subroutine ann_best_symfunc
subroutine cal_symfunc_diversity(n_tot,his,ann,disparity)
    use mod_ann, only: typ_ann
    implicit none
    integer, intent(in):: n_tot
    real(8), intent(in):: his(1000,n_tot)
    type(typ_ann), intent(inout):: ann
    real(8), intent(inout):: disparity
end subroutine cal_symfunc_diversity
subroutine gbounds_distro(ann,atoms_arr,strmess)
    use mod_atoms, only: typ_atoms_arr
    use mod_ann, only: typ_ann
    implicit none
    type(typ_ann), intent(inout):: ann
    type(typ_atoms_arr), intent(inout):: atoms_arr
    character(*), intent(in):: strmess
end subroutine gbounds_distro
! ./src/ann_check_symmetry_function.F90 :
subroutine ann_check_symmetry_function(parini)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr, typ_symfunc_arr, typ_symfunc
    use mod_atoms, only: typ_atoms_arr
    implicit none
    type(typ_parini), intent(in):: parini
end subroutine ann_check_symmetry_function
! ./src/ann_ekf_behler.F90 :
subroutine ekf_behler(parini,ann_arr,symfunc_train,symfunc_valid,atoms_train,atoms_valid,ekf)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr, typ_symfunc_arr, typ_ekf
    use mod_atoms, only: typ_atoms, typ_atoms_arr
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_symfunc_arr), intent(inout):: symfunc_train, symfunc_valid
    type(typ_atoms_arr), intent(inout):: atoms_train
    type(typ_atoms_arr), intent(inout):: atoms_valid
    type(typ_ekf), intent(inout):: ekf
end subroutine ekf_behler
! ./src/ann_ekf_rivals.F90 :
subroutine ekf_rivals(parini,ann_arr,symfunc_train,symfunc_valid,atoms_train,atoms_valid,ekf)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr, typ_symfunc_arr, typ_ekf
    use mod_atoms, only: typ_atoms, typ_atoms_arr
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_symfunc_arr), intent(inout):: symfunc_train, symfunc_valid
    type(typ_atoms_arr), intent(inout):: atoms_train
    type(typ_atoms_arr), intent(inout):: atoms_valid
    type(typ_ekf), intent(inout):: ekf
end subroutine ekf_rivals
subroutine ekf_rivals_tmp(parini,ann_arr,symfunc_train,symfunc_valid,atoms_train,atoms_valid,ekf)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr, typ_symfunc_arr, typ_ekf
    use mod_atoms, only: typ_atoms, typ_atoms_arr
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_symfunc_arr), intent(inout):: symfunc_train, symfunc_valid
    type(typ_atoms_arr), intent(inout):: atoms_train
    type(typ_atoms_arr), intent(inout):: atoms_valid
    type(typ_ekf), intent(inout):: ekf
end subroutine ekf_rivals_tmp
subroutine set_ref_energy(parini,atoms_train,atoms_ref,ind)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, typ_atoms_arr
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms_arr), intent(in):: atoms_train
    type(typ_atoms_arr), intent(inout):: atoms_ref
    integer, intent(out):: ind(200)
end subroutine set_ref_energy
subroutine analyze_epoch_init(parini,atoms_train,ann_arr)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms_arr
    use mod_ann, only: typ_ann_arr
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms_arr), intent(in):: atoms_train
    type(typ_ann_arr), intent(inout):: ann_arr
end subroutine analyze_epoch_init
subroutine analyze_epoch_print(parini,iter,atoms_train,ann_arr)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms_arr
    use mod_ann, only: typ_ann_arr
    implicit none
    type(typ_parini), intent(in):: parini
    integer, intent(in):: iter
    type(typ_atoms_arr), intent(in):: atoms_train
    type(typ_ann_arr), intent(inout):: ann_arr
end subroutine analyze_epoch_print
! ./src/ann_evaluate.F90 :
subroutine ann_evaluate_subtask(parini)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann, typ_symfunc_arr
    use mod_atoms, only: typ_atoms_arr
    implicit none
    type(typ_parini), intent(in):: parini
end subroutine ann_evaluate_subtask
! ./src/ann_gen_symmetry_function.F90 :
subroutine ann_gen_symmetry_function(parini)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr, typ_symfunc_arr, typ_symfunc
    use mod_atoms, only: typ_atoms_arr
    implicit none
    type(typ_parini), intent(in):: parini
end subroutine ann_gen_symmetry_function
! ./src/ann_io.F90 :
subroutine read_input_ann(parini,iproc,ann_arr)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    implicit none
    type(typ_parini), intent(in):: parini
    integer, intent(in):: iproc
    type(typ_ann_arr), intent(inout):: ann_arr
end subroutine read_input_ann
subroutine read_symmetry_functions(parini,iproc,ifile,ann,rcut)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann
    implicit none
    type(typ_parini), intent(in):: parini
    integer, intent(in):: iproc, ifile
    type(typ_ann), intent(inout):: ann
    real(8), intent(out):: rcut
end subroutine read_symmetry_functions
subroutine set_radial_atomtype(parini,sat1,ityp)
    use mod_parini, only: typ_parini
    implicit none
    type(typ_parini), intent(in):: parini
    character(*), intent(in):: sat1
    integer, intent(out):: ityp(1)
end subroutine set_radial_atomtype
subroutine set_angular_atomtype(parini,sat1,sat2,ityp)
    use mod_parini, only: typ_parini
    implicit none
    type(typ_parini), intent(in):: parini
    character(*), intent(in):: sat1, sat2
    integer, intent(out):: ityp(2)
end subroutine set_angular_atomtype
subroutine write_ann(parini,filename,ann)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann
    implicit none
    type(typ_parini), intent(in):: parini
    character(*):: filename
    type(typ_ann), intent(in):: ann
end subroutine write_ann
subroutine read_ann(parini,ann_arr)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
end subroutine read_ann
subroutine read_data(parini,filename_list,atoms_arr)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms_all, typ_atoms_arr
    implicit none
    type(typ_parini), intent(in):: parini
    character(*), intent(in):: filename_list
    type(typ_atoms_arr), intent(inout):: atoms_arr
end subroutine read_data
! ./src/ann_lm.F90 :
subroutine ann_lm(parini,ann_arr,atoms_train,atoms_valid,symfunc_train,symfunc_valid,ekf)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr, typ_symfunc_arr, typ_ekf
    use mod_atoms, only: typ_atoms_arr
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr):: ann_arr
    type(typ_ekf):: ekf
    type(typ_atoms_arr):: atoms_train
    type(typ_atoms_arr):: atoms_valid
    type(typ_symfunc_arr):: symfunc_train
    type(typ_symfunc_arr):: symfunc_valid
end subroutine ann_lm
subroutine fcn_least_squares(m,n,x,fvec,fjac,ldfjac,iflag,parini,ann_arr,atoms_train,atoms_valid,symfunc_train,symfunc_valid,ekf)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr, typ_symfunc_arr, typ_ekf
    use mod_atoms, only: typ_atoms, typ_atoms_arr
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_atoms_arr), intent(inout):: atoms_train, atoms_valid
    type(typ_symfunc_arr), intent(inout):: symfunc_train, symfunc_valid
    type(typ_ekf), intent(inout):: ekf
    integer:: m, n, ldfjac, iflag
    real(8):: x(n), fvec(m), fjac(ldfjac,n)
end subroutine fcn_least_squares
! ./src/ann_pot_atom.F90 :
subroutine cal_ann_atombased(parini,atoms,symfunc,ann_arr,ekf)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_ann, only: typ_ann_arr, typ_symfunc, typ_ekf
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_symfunc), intent(inout):: symfunc
    type(typ_ekf), intent(inout):: ekf
end subroutine cal_ann_atombased
! ./src/ann_pot_eem1.F90 :
subroutine cal_ann_eem1(parini,atoms,symfunc,ann_arr,ekf)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_ann, only: typ_ann_arr, typ_symfunc, typ_ekf
    use mod_electrostatics, only: typ_ewald_p3d
    use mod_linked_lists, only: typ_pia_arr
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_symfunc), intent(inout):: symfunc
    type(typ_ekf), intent(inout):: ekf
end subroutine cal_ann_eem1
subroutine get_qat_from_chi(parini,ann_arr,atoms,ewald_p3d,a)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_atoms, only: typ_atoms
    use mod_electrostatics, only: typ_ewald_p3d
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_atoms), intent(inout):: atoms
    type(typ_ewald_p3d), intent(inout):: ewald_p3d
    real(8), intent(inout):: a(atoms%nat+1,atoms%nat+1)
end subroutine get_qat_from_chi
subroutine get_qat_from_chi_dir(parini,ann_arr,atoms,a)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_atoms), intent(inout):: atoms
    real(8), intent(inout):: a(atoms%nat+1,atoms%nat+1)
end subroutine get_qat_from_chi_dir
subroutine cal_electrostatic_eem1(parini,str_job,atoms,ann_arr,epot_c,a,ewald_p3d)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_ann, only: typ_ann_arr
    use mod_electrostatics, only: typ_ewald_p3d
    implicit none
    type(typ_parini), intent(in):: parini
    character(*), intent(in):: str_job
    type(typ_atoms), intent(inout):: atoms
    type(typ_ann_arr), intent(inout):: ann_arr
    real(8), intent(out):: epot_c
    real(8), intent(inout):: a(atoms%nat+1,atoms%nat+1)
    type(typ_ewald_p3d), intent(inout):: ewald_p3d
end subroutine cal_electrostatic_eem1
subroutine cal_electrostatic_ann(parini,atoms,ann_arr,a,ewald_p3d)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_ann, only: typ_ann_arr
    use mod_electrostatics, only: typ_ewald_p3d
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
    type(typ_ann_arr), intent(inout):: ann_arr
    real(8), intent(in):: a(atoms%nat+1,atoms%nat+1)
    type(typ_ewald_p3d), intent(inout):: ewald_p3d
end subroutine cal_electrostatic_ann
subroutine charge_analysis(parini,atoms,ann_arr)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_ann, only: typ_ann_arr
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(in):: atoms
    type(typ_ann_arr), intent(inout):: ann_arr
end subroutine charge_analysis
subroutine get_qat_from_chi_iter(parini,ann_arr,atoms,a)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_atoms), intent(inout):: atoms
    real(8), intent(inout):: a(atoms%nat+1,atoms%nat+1)
end subroutine get_qat_from_chi_iter
subroutine cal_ugradient(parini,ewald_p3d,ann_arr,atoms,g,qtot)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_atoms, only: typ_atoms
    use mod_electrostatics, only: typ_ewald_p3d
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ewald_p3d),intent(inout):: ewald_p3d
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_atoms), intent(inout):: atoms
    real(8), intent(out):: g(atoms%nat), qtot
end subroutine cal_ugradient
subroutine get_qat_from_chi_operator(parini,ewald_p3d,ann_arr,atoms)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_atoms, only: typ_atoms
    use mod_electrostatics, only: typ_ewald_p3d
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_atoms), intent(inout):: atoms
    type(typ_ewald_p3d),intent(inout):: ewald_p3d
end subroutine get_qat_from_chi_operator
! ./src/ann_pot_eem2.F90 :
subroutine cal_ann_eem2(parini,atoms,symfunc,ann_arr,ekf)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_ann, only: typ_ann_arr, typ_symfunc, typ_ekf
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_symfunc), intent(inout):: symfunc
    type(typ_ekf), intent(inout):: ekf
end subroutine cal_ann_eem2
subroutine get_qat_from_chi2(parini,ann_arr,atoms,a)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_atoms), intent(inout):: atoms
    real(8), intent(inout):: a(atoms%nat+1,atoms%nat+1)
end subroutine get_qat_from_chi2
subroutine cal_electrostatic_eem2(parini,str_job,atoms,ann_arr,epot_es,a)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_ann, only: typ_ann_arr
    implicit none
    type(typ_parini), intent(in):: parini
    character(*), intent(in):: str_job
    type(typ_atoms), intent(in):: atoms
    type(typ_ann_arr), intent(inout):: ann_arr
    real(8), intent(out):: epot_es
    real(8), intent(out):: a(atoms%nat+1,atoms%nat+1)
end subroutine cal_electrostatic_eem2
! ./src/ann_pot_main.F90 :
subroutine cal_ann_main(parini,atoms,symfunc,ann_arr,ekf)
    use mod_tightbinding, only: typ_partb
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_ann, only: typ_ann_arr, typ_symfunc, typ_ekf
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_symfunc), intent(inout):: symfunc
    type(typ_ekf), intent(inout):: ekf
end subroutine cal_ann_main
! ./src/ann_pot_tb.F90 :
subroutine cal_ann_tb(parini,partb,atoms,ann_arr,symfunc,ekf)
    use mod_parini, only: typ_parini
    use mod_tightbinding, only: typ_partb
    use mod_potl, only: potl_typ
    use mod_atoms, only: typ_atoms
    use mod_ann, only: typ_ann_arr, typ_symfunc, typ_ekf
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
    type(typ_ekf), intent(inout):: ekf
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_symfunc), intent(inout):: symfunc
    type(typ_partb), intent(inout):: partb
    real(8):: hgen_der(4,1:atoms%nat,1:atoms%nat)  , ttxyz !derivative of 
end subroutine cal_ann_tb
subroutine lenoskytb_ann(partb,atoms,natsi,count_md)
    use mod_atoms, only: typ_atoms
    use mod_potl, only: potl_typ
    use mod_tightbinding, only: typ_partb, lenosky
    implicit none
    type(typ_partb), intent(inout):: partb
    type(typ_atoms), intent(inout):: atoms
    integer, intent(in):: natsi
    real(8), intent(inout):: count_md
end subroutine lenoskytb_ann
! ./src/ann_process.F90 :
subroutine cal_architecture(ann,epot)
    use mod_ann, only: typ_ann
    implicit none
    type(typ_ann), intent(inout):: ann
    real(8), intent(inout):: epot
end subroutine cal_architecture
subroutine cal_architecture_1hiddenlayer(ann,epot)
    use mod_ann, only: typ_ann
    implicit none
    type(typ_ann), intent(inout):: ann
    real(8), intent(inout):: epot
end subroutine cal_architecture_1hiddenlayer
subroutine cal_architecture_2hiddenlayer(ann,epot)
    use mod_ann, only: typ_ann
    implicit none
    type(typ_ann), intent(inout):: ann
    real(8), intent(inout):: epot
end subroutine cal_architecture_2hiddenlayer
subroutine cal_architecture_der(ann,epot)
    use mod_ann, only: typ_ann
    implicit none
    type(typ_ann), intent(inout):: ann
    real(8), intent(inout):: epot
end subroutine cal_architecture_der
subroutine cal_architecture_der_1hiddenlayer(ann,epot)
    use mod_ann, only: typ_ann
    implicit none
    type(typ_ann), intent(inout):: ann
    real(8), intent(inout):: epot
end subroutine cal_architecture_der_1hiddenlayer
subroutine cal_architecture_der_2hiddenlayer(ann,epot)
    use mod_ann, only: typ_ann
    implicit none
    type(typ_ann), intent(inout):: ann
    real(8), intent(inout):: epot
end subroutine cal_architecture_der_2hiddenlayer
! ./src/ann_symfunc_atom_behler.F90 :
subroutine symmetry_functions_driver(parini,ann_arr,atoms,symfunc)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr, typ_symfunc
    use mod_atoms, only: typ_atoms
    use mod_linked_lists, only: typ_pia_arr !,typ_linked_lists
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_atoms), intent(in):: atoms
    type(typ_symfunc), intent(inout):: symfunc
end subroutine symmetry_functions_driver
subroutine symmetry_functions_g02_atom(ann_arr,pia,ib,iat,isat,jsat,symfunc)
    use mod_ann, only: typ_ann_arr, typ_symfunc
    use mod_linked_lists, only: typ_pia
    implicit none
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_pia), intent(in):: pia
    integer, intent(in):: ib, iat, isat, jsat
    type(typ_symfunc), intent(inout):: symfunc
end subroutine symmetry_functions_g02_atom
subroutine symmetry_functions_g04_atom(ann_arr,isat,iat,jsat,jat_maincell,ksat,kat_maincell,rij,rik,rjk,drij,drik,drjk,fcij,fcdij,fcik,fcdik,fcjk,fcdjk,symfunc)
    use mod_ann, only: typ_ann_arr, typ_symfunc
    implicit none
    type(typ_ann_arr), intent(inout):: ann_arr
    integer, intent(in):: isat, iat, jsat, jat_maincell, ksat, kat_maincell
    real(8), intent(in):: rij, rik, rjk, drij(3), drik(3), drjk(3), fcij, fcdij, fcik, fcdik, fcjk, fcdjk
    type(typ_symfunc), intent(inout):: symfunc
end subroutine symmetry_functions_g04_atom
subroutine symmetry_functions_g05_atom(ann_arr,piaij,piaik,ibij,ibik,iat,isat,jsat,ksat,symfunc)
    use mod_ann, only: typ_ann_arr, typ_symfunc
    use mod_linked_lists, only: typ_pia
    implicit none
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_pia), intent(in):: piaij, piaik
    integer, intent(in):: ibij, ibik, isat, iat, jsat, ksat
    type(typ_symfunc), intent(inout):: symfunc
end subroutine symmetry_functions_g05_atom
subroutine symmetry_functions_g06_atom(ann,iat,jat_maincell,r,dr,fc,fcd)
    use mod_ann, only: typ_ann
    implicit none
    type(typ_ann), intent(inout):: ann
    integer, intent(in):: iat, jat_maincell
    real(8), intent(in):: r, dr(3), fc, fcd
end subroutine symmetry_functions_g06_atom
subroutine symmetry_functions_g02(ann,iat,atoms,i0)
    use mod_ann, only: typ_ann
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_ann), intent(inout):: ann
    integer, intent(in):: iat
    type(typ_atoms), intent(in):: atoms
    integer, intent(inout):: i0
end subroutine symmetry_functions_g02
subroutine symmetry_functions_g04(ann,iat,atoms,i0)
    use mod_ann, only: typ_ann
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_ann), intent(inout):: ann
    integer, intent(in):: iat
    type(typ_atoms), intent(in):: atoms
    integer, intent(inout):: i0
end subroutine symmetry_functions_g04
subroutine symmetry_functions_g05(ann,iat,atoms,i0)
    use mod_ann, only: typ_ann
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_ann), intent(inout):: ann
    integer, intent(in):: iat
    type(typ_atoms), intent(in):: atoms
    integer, intent(inout):: i0
end subroutine symmetry_functions_g05
subroutine symmetry_functions_g06(ann,iat,atoms,i0)
    use mod_ann, only: typ_ann
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_ann), intent(inout):: ann
    integer, intent(in):: iat
    type(typ_atoms), intent(in):: atoms
    integer, intent(inout):: i0
end subroutine symmetry_functions_g06
function cutoff_function(r, rc) result(fc)
    implicit none
    real(8), intent(in):: r, rc
    real(8):: fc, pi
end function cutoff_function
function cutoff_function_der(r, rc) result(fcd)
    implicit none
    real(8), intent(in):: r, rc
    real(8):: fcd, pi
end function cutoff_function_der
! ./src/ann_symfunc_atom_stefan.F90 :
subroutine symmetry_functions_driver_stefan(parini,ann_arr,atoms,symfunc)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr, typ_symfunc
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_atoms), intent(inout):: atoms
    type(typ_symfunc), intent(inout):: symfunc
end subroutine symmetry_functions_driver_stefan
subroutine fingerprint_power(nat, npl, alat, rxyz, rcov, fpall)
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
    use mod_ann, only: typ_ann_arr, typ_symfunc
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_atoms), intent(inout):: atoms
    type(typ_symfunc), intent(inout):: symfunc
    logical, intent(in):: apply_gbounds
end subroutine symmetry_functions
! ./src/ann_symfunc_pair_behler.F90 :
subroutine symmetry_functions_driver_bond(parini,ann_arr,atoms,symfunc)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr, typ_symfunc
    use mod_atoms, only: typ_atoms
    use mod_linked_lists, only: typ_pia_arr !,typ_linked_lists
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_atoms), intent(in):: atoms
    type(typ_symfunc), intent(inout):: symfunc
end subroutine symmetry_functions_driver_bond
subroutine symmetry_functions_driver_bond_tmp(ann_arr,atoms)
    use mod_ann, only: typ_ann_arr
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_atoms), intent(in):: atoms
end subroutine symmetry_functions_driver_bond_tmp
subroutine symmetry_functions_g01_bond(ann_arr,ib,pia,symfunc)
    use mod_linked_lists, only: typ_pia
    use mod_ann, only: typ_ann_arr, typ_symfunc
    implicit none
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_pia), intent(in):: pia
    integer, intent(in):: ib
    type(typ_symfunc), intent(inout):: symfunc
end subroutine symmetry_functions_g01_bond
subroutine symmetry_functions_g02_bond(ann_arr,iat,jat,rij,drij,fcij,fcdij,rik,drik,fcik,fcdik,rjk,drjk,fcjk,fcdjk)
    use mod_ann, only: typ_ann_arr
    implicit none
    type(typ_ann_arr), intent(inout):: ann_arr
    integer, intent(in):: iat, jat
    real(8), intent(in):: fcij, fcik, fcjk, fcdij, fcdik, fcdjk
    real(8), intent(in):: rij, rik, rjk
    real(8), intent(in):: drij(3), drik(3), drjk(3)
end subroutine symmetry_functions_g02_bond
subroutine symmetry_functions_g04_bond(ann_arr,iat,jat,rij,drij,fcij,fcdij,rik,drik,fcik,fcdik,rjk,drjk,fcjk,fcdjk)
    use mod_ann, only: typ_ann_arr
    implicit none
    type(typ_ann_arr), intent(inout):: ann_arr
    integer, intent(in):: iat, jat
    real(8), intent(in):: fcij, fcik, fcjk, fcdij, fcdik, fcdjk
    real(8), intent(in):: rij, rik, rjk
    real(8), intent(in):: drij(3), drik(3), drjk(3)
end subroutine symmetry_functions_g04_bond
! ./src/ann_train.F90 :
subroutine ann_train(parini)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr, typ_symfunc_arr, typ_ekf
    use mod_atoms, only: typ_atoms_arr
    implicit none
    type(typ_parini), intent(in):: parini
end subroutine ann_train
subroutine apply_gbounds_atom(parini,ann_arr,atoms_arr,symfunc_arr)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr, typ_symfunc_arr
    use mod_atoms, only: typ_atoms_arr
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(in):: ann_arr
    type(typ_atoms_arr), intent(inout):: atoms_arr
    type(typ_symfunc_arr), intent(inout):: symfunc_arr
end subroutine apply_gbounds_atom
subroutine apply_gbounds_bond(ann_arr,atoms_arr,symfunc_arr)
    use mod_ann, only: typ_ann_arr, typ_symfunc_arr
    use mod_atoms, only: typ_atoms_arr
    implicit none
    type(typ_ann_arr), intent(in):: ann_arr
    type(typ_atoms_arr), intent(inout):: atoms_arr
    type(typ_symfunc_arr), intent(inout):: symfunc_arr
end subroutine apply_gbounds_bond
subroutine prepare_atoms_arr(parini,ann_arr,atoms_arr)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_atoms, only: typ_atoms_arr
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(in):: ann_arr
    type(typ_atoms_arr), intent(inout):: atoms_arr
end subroutine prepare_atoms_arr
subroutine set_ebounds(ann_arr,atoms_train,atoms_valid,symfunc_train,symfunc_valid)
    use mod_ann, only: typ_ann_arr, typ_symfunc_arr
    use mod_atoms, only: typ_atoms_arr
    implicit none
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_atoms_arr), intent(in):: atoms_train, atoms_valid
    type(typ_symfunc_arr), intent(inout):: symfunc_train, symfunc_valid
end subroutine set_ebounds
subroutine ann_evaluate(parini,iter,ann_arr,symfunc_arr,atoms_arr,ifile,partb)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr, typ_symfunc_arr, typ_ekf
    use mod_atoms, only: typ_atoms, typ_atoms_arr
    use mod_tightbinding, only: typ_partb
    implicit none
    type(typ_parini), intent(in):: parini
    integer, intent(in):: iter
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_symfunc_arr), intent(inout):: symfunc_arr
    type(typ_atoms_arr), intent(inout):: atoms_arr
    integer, intent(in):: ifile
    type(typ_partb), optional, intent(inout):: partb
end subroutine ann_evaluate
subroutine eval_cal_ann_main(parini,atoms,symfunc,ann_arr)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_tightbinding, only: typ_partb
    use mod_ann, only: typ_ann_arr, typ_symfunc, typ_ekf
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_symfunc), intent(inout):: symfunc
end subroutine eval_cal_ann_main
subroutine set_gbounds(parini,ann_arr,atoms_arr,strmess,symfunc_arr)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr, typ_symfunc_arr
    use mod_atoms, only: typ_atoms_arr
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_atoms_arr), intent(inout):: atoms_arr
    character(*), intent(in):: strmess
    type(typ_symfunc_arr), intent(inout):: symfunc_arr
end subroutine set_gbounds
subroutine write_symfunc(parini,iconf,atoms_arr,strmess,symfunc_arr)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr, typ_symfunc_arr
    use mod_atoms, only: typ_atoms_arr
    implicit none
    type(typ_parini), intent(in):: parini
    integer, intent(in):: iconf
    type(typ_atoms_arr), intent(inout):: atoms_arr
    character(*), intent(in):: strmess
    type(typ_symfunc_arr), intent(inout):: symfunc_arr
end subroutine write_symfunc
subroutine read_symfunc(parini,iconf,ann_arr,atoms_arr,strmess,symfunc_arr)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr, typ_symfunc_arr
    use mod_atoms, only: typ_atoms_arr
    use mod_linked_lists, only: typ_pia_arr
    implicit none
    type(typ_parini), intent(in):: parini
    integer, intent(in):: iconf
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_atoms_arr), intent(inout):: atoms_arr
    character(*), intent(in):: strmess
    type(typ_symfunc_arr), intent(inout):: symfunc_arr
end subroutine read_symfunc
subroutine save_gbounds(parini,ann_arr,atoms_arr,strmess,symfunc_arr)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr, typ_symfunc_arr
    use mod_atoms, only: typ_atoms_arr
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_atoms_arr), intent(inout):: atoms_arr
    character(*), intent(in):: strmess
    type(typ_symfunc_arr), intent(inout):: symfunc_arr
end subroutine save_gbounds
subroutine convert_x_ann(n,x,ann)
    use mod_ann, only: typ_ann
    implicit none
    integer, intent(in):: n
    real(8), intent(in):: x(n)
    type(typ_ann), intent(inout):: ann
end subroutine convert_x_ann
subroutine convert_ann_epotd(ann,n,epotd)
    use mod_ann, only: typ_ann
    implicit none
    type(typ_ann), intent(in):: ann
    integer, intent(in):: n
    real(8), intent(inout):: epotd(n)
end subroutine convert_ann_epotd
subroutine randomize_data_order(atoms_arr)
    use mod_atoms, only: typ_atoms_arr, typ_atoms
    implicit none
    type(typ_atoms_arr), intent(inout):: atoms_arr
end subroutine randomize_data_order
! ./src/ann_weights_init.F90 :
subroutine set_annweights(parini,ekf)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ekf
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ekf), intent(inout):: ekf
end subroutine set_annweights
! ./src/basic_atoms.F90 :
subroutine atom_allocate(atoms,nat,natim,nfp)
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_atoms), intent(inout):: atoms
    integer, intent(in):: nat, natim, nfp
end subroutine atom_allocate
subroutine atom_deallocate(atoms)
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_atoms), intent(inout):: atoms
end subroutine atom_deallocate
subroutine atom_allocate_old(atoms,nat,natim,nfp,sat,vat,amass,fat,bemoved,qat,zat,rcov,typat)
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_atoms), intent(inout):: atoms
    integer, intent(in):: nat, natim, nfp
    logical, optional, intent(in):: sat, vat, amass
    logical, optional, intent(in):: fat, bemoved, qat, zat, rcov, typat
end subroutine atom_allocate_old
subroutine atom_deallocate_old(atoms,sat,rat,ratim,vat,amass,fat,bemoved,qat,zat,rcov,fp,typat)
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_atoms), intent(inout):: atoms
    logical, optional, intent(in):: sat, rat, ratim, vat, amass
    logical, optional, intent(in):: fat, bemoved, qat, zat, rcov, fp, typat
end subroutine atom_deallocate_old
subroutine atom_all_allocate(atoms_all,ratall,fatall,epotall,fpall,qtotall)
    use mod_atoms, only: typ_atoms_all
    implicit none
    type(typ_atoms_all), intent(inout):: atoms_all
    logical, optional, intent(in):: ratall, fatall, epotall, fpall, qtotall
end subroutine atom_all_allocate
subroutine atom_all_deallocate(atoms_all,ratall,fatall,epotall,fpall,qtotall)
    use mod_atoms, only: typ_atoms_all
    implicit none
    type(typ_atoms_all), intent(inout):: atoms_all
    logical, optional, intent(in):: ratall, fatall, epotall, fpall, qtotall
end subroutine atom_all_deallocate
subroutine atom_copy(at_inp,at_out,str_message)
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_atoms), intent(in):: at_inp
    type(typ_atoms), intent(inout):: at_out
    character(*):: str_message
end subroutine atom_copy
subroutine atom_copy_old(at_inp,at_out,str_message,sat,rat,ratim,vat,amass,fat,bemoved,qat,zat,rcov,fp,typat)
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_atoms), intent(in):: at_inp
    type(typ_atoms), intent(inout):: at_out
    character(*):: str_message
    logical, optional, intent(in):: sat, rat, ratim, vat, amass
    logical, optional, intent(in):: fat, bemoved, qat, zat, rcov, fp, typat
end subroutine atom_copy_old
subroutine atom_build_periodic_images(atoms,rcut)
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_atoms), intent(inout):: atoms
    real(8), intent(in):: rcut
end subroutine atom_build_periodic_images
subroutine set_typat(atoms)
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_atoms), intent(inout):: atoms
end subroutine set_typat
subroutine set_ndof(atoms)
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_atoms), intent(inout):: atoms
end subroutine set_ndof
subroutine bemoved2string(bemoved,str_motion)
    implicit none
    logical:: bemoved(3)
    character(3):: str_motion
end subroutine bemoved2string
subroutine string2bemoved(str_motion,bemoved)
    implicit none
    character(3):: str_motion
    logical:: bemoved(3)
end subroutine string2bemoved
subroutine atoms_all_writexyz(filename,fn_position,atoms_all,strkey)
    use mod_atoms, only: typ_atoms_all
    implicit none
    character(*), intent(in):: filename, fn_position, strkey
    type(typ_atoms_all), intent(in):: atoms_all
end subroutine atoms_all_writexyz
subroutine atom_normalizevector(nat,bemoved,vxyz)
    implicit none
    integer, intent(in):: nat
    logical, intent(in):: bemoved(3,nat)
    real(8), intent(inout):: vxyz(3,nat)
end subroutine atom_normalizevector
function atom_ddot(nat,vec1,vec2,bemoved) result(res)
    implicit none
    integer, intent(in):: nat
    real(8), intent(in):: vec1(3,nat), vec2(3,nat)
    logical, intent(in):: bemoved(3,nat)
    real(8):: res
end function atom_ddot
subroutine atom_calnorm(nat,bemoved,vec,vnrm)
    implicit none
    integer, intent(in):: nat
    logical, intent(in):: bemoved(3,nat)
    real(8), intent(in):: vec(3,nat)
    real(8), intent(out):: vnrm
end subroutine atom_calnorm
subroutine atom_calmaxforcecomponent(nat,bemoved,vec,vmax)
    implicit none
    integer, intent(in):: nat
    logical, intent(in):: bemoved(3,nat)
    real(8), intent(in):: vec(3,nat)
    real(8), intent(out):: vmax
end subroutine atom_calmaxforcecomponent
subroutine checkallatomstobeincell(nat,rat,cellvec,allatomsincell)
    implicit none
    integer:: iat,nat
    real(8):: rat(3,nat),cellvec(3,3)
    logical:: allatomsincell
end subroutine checkallatomstobeincell
subroutine determinexyzminmax(nat,rat,cellvec,xmin,ymin,zmin,xmax,ymax,zmax)
    implicit none
    integer:: iat,nat
    real(8):: rat(3,nat),cellvec(3,3)
    real(8):: xmin,ymin,zmin,xmax,ymax,zmax,x,y,z
end subroutine determinexyzminmax
subroutine set_rcov(atoms)
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_atoms), intent(inout) :: atoms
end subroutine set_rcov
subroutine set_qat(atoms)
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_atoms), intent(inout) :: atoms
end subroutine set_qat
subroutine sat_to_iatom(sat,iatom)
    implicit none
    character(*), intent(in) :: sat
    integer, intent(out) :: iatom
end subroutine sat_to_iatom
subroutine iatom_to_sat(iatom,sat)
    implicit none
    integer, intent(in) :: iatom
    character(*), intent(out) :: sat
end subroutine iatom_to_sat
! ./src/basic.F90 :
subroutine elim_white_space(string)
    implicit none
    character(256), intent(inout):: string
end subroutine elim_white_space
function delta_kronecker(i,j) result(delta)
    implicit none
    integer, intent(in):: i, j
    real(8):: delta
end function delta_kronecker
subroutine backtocell_alborz(nat,latvec,rxyz_red)
    implicit none
    integer, intent(in):: nat
    real(8), intent(in):: latvec(3,3)
    real(8), intent(inout):: rxyz_red(3,nat)
end subroutine backtocell_alborz
subroutine getvol_alborz(cellvec,vol)
    implicit none
    real(8), intent(in):: cellvec(3,3)
    real(8), intent(out):: vol
end subroutine getvol_alborz
subroutine pbc_distance1_alborz(cellvec,xred_1,xred_2,distance2,dxyz)
    implicit none
    real(8), intent(in):: cellvec(3,3), xred_1(3), xred_2(3)
    real(8), intent(inout):: distance2, dxyz(3)
end subroutine pbc_distance1_alborz
subroutine n_rep_dim_alborz(cellvec,rcut,nec1,nec2,nec3)
    implicit none
    real(8), intent(in):: cellvec(3,3), rcut
    integer, intent(out):: nec1, nec2, nec3
end subroutine n_rep_dim_alborz
subroutine nveclatvec_alborz(cellvec,vn)
    implicit none
    real(8), intent(in) :: cellvec(3,3)
    real(8), intent(out):: vn(3,3)
end subroutine nveclatvec_alborz
subroutine dist2plane_alborz(r1,vn,r0,dist)
    implicit none
    real(8), intent(in):: r1(3), vn(3), r0(3)
    real(8), intent(out):: dist
end subroutine dist2plane_alborz
subroutine write_atomic_file_ascii_alborz(filename,nat,xred,latvec0,energy,pressure,printval1,printval2,kinds)
implicit none
integer:: nat,natin,iat
character(40):: filename,units
real(8):: pos(3,nat),xred(3,nat),latvec(3,3),latvec0(3,3),dproj(6),rotmat(3,3),v(3,3),ucvol
real(8):: energy, etotal, enthalpy, enthalpy_at,pressure,printval1,printval2
integer:: Kinds(nat)
end subroutine write_atomic_file_ascii_alborz
subroutine dproj2latvec_alborz(dproj,cellvec)
    implicit none
    real(8), intent(in):: dproj(6)
    real(8), intent(out):: cellvec(3,3)
end subroutine dproj2latvec_alborz
subroutine latvec2dproj_alborz(dproj,latvec,rotmat,rxyz,nat)
    implicit none
    integer, intent(in):: nat
    real(8),intent(inout):: dproj(6), latvec(3,3), rotmat(3,3), rxyz(3,nat)
end subroutine latvec2dproj_alborz
subroutine cross_product_alborz(a,b,c)
    implicit none
    real(8), intent(in):: a(3), b(3)
    real(8), intent(out):: c(3)
end subroutine cross_product_alborz
subroutine rotation_alborz(angle,axe,rotmat)
    implicit none
    real(8), intent(in):: angle
    real(8), intent(in):: axe(3)
    real(8), intent(out):: rotmat(3,3)
end subroutine rotation_alborz
subroutine fxyz_cart2int_alborz(nat,v_cart,cv,v_int)
    implicit none
    integer, intent(in):: nat
    real(8), intent(in):: v_cart(3,nat), cv(3,3)
    real(8), intent(out):: v_int(3,nat)
end subroutine fxyz_cart2int_alborz
subroutine fxyz_red2cart(nat,fint,cv,fcart)
    implicit none
    integer, intent(in):: nat
    real(8), intent(in):: fint(3,nat), cv(3,3)
    real(8), intent(out):: fcart(3,nat)
end subroutine fxyz_red2cart
subroutine count_words(str,n)
    implicit none
    character(*), intent(in):: str
    integer, intent(out):: n
end subroutine count_words
subroutine count_substring(str1,str2,n)
    implicit none
    character(*), intent(in) :: str1, str2
    integer, intent(out):: n
end subroutine count_substring
! ./src/basic_utilities.F90 :
subroutine elim_moment_alborz(nat,atomic_vector)
  implicit none
  integer, intent(in):: nat
  real(8), intent(inout):: atomic_vector(3,nat)
end subroutine elim_moment_alborz
subroutine elim_moment_mass(nat,atomic_vector,atomic_mass)
  implicit none
  integer, intent(in):: nat
  real(8), intent(inout):: atomic_vector(3,nat)
  real(8), intent(in):: atomic_mass(nat)
end subroutine elim_moment_mass
subroutine elim_torque_reza_alborz(nat,rat0,fat)
  implicit none
  integer, intent(in) :: nat
  real(8), dimension(3*nat), intent(in) :: rat0
  real(8), dimension(3*nat), intent(inout) :: fat
  real(8), dimension(3*nat) :: rat
  real(8), dimension(3*nat,3) :: vrot
end subroutine elim_torque_reza_alborz
subroutine mycross(a,b,c)
  implicit none
  real(8), dimension(3), intent(in):: a,b
  real(8), dimension(3), intent(out):: c
end subroutine mycross
subroutine moment_of_inertia_alborz(nat,rat,teneria,evaleria)
  implicit none
  integer, intent(in) :: nat
  real(8), dimension(3,nat), intent(in) :: rat
  real(8), dimension(3), intent(out) :: evaleria
  real(8), dimension(3,3), intent(out) :: teneria
end subroutine moment_of_inertia_alborz
subroutine normalizevector_alborz(n,v)
    implicit none
    integer, intent(in):: n
    real(8), intent(inout):: v(n)
end subroutine normalizevector_alborz
subroutine calnorm(n,v,vnrm)
    implicit none
    integer, intent(in):: n
    real(8), intent(in):: v(n)
    real(8), intent(out):: vnrm
end subroutine calnorm
subroutine calmaxforcecomponent(n,v,vmax)
    implicit none
    integer, intent(in):: n
    real(8), intent(in):: v(n)
    real(8), intent(out):: vmax
end subroutine calmaxforcecomponent
subroutine rxyz_cart2int_alborz(nat,latvec,rxyzcart,rxyzint)
    implicit none
    integer:: nat,iat
    real(8):: rxyzint(3,nat), rxyzcart(3,nat), latvec(3,3), latvecinv(3,3)
end subroutine rxyz_cart2int_alborz
subroutine rxyz_int2cart_alborz(nat,cellvec,rat_int,rat_cart)
    implicit none
    integer, intent(in):: nat
    real(8), intent(in):: cellvec(3,3), rat_int(3,nat)
    real(8), intent(inout):: rat_cart(3,nat)
end subroutine rxyz_int2cart_alborz
subroutine invertmat_alborz(a,ainv)
    implicit none
    real(8),intent(in):: a(3,3)
    real(8),intent(out):: ainv(3,3)
end subroutine invertmat_alborz
subroutine convertupper(str)
    character(*), intent(inout):: str
end subroutine convertupper
subroutine convertlower(str)
    character(*), intent(inout) :: str
end subroutine convertlower
subroutine check_whether_time_exceeded
    implicit none
end subroutine check_whether_time_exceeded
subroutine expdist(n,x)
    implicit none
    integer, intent(in):: n
    real(8), intent(out):: x(n)
end subroutine expdist
subroutine gausdist_alborz(n,x)
    implicit none
    integer, intent(in) ::n
    real(8), intent(out) :: x(n)
end subroutine gausdist_alborz
subroutine randdist(a,n,x)
    implicit none
    real(8), intent(in) :: a
    integer, intent(in):: n
    real(8), intent(out) :: x(n)
end subroutine randdist
subroutine hunt2(n,x,p,ip)
    implicit none
    integer, intent(in):: n
    real(8), intent(in):: x(n), p
    integer, intent(out):: ip
end subroutine hunt2
subroutine hpsort(n,ra)
    real*8 ::ra(n)
end subroutine hpsort
! ./src/cell_linkedlists.F90 :
subroutine linkedlists_init(parini,atoms,cell,linked_lists)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_electrostatics, only: typ_linked_lists
    implicit none
    type(typ_parini), intent(in):: parini 
    type(typ_atoms), intent(in):: atoms
    real(8), intent(out):: cell(3)
    type(typ_linked_lists), intent(inout):: linked_lists
end subroutine linkedlists_init
subroutine linkedlists_final(linked_lists)
    use mod_electrostatics, only: typ_linked_lists
    implicit none
    type(typ_linked_lists), intent(inout):: linked_lists
end subroutine linkedlists_final
subroutine prepprimelast(atoms,linked_lists,cell)
    use mod_atoms, only: typ_atoms
    use mod_electrostatics, only: typ_linked_lists
    implicit none
    type(typ_atoms), intent(in):: atoms
    type(typ_linked_lists), intent(inout):: linked_lists
    real(8):: cell(3)
end subroutine prepprimelast
subroutine make_list_new(parini,atoms,linked_lists,cell)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_electrostatics, only: typ_linked_lists
    implicit none
    type(typ_parini), intent(in):: parini 
    type(typ_atoms), intent(in):: atoms
    type(typ_linked_lists), intent(inout):: linked_lists
    real(8), intent(in):: cell(3)
end subroutine make_list_new
subroutine determine_sclimitsphere(linked_lists)
    use mod_electrostatics, only: typ_linked_lists
    implicit none
    type(typ_linked_lists), intent(inout):: linked_lists
end subroutine determine_sclimitsphere
subroutine call_linkedlist(parini,atoms,linked_lists,pia_arr)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, type_pairs
    use mod_linked_lists, only: typ_linked_lists, typ_pia_arr
    implicit none
    type(typ_parini), intent(in):: parini 
    type(typ_atoms), intent(in):: atoms 
    type(typ_linked_lists), intent(inout):: linked_lists
    type(typ_pia_arr), intent(inout):: pia_arr
end subroutine call_linkedlist
! ./src/es_coulomb_p3d_bias.F90 :
subroutine erfc_surface_zero(parini,atoms,ewald_p3d,nlayer)
    use mod_electrostatics, only: typ_ewald_p3d
    use mod_atoms, only: typ_atoms
    use mod_electrostatics, only: typ_linked_lists
    use mod_parini, only: typ_parini
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ewald_p3d), intent(inout):: ewald_p3d
    type(typ_atoms), intent(inout):: atoms
    integer:: nimat !number of image atoms.
    integer:: nlayer, igpx,igpy,igpz,mx,my,mz, mlimnlayer
end subroutine erfc_surface_zero
subroutine sollaplaceq(poisson_p3d,hz,cell,vl,vu)
    use mod_electrostatics, only: typ_poisson_p3d
    implicit none
    type(typ_poisson_p3d), intent(inout):: poisson_p3d
    real(8):: vl, vu , zlmzu , sinhzlmzu, zlmzuinv
    real(8):: cell(3)
    real(8):: hz , vlmvu, vlzumvuzl 
end subroutine sollaplaceq
 subroutine calculate_force_ener_plane(atoms,ewald_p3d,epot)
    use mod_electrostatics, only: typ_ewald_p3d
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_ewald_p3d), intent(inout):: ewald_p3d
    type(typ_atoms), intent(inout):: atoms
    real(8):: x,y,z ,t,tl ,epot ,t1,t2
    real(8):: fatp(3,atoms%nat) 
end subroutine calculate_force_ener_plane
subroutine LG_weight(nlx,nly,nlz,hx,hy,hz,wx,wy,wz)
    implicit none
    integer:: nlx ,nly, nlz !number of point for Lagrange interpolation
    real(8):: hx ,hy ,hz , hxp , hyp, hzp
    real(8):: wx(nlx), wy(nly), wz(nlz) 
end subroutine lg_weight
subroutine LGW(n, w, h, x, LGx, DLGx, ix1, nbgp)
implicit none
    integer:: ix, jx,kx, ix1,ixo ,n ,n2 ,nbgp
    real(8):: w(n), q(n),LGx(n),DLGx(n),h ,x ,diffx,x1,protot,t
end subroutine lgw
subroutine LGW4(n, w, h, x, LGx, DLGx, ix1, nbgp)
implicit none
    integer:: ix, jx, ix1,ixo ,n ,n2 ,nbgp
    real(8):: w(n), q(n), qinv(n), LGx(n), DLGx(n), h ,x ,diffx,x1,protot
end subroutine lgw4
subroutine surface_charge(ewald_p3d,pot_short,vl,vu)
    use mod_electrostatics, only: typ_ewald_p3d
    implicit none
    type(typ_ewald_p3d), intent(inout):: ewald_p3d
    real(8):: t, tt ,density(ewald_p3d%poisson_p3d%ngpx,ewald_p3d%poisson_p3d%ngpy,2),vl,vu
    real(8)::hgzinv,pi,pot_layerl,pot_layeru,pot_short(ewald_p3d%poisson_p3d%ngpx,ewald_p3d%poisson_p3d%ngpy,2,4)
end subroutine surface_charge
subroutine determine_limitsphere(ewald_p3d,mboundg,mboundgy,nbgpx,nbgpy,nbgpz)
    use mod_electrostatics, only: typ_ewald_p3d
    implicit none
    type(typ_ewald_p3d), intent(inout):: ewald_p3d
    integer:: ix, iy, iz, mboundg(2,-nbgpy:nbgpy,-nbgpz:nbgpz), mboundgy(2,-nbgpz:nbgpz)
    integer:: nbgpx, nbgpy, nbgpz
end subroutine determine_limitsphere
! ./src/es_coulomb_p3d.F90 :
subroutine construct_ewald_p3d(parini,atoms,ewald_p3d)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_electrostatics, only: typ_ewald_p3d
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(in):: atoms
    type(typ_ewald_p3d), intent(inout):: ewald_p3d
end subroutine construct_ewald_p3d
subroutine destruct_ewald_p3d(parini,atoms,ewald_p3d)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_electrostatics, only: typ_ewald_p3d
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(in):: atoms
    type(typ_ewald_p3d), intent(inout):: ewald_p3d
end subroutine destruct_ewald_p3d
subroutine calculate_forces_energy(parini,ewald_p3d,atoms)
    use mod_electrostatics, only: typ_ewald_p3d
    use mod_atoms, only: typ_atoms
    use mod_parini, only: typ_parini
    implicit none
    type(typ_ewald_p3d), intent(inout):: ewald_p3d
    type(typ_atoms), intent(inout):: atoms
    type(typ_parini), intent(in):: parini
end subroutine calculate_forces_energy
subroutine calparam(parini,atoms,ewald_p3d_rough,ewald_p3d)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_electrostatics, only: typ_ewald_p3d
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(in):: atoms
    type(typ_ewald_p3d), intent(in):: ewald_p3d_rough
    type(typ_ewald_p3d), intent(inout):: ewald_p3d
end subroutine calparam
subroutine determine_glimitsphere(ewald_p3d)
    use mod_electrostatics, only: typ_ewald_p3d
    implicit none
    type(typ_ewald_p3d), intent(inout):: ewald_p3d
end subroutine determine_glimitsphere
subroutine putgaussgrid(parini,atoms,ewald_p3d,gausswidth)
    use mod_atoms, only: typ_atoms
    use mod_electrostatics, only: typ_ewald_p3d
    use mod_parini, only: typ_parini
    implicit none
    type(typ_atoms), intent(in):: atoms
    type(typ_ewald_p3d), intent(inout):: ewald_p3d
    real(8), intent(in):: gausswidth(atoms%nat)
    type(typ_parini), intent(in):: parini
end subroutine putgaussgrid
subroutine longerange_forces(atoms,ewald_p3d,gausswidth)
    use mod_atoms, only: typ_atoms
    use mod_electrostatics, only: typ_ewald_p3d
    implicit none
    type(typ_atoms), intent(inout):: atoms
    type(typ_ewald_p3d), intent(inout):: ewald_p3d
    real(8), intent(in):: gausswidth(atoms%nat)
end subroutine longerange_forces
! ./src/es_coulomb_spline.F90 :
subroutine build_shortrange_spline(shortrange,spline,rcut,a)
    use mod_shortrange, only: typ_shortrange
    use mod_spline, only: typ_spline
    implicit none
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
    implicit none
    external:: cal_f_fd_fdd
    real(16), intent(in):: rcutq !first and second cutoff for the function
    real(16), intent(in):: hspq
    real(16), intent(in):: aq !prefacor in exponent of exponential function
    integer, intent(in):: nsp
    real(16), intent(out):: fdspq(0:3,0:nsp)
    real(16), intent(out):: fspq(0:4,0:nsp-1)
end subroutine build_spline
subroutine erf_over_r(r,a,hsp,func,funcder,funcsecder)
    implicit none 
    real(16), intent(in):: r
    real(16), intent(in):: a
    real(16), intent(in):: hsp
    real(16), intent(out):: func
    real(16), intent(out):: funcder
    real(16), intent(out):: funcsecder
end subroutine erf_over_r
subroutine one_over_r6(r,a,hsp,func,funcder,funcsecder)
    implicit none 
    real(16), intent(in):: r
    real(16), intent(in):: a
    real(16), intent(in):: hsp
    real(16), intent(out):: func
    real(16), intent(out):: funcder
    real(16), intent(out):: funcsecder
end subroutine one_over_r6
subroutine one_over_r8(r,a,hsp,func,funcder,funcsecder)
    implicit none 
    real(16), intent(in):: r
    real(16), intent(in):: a
    real(16), intent(in):: hsp
    real(16), intent(out):: func
    real(16), intent(out):: funcder
    real(16), intent(out):: funcsecder
end subroutine one_over_r8
subroutine exp_ar(r,a,hsp,func,funcder,funcsecder)
    implicit none 
    real(16), intent(in):: r
    real(16), intent(in):: a
    real(16), intent(in):: hsp
    real(16), intent(out):: func
    real(16), intent(out):: funcder
    real(16), intent(out):: funcsecder
end subroutine exp_ar
! ./src/es_hartee_bps.F90 :
subroutine cal_hartree_pot_bps(ewald_p3d,atoms,ehartree)
    use mod_atoms, only: typ_atoms
    use mod_electrostatics, only: typ_ewald_p3d
    implicit none
    type(typ_ewald_p3d),intent(inout):: ewald_p3d
    type(typ_atoms), intent(inout):: atoms
    real(8), intent(out):: ehartree
end subroutine cal_hartree_pot_bps
subroutine construct_ewald_bps(parini,atoms,ewald_p3d)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_electrostatics, only: typ_ewald_p3d
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(in):: atoms
    type(typ_ewald_p3d), intent(inout):: ewald_p3d
end subroutine construct_ewald_bps
subroutine destruct_ewald_bps(ewald_p3d)
    use mod_electrostatics, only: typ_ewald_p3d
    implicit none
    type(typ_ewald_p3d), intent(inout):: ewald_p3d
end subroutine destruct_ewald_bps
subroutine set_ngp_bps(ewald_p3d_rough,ewald_p3d)
    use mod_electrostatics, only: typ_ewald_p3d
    implicit none
    type(typ_ewald_p3d), intent(in):: ewald_p3d_rough
    type(typ_ewald_p3d), intent(inout):: ewald_p3d
end subroutine set_ngp_bps
! ./src/es_hartee_fourier.F90 :
subroutine kwald(iverbose,nat,rat,ratred,qat,cv,gwsq,ecut,ehartree,fat,eqd,stress,celldv)
    implicit none
    integer, intent(in):: iverbose, nat
    real(8), intent(in):: rat(3,nat), ratred(3,nat), qat(nat)
    real(8), intent(in):: cv(3,3), gwsq(nat), ecut
    real(8), intent(out):: fat(3,nat), eqd(nat), ehartree, stress(3,3), celldv(3,3)
end subroutine kwald
subroutine kwald_samare(iverbose,nat,rat,ratred,qat,cv,alphasq,ecut,ehartree,fat,eqd,stress,celldv)
    implicit none
    integer, intent(in):: iverbose, nat
    real(8), intent(in):: rat(3,nat), ratred(3,nat), qat(nat)
    real(8), intent(in):: cv(3,3), alphasq, ecut
    real(8), intent(out):: fat(3,nat), eqd(nat), ehartree, stress(3,3), celldv(3,3)
end subroutine kwald_samare
! ./src/es_hartee_main.F90 :
subroutine get_hartree(parini,ewald_p3d,atoms,gausswidth,ehartree,g)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_electrostatics, only: typ_ewald_p3d
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ewald_p3d),intent(inout):: ewald_p3d
    type(typ_atoms), intent(inout):: atoms
    real(8), intent(in):: gausswidth(atoms%nat)
    real(8), intent(out):: ehartree, g(atoms%nat)
end subroutine get_hartree
subroutine get_g_from_pot(parini,atoms,ewald_p3d,gausswidth,g)
    use mod_atoms, only: typ_atoms
    use mod_electrostatics, only: typ_ewald_p3d
    use mod_parini, only: typ_parini
    implicit none
    type(typ_atoms), intent(in):: atoms
    type(typ_ewald_p3d), intent(inout):: ewald_p3d
    type(typ_parini), intent(in):: parini
    real(8):: g(atoms%nat) 
    real(8), intent(in):: gausswidth(atoms%nat)
end subroutine get_g_from_pot
subroutine real_part(parini,atoms,gausswidth,alpha,epotreal,gg,stress)
    use mod_parini, only: typ_parini
    use mod_linked_lists, only: typ_linked_lists
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
    real(8), intent(in):: gausswidth(atoms%nat),alpha
    real(8):: gg(atoms%nat),rr
    real(8):: cell(3) , vol, stress(3,3)
    real(8)::epotreal,alphatwoinv,ralphasq,rbetasq,rbetainv,alphasq,betainv
end subroutine real_part
! ./src/es_hartee_p3d.F90 :
subroutine ps2dp1df_construction(poisson_p3d)
    use mod_electrostatics, only: typ_poisson_p3d
    implicit none
    type(typ_poisson_p3d), intent(inout):: poisson_p3d
end subroutine ps2dp1df_construction
subroutine ps2dp1df_destruction(poisson_p3d)
    use mod_electrostatics, only: typ_poisson_p3d
    implicit none
    type(typ_poisson_p3d), intent(inout):: poisson_p3d
end subroutine ps2dp1df_destruction
subroutine calculate_potener_pot(poisson_p3d,cell,hx,hy,hz,epot,beta)
    use mod_electrostatics, only: typ_poisson_p3d
    implicit none
    type(typ_poisson_p3d), intent(inout):: poisson_p3d
    real(8):: cell(3) !cell array contains size of the simulation box.
    real(8):: hx, hy, hz
    real(8):: epot
    real(8), optional:: beta !beta is proportion to dipole moment as it is in paper.
end subroutine calculate_potener_pot
subroutine solsyslinequ(poisson_p3d,hz,cell,beta_arg)
    use mod_electrostatics, only: typ_poisson_p3d
    implicit none
    type(typ_poisson_p3d), intent(inout):: poisson_p3d
    real(8):: hz, cell(3)
    real(8), optional:: beta_arg !beta_arg is proportion to dipole moment as it is in paper.
    real(8):: d(poisson_p3d%ngpz+2*8) !nem was replaced by 8 to be able to compile interface_mod.F90
    real(8):: e1(poisson_p3d%ngpz), e2(poisson_p3d%ngpz-1), c(poisson_p3d%ngpz)
end subroutine solsyslinequ
subroutine fdcoeff(ngpz,e1,e2,g,gsq,hz,hzsq)
    implicit none
    integer::ngpz,ngpzm1
    real(8)::e1(ngpz) !Diagonal elements of the matrix
    real(8)::e2(ngpz-1) !Offdiagonal elements of the matrix
    real(8)::gsqhzsq,a,hz,hzsq,gsq,g,diagonal_fl,diagonal,offdiagonal
end subroutine fdcoeff
subroutine prepare00(ngpz,nem,f,c,hz)
    implicit none
    integer::ngpz,i,nem,n
    real(8)::f(ngpz+2*nem),c(ngpz),eta(6),hz 
end subroutine prepare00
subroutine prepare(ngpz,nem,f,c,gsq,hz,hzsq)
    implicit none
    integer::ngpz,i,nem,n
    real(8)::f(ngpz+2*nem),c(ngpz),eta(6),gsq,hz,hzsq,a
end subroutine prepare
subroutine prepcoeff(hz,eta,coefftot1,coefftoti,coefftotn)
    implicit none
    real(8)::hz,coeff(16,8),coefftot1(17),coefftoti(17),coefftotn(17),eta(6)
end subroutine prepcoeff
subroutine calbeta(hzsq,ngpz,analc00,beta)
    implicit none
    integer::ngpz,iz
    real(8)::analc00(ngpz),hzsq,beta
end subroutine calbeta
! ./src/forcefield.F90 :
subroutine forcefield_init(parini,atoms)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_electrostatics, only: typ_ewald_p3d
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
end subroutine forcefield_init
subroutine calculate_forces_energy_ff(parini,atoms)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
end subroutine calculate_forces_energy_ff
subroutine forcefield_final(parini,atoms)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(in):: atoms
end subroutine forcefield_final
! ./src/genconf_diatomic.F90 :
subroutine genconf_diatomic(parini,genconf)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms_all, typ_file_info
    use mod_genconf, only: typ_genconf
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_genconf), intent(in):: genconf
end subroutine genconf_diatomic
! ./src/genconf_random.F90 :
subroutine genrandom(parini,genconf)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, typ_file_info
    use mod_genconf, only: typ_genconf
    implicit none
    type(typ_parini), intent(inout):: parini
    type(typ_genconf), intent(in):: genconf
end subroutine genrandom
! ./src/genconf_rangrow.F90 :
subroutine rangrow(parini,genconf)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, typ_file_info
    use mod_genconf, only: typ_genconf
    implicit none
    type(typ_parini), intent(inout):: parini
    type(typ_genconf), intent(in):: genconf
end subroutine rangrow
! ./src/genconf_trimer.F90 :
subroutine genconf_trimer(parini,genconf)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms_all, typ_file_info
    use mod_genconf, only: typ_genconf
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_genconf), intent(in):: genconf
end subroutine genconf_trimer
! ./src/hung.F90 :
subroutine hung(N,A,F,Z)
    implicit none
      integer:: n
      real(8)::  A(n,n),Z,U(n),V(n)
      integer F(N),FB(n), RC(n)
end subroutine hung
subroutine INCR_inalborz(n,F,J,FB,RC)
implicit none
      integer:: n,I,J,JJ,  F(n),FB(n),RC(n)
end subroutine incr_inalborz
subroutine INIT_inalborz(N,A,F,M,U,V,FB,P)
    implicit none
      integer:: n,m, F(n),FB(n),P(n)
      real(8) A(n,n) , U(n),V(n)
end subroutine init_inalborz
subroutine PATH_inalborz(N,A,II,F,JJ,U,V,FB,RC)
      implicit none
      integer:: N 
      real(8)::  A(n,n),U(n),V(N),PI(n), IA, MIN
      integer:: F(N),LR(n),UC(n)
      integer:: FB(n),RC(n)
      integer::  i,j,k,L,ii,jj,NUC,NLR,R
end subroutine path_inalborz
! ./src/io_acf.F90 :
subroutine acf_write(file_info,atoms,atoms_all,strkey)
    use mod_atoms, only: typ_file_info, typ_atoms, typ_atoms_all
    implicit none
    type(typ_file_info), intent(inout):: file_info
    type(typ_atoms), optional, intent(in):: atoms
    type(typ_atoms_all), optional, intent(in):: atoms_all
    character(*), optional, intent(in):: strkey
end subroutine acf_write
subroutine acf_write_new(file_info,atoms_arr,strkey)
    use mod_atoms, only: typ_file_info, typ_atoms, typ_atoms_arr
    implicit none
    type(typ_file_info), intent(inout):: file_info
    type(typ_atoms_arr),intent(in):: atoms_arr
    character(*), optional, intent(in):: strkey
end subroutine acf_write_new
subroutine rotate4acf(nat,rat,cv,cvrot)
    implicit none
    integer, intent(in):: nat
    real(8), intent(inout):: rat(3,nat)
    real(8), intent(in):: cv(3,3), cvrot(3,3)
end subroutine rotate4acf
subroutine acf_force_write(file_info,atoms,atoms_all,strkey)
    use mod_atoms, only: typ_file_info, typ_atoms, typ_atoms_all
    implicit none
    type(typ_file_info), intent(inout):: file_info
    type(typ_atoms), optional, intent(in):: atoms
    type(typ_atoms_all), optional, intent(in):: atoms_all
    character(*), optional, intent(in):: strkey
end subroutine acf_force_write
subroutine acf_read(parini,filename,nconfmax,atoms,atoms_all)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, typ_atoms_all
    implicit none
    type(typ_parini), intent(in):: parini
    character(*), intent(in):: filename
    integer, intent(in):: nconfmax
    type(typ_atoms), optional, intent(inout):: atoms
    type(typ_atoms_all), optional, intent(inout):: atoms_all
end subroutine acf_read
subroutine acf_read_new(parini,filename,nconfmax,atoms_arr)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, typ_atoms_arr
    implicit none
    type(typ_parini), intent(in):: parini
    character(*), intent(in):: filename
    integer, intent(in):: nconfmax
    type(typ_atoms_arr), intent(inout):: atoms_arr
end subroutine acf_read_new
subroutine str_parse(str,line,atoms_t,c5)
    use mod_atoms, only: typ_atoms, typ_atoms_all
    implicit none
    character(*), intent(in):: str
    integer, intent(in):: line
    type(typ_atoms), intent(inout):: atoms_t
    character(*), intent(inout):: c5
end subroutine str_parse
subroutine str_motion2bemoved(str_motion,bemoved)
    implicit none
    character(*), intent(in):: str_motion
    logical, intent(inout):: bemoved(3)
end subroutine str_motion2bemoved
! ./src/io_cube.F90 :
subroutine cube_read(filename,atoms,poisson)
    use mod_atoms, only: typ_atoms
    use mod_electrostatics, only: typ_poisson
    implicit none
    character(*), intent(in):: filename
    type(typ_atoms), intent(out):: atoms
    type(typ_poisson), intent(out):: poisson
end subroutine cube_read
subroutine cube_write(filename,atoms,poisson,rho_or_pot)
    use mod_atoms, only: typ_atoms
    use mod_electrostatics, only: typ_poisson
    implicit none
    character(*), intent(in):: filename
    type(typ_atoms), intent(in):: atoms
    type(typ_poisson), intent(in):: poisson
    character(*), intent(in):: rho_or_pot
end subroutine cube_write
! ./src/io_vasp.F90 :
subroutine write_poscar(filename,nat,rat,latvec,ntypat,natarr,comment,vasp5,comment2,atom_motion)
    implicit none
    integer:: nat, ntypat, natarr(128)
    character(*):: filename
    real(8):: rat(3,nat),latvec(3,3)
    character(*):: comment, comment2
    logical:: vasp5, atom_motion(3,nat)
end subroutine write_poscar
! ./src/io_xyz.F90 :
subroutine writexyz(filename,fn_position,nat,rat,bemoved,sat,cellvec,boundcond,comment)
    implicit none
    integer, intent(in):: nat
    real(8), intent(in):: rat(3,nat)
    logical, intent(in):: bemoved(3,nat)
    character(5), intent(in):: sat(nat)
    real(8), intent(in):: cellvec(3,3)
    character(*), intent(in):: filename, fn_position, boundcond, comment
end subroutine writexyz
subroutine readxyznat(filename,nat)
    implicit none
    integer:: nat
    character(*):: filename
end subroutine readxyznat
subroutine readxyz(filename,nat,rat,sat,comment1,comment2,atom_motion)
    implicit none
    integer:: nat
    real(8):: rat(3,nat)
    character(5):: sat(nat)
    character(*):: filename, comment1, comment2
    logical:: atom_motion(3,nat)
end subroutine readxyz
! ./src/lenosky_tightbinding.F90 :
subroutine lenoskytb_alborz(atoms,natsi,count_md)
    use mod_atoms, only: typ_atoms
    use mod_potl, only: potl_typ
    use mod_tightbinding, only: typ_partb, lenosky
    implicit none
    type(typ_atoms), intent(inout):: atoms
    integer, intent(in):: natsi
    real(8), intent(inout):: count_md
end subroutine lenoskytb_alborz
subroutine lenoskytb_init(partb,atoms,natsi)
    use mod_tightbinding, only: typ_partb, lenosky
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_partb), intent(inout):: partb
    type(typ_atoms), intent(in):: atoms
    integer, intent(in):: natsi
end subroutine lenoskytb_init
subroutine totalenergy(partb,atoms,natsi,pplocal)
    use mod_atoms, only: typ_atoms
    use mod_potl, only: potl_typ
    use mod_tightbinding, only: typ_partb, lenosky
    implicit none
    type(typ_partb), intent(inout):: partb
    type(typ_atoms), intent(inout):: atoms
    integer, intent(in):: natsi
    type(potl_typ), intent(inout):: pplocal
end subroutine totalenergy
subroutine pairenergy(partb,atoms,pplocal,natsi)
    use mod_tightbinding, only: typ_partb
    use mod_atoms, only: typ_atoms
    use mod_potl, only: potl_typ
    use mod_frame, only: clsframepp_type
    implicit none
    type(typ_partb), intent(inout):: partb
    type(typ_atoms), intent(inout):: atoms
    integer, intent(in):: natsi
    type(potl_typ), intent(in):: pplocal
end subroutine pairenergy
subroutine lenoskytb_final(partb)
    use mod_tightbinding, only: typ_partb
    implicit none
    type(typ_partb), intent(inout):: partb
end subroutine lenoskytb_final
subroutine VECT_SUBTRACT_F90(a,b,c) 
    implicit none
    real(8), intent(in):: a(3), b(3)
    real(8), intent(out):: c(3)
end subroutine vect_subtract_f90
subroutine APPLY_PBC_F90(a)
    implicit none
    real(8), intent(inout) :: a(3)
end subroutine apply_pbc_f90
subroutine CELLDIST_F90(p1,p2,a,d)
    implicit none
    real(8), intent(in):: p1(3), p2(3), a(3,3)
    real(8), intent(out):: d
end subroutine celldist_f90
subroutine CELLGRAD_F90(p1,p2,a,g)
    implicit none
    real(8), intent(in):: p1(3), p2(3), a(3,3)
    real(8), intent(out):: g(3)
end subroutine cellgrad_f90
subroutine radelmgeneralsp(r,radar,dradar,atomtypei,atomtypej,pplocal)
    use mod_potl, only: potl_typ
    implicit none
    type(potl_typ), intent(in):: pplocal
    real(8), intent(in):: r
    real(8), intent(out):: radar(0:3), dradar(0:3)
    integer, intent(in):: atomtypei, atomtypej
end subroutine radelmgeneralsp
subroutine clssplint(s,x,y,deriv,extype)
    use mod_splinetb, only: NSPMAX, spline_typ
    implicit none
    type(spline_typ), intent(in) :: s
    integer, intent(in)::  extype 
    real(8), intent(in) :: x
    real(8), intent(out) :: y, deriv
end subroutine clssplint
subroutine eselfgeneral(eself)
    implicit none
    real(8), intent(inout):: eself(0:3)
end subroutine eselfgeneral
subroutine prmst38c(partb,pplocal)
    use mod_tightbinding, only: typ_partb
    use mod_potl, only: potl_typ
    implicit none
    type(typ_partb):: partb
    type(potl_typ):: pplocal
end subroutine prmst38c
subroutine clsfread_spline(unit,s)
    use mod_splinetb, only: NSPMAX, spline_typ
    implicit none
    integer:: unit
    type(spline_typ), intent(out):: s
end subroutine clsfread_spline
subroutine clsspline(s)
    use mod_splinetb, only: NSPMAX, spline_typ
    implicit none
    type(spline_typ), intent(inout):: s
end subroutine clsspline
! ./src/lmder_modified.F90 :
subroutine lmder_modified(fcn,m,n,x,fvec,fjac,ldfjac,ftol,xtol,gtol,maxfev,diag,mode,factor,nprint,info,nfev,njev,ipvt,qtf,wa1,wa2,wa3,wa4,parini,ann_arr,atoms_train,atoms_valid,symfunc_train,symfunc_valid,ekf)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr, typ_symfunc_arr, typ_ekf
    use mod_atoms, only: typ_atoms_arr
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_atoms_arr), intent(in):: atoms_train, atoms_valid
    type(typ_symfunc_arr), intent(in):: symfunc_train, symfunc_valid
    type(typ_ekf), intent(inout):: ekf
      external fcn
      integer:: m,n,ldfjac,maxfev,mode,nprint,info,nfev,njev
      integer:: ipvt(n)
      real(8):: ftol,xtol,gtol,factor
      real(8):: x(n),fvec(m),fjac(ldfjac,n),diag(n),qtf(n), wa1(n),wa2(n),wa3(n),wa4(m)
end subroutine lmder_modified
! ./src/md.F90 :
subroutine dynamics(parini)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, typ_file_info
    implicit none
    type(typ_parini), intent(inout):: parini
end subroutine dynamics
subroutine md_nve(parini,atoms)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, typ_file_info
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms):: atoms, atoms_old
end subroutine md_nve
subroutine md_nph(parini,atoms)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, typ_file_info
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms):: atoms
end subroutine md_nph
subroutine set_velocities(atoms, ekin_arg)
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_atoms), intent(inout):: atoms
    real(8), optional::ekin_arg
end subroutine set_velocities
subroutine ekin_temprature(atoms,temp,vcm,rcm,totmass) 
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_atoms):: atoms
    real(8) :: vcm(3), rcm(3), t1, tmp
    real(8) :: aboltzmann, totmass, temp 
end subroutine ekin_temprature
! ./src/md_NVT.F90 :
subroutine md_nvt_langevin(parini,atoms)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, typ_file_info
    implicit none
    type(typ_parini), intent(inout):: parini
    type(typ_atoms):: atoms
    real(8):: eta(3,atoms%nat)
    real(8):: langev(atoms%nat), forces_langevin(3,atoms%nat)
    real(8):: rat_next(3,atoms%nat), vat_old(3,atoms%nat)
end subroutine md_nvt_langevin
subroutine md_nvt_nose_hoover_cp(parini,atoms)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, typ_file_info
    implicit none
    type(typ_parini), intent(inout):: parini
    type(typ_atoms):: atoms
    real(8):: forces_nosehoover(3,atoms%nat)
    real(8):: rat_next(3,atoms%nat), vat_old(3,atoms%nat)
end subroutine md_nvt_nose_hoover_cp
subroutine md_nvt_nose_hoover_chain(parini,atoms)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, typ_file_info
    implicit none
    type(typ_parini), intent(inout):: parini
    type(typ_atoms):: atoms
    real(8):: forces_nosehoover(3,atoms%nat)
    real(8):: rat_next(3,atoms%nat), vat_old(3,atoms%nat)
end subroutine md_nvt_nose_hoover_chain
subroutine set_langevin_randforce(eta,nat)
    implicit none
    integer :: nat
    real(8) ::eta(3,nat), sum1, sum2, sum3
end subroutine set_langevin_randforce
subroutine back_to_cell(atoms)
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_atoms):: atoms
end subroutine back_to_cell
subroutine plane_repulsion(atoms)
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_atoms):: atoms
end subroutine plane_repulsion
subroutine thermostat_evolution(atoms,zeta_next,zeta,zeta_prev,dzeta,mass_q,kt,ntherm,imd)
    use mod_atoms, only: typ_atoms, typ_file_info
    implicit none
    type(typ_atoms):: atoms
    integer :: ntherm, imd 
    real(8) :: kt, t1
    real(8) :: zeta_next(3,atoms%nat,ntherm), zeta(3,atoms%nat,ntherm),zeta_prev(3,atoms%nat,ntherm)
    real(8) :: dzeta(3,atoms%nat,ntherm), mass_q(ntherm)
    real(8) :: force_therm(3,atoms%nat,ntherm)
end subroutine thermostat_evolution
subroutine thermostat_evolution_2(atoms,zeta_next,zeta,zeta_prev,dzeta,mass_q,kt,ntherm,imd)
    use mod_atoms, only: typ_atoms, typ_file_info
    implicit none
    type(typ_atoms):: atoms
    integer :: ntherm, imd 
    real(8) :: kt, t1, tt
    real(8) :: zeta_next(ntherm), zeta(ntherm),zeta_prev(ntherm)
    real(8) :: dzeta(ntherm), mass_q(ntherm)
    real(8) :: force_therm(ntherm)
end subroutine thermostat_evolution_2
subroutine md_nvt_nose_hoover(parini,atoms)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, typ_file_info
    implicit none
    type(typ_parini), intent(inout):: parini
    type(typ_atoms):: atoms
    real(8):: eta(3,atoms%nat)
    real(8)::  forces_nose(3,atoms%nat)
    real(8):: rat_next(3,atoms%nat), vat_old(3,atoms%nat)
end subroutine md_nvt_nose_hoover
subroutine get_atomic_mass(atoms,totmass)
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_atoms):: atoms
    real(8):: totmass
end subroutine get_atomic_mass
! ./src/minhopp_allocation.F90 :
subroutine allocate_minhopp_arrays1(nproc)
    implicit none
    integer, intent(in):: nproc
end subroutine allocate_minhopp_arrays1
subroutine allocate_minhopp_arrays2(nat,nproc)
    implicit none
    integer, intent(in):: nat, nproc
end subroutine allocate_minhopp_arrays2
subroutine deallocate_minhopp_arrays
    implicit none
end subroutine deallocate_minhopp_arrays
! ./src/minhopp.F90 :
subroutine minimahopping(parini)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, typ_atoms_arr
    use mod_opt, only: typ_paropt
    implicit none
    type(typ_parini), intent(in):: parini
end subroutine minimahopping
subroutine init_minimahopping(parini,atoms_curr,atoms_hopp,atoms_allproc,atoms_locmin,paropt,paropt_prec)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, typ_atoms_arr
    use mod_opt, only: typ_paropt
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_paropt), intent(inout):: paropt_prec, paropt
    type(typ_atoms), intent(inout):: atoms_curr, atoms_hopp
    type(typ_atoms_arr), intent(inout):: atoms_allproc
    type(typ_atoms_arr), intent(inout):: atoms_locmin
end subroutine init_minimahopping
subroutine final_minimahopping(parini,atoms_curr,atoms_hopp,atoms_allproc,atoms_locmin,paropt,paropt_prec)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, typ_atoms_arr
    use mod_opt, only: typ_paropt
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_paropt), intent(inout):: paropt_prec, paropt
    type(typ_atoms), intent(inout):: atoms_curr, atoms_hopp
    type(typ_atoms_arr), intent(inout):: atoms_allproc
    type(typ_atoms_arr), intent(inout):: atoms_locmin
end subroutine final_minimahopping
subroutine set_amass(atoms_hopp)
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_atoms), intent(inout):: atoms_hopp
end subroutine set_amass
subroutine relax_minhopp(parini,atoms,paropt_prec,paropt)
    use mod_parini, only: typ_parini
    use mod_opt, only: typ_paropt
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
    type(typ_paropt), intent(inout):: paropt_prec, paropt
end subroutine relax_minhopp
subroutine print_minhopp_parameters
    implicit none
end subroutine print_minhopp_parameters
subroutine read_earr
    implicit none
end subroutine read_earr
subroutine readnat(atoms_curr)
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_atoms):: atoms_curr
end subroutine readnat
subroutine read_poscur(atoms_curr,atoms_allproc)
    use mod_atoms, only: typ_atoms, typ_atoms_arr
    implicit none
    type(typ_atoms):: atoms_curr
    type(typ_atoms_arr):: atoms_allproc
end subroutine read_poscur
subroutine read_minhopp_parameters 
    implicit none
end subroutine read_minhopp_parameters
subroutine minhopp_newrun_initialization(atoms_curr,atoms_locmin)
    use mod_atoms, only: typ_atoms, typ_atoms_arr
    implicit none
    type(typ_atoms), intent(in):: atoms_curr
    type(typ_atoms_arr), intent(inout):: atoms_locmin
end subroutine minhopp_newrun_initialization
subroutine read_poslow(atoms_locmin)
    use mod_atoms, only: typ_atoms_arr
    implicit none
    type(typ_atoms_arr), intent(inout):: atoms_locmin
end subroutine read_poslow
subroutine send_minimum_to_all(atoms_curr)
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_atoms), intent(in):: atoms_curr
end subroutine send_minimum_to_all
subroutine send_minhopp_parameters_to_all(atoms_curr)
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_atoms), intent(in):: atoms_curr
end subroutine send_minhopp_parameters_to_all
subroutine mdescape(parini,atoms_hopp)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, typ_file_info
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms_hopp
end subroutine mdescape
subroutine collect_data_from_all_processors(ntry,atoms_curr,atoms_allproc,atoms_locmin)
    use mod_atoms, only: typ_atoms, typ_atoms_arr
    implicit none
    integer, intent(in):: ntry
    type(typ_atoms), intent(in):: atoms_curr
    type(typ_atoms_arr), intent(inout):: atoms_allproc, atoms_locmin
end subroutine collect_data_from_all_processors
subroutine request_receive(atoms_allproc)
    use mod_atoms, only: typ_atoms_arr
    implicit none
    type(typ_atoms_arr), intent(in):: atoms_allproc
end subroutine request_receive
subroutine test_receive(atoms_allproc,atoms_locmin)
    use mod_atoms, only: typ_atoms_arr
    implicit none
    type(typ_atoms_arr), intent(inout):: atoms_allproc, atoms_locmin
end subroutine test_receive
subroutine cancel_excessive_irecv
    implicit none
end subroutine cancel_excessive_irecv
subroutine insert_alborz(kepos,epos)
    implicit none
    integer, intent(in):: kepos
    real(8), intent(in):: epos
end subroutine insert_alborz
subroutine save_low_conf_alborz(atoms,atoms_locmin)
    use mod_atoms, only: typ_atoms, typ_atoms_arr
    implicit none
    type(typ_atoms), intent(in):: atoms
    type(typ_atoms_arr), intent(inout):: atoms_locmin
end subroutine save_low_conf_alborz
subroutine velopt(parini,atoms)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
end subroutine velopt
subroutine soften(parini,nstep,atoms0,count_soften,count_soften_tot)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, typ_file_info
    implicit none
    type(typ_parini), intent(in):: parini
    integer, intent(in):: nstep
    type(typ_atoms), intent(inout):: atoms0
    real(8), intent(inout):: count_soften, count_soften_tot
end subroutine soften
subroutine write_minhopp(atoms_allproc,atoms_locmin)
    use mod_atoms, only: typ_atoms_arr, typ_file_info
    implicit none
    type(typ_atoms_arr), intent(inout):: atoms_allproc, atoms_locmin
end subroutine write_minhopp
subroutine write_minhopp_parameters
    implicit none
end subroutine write_minhopp_parameters
subroutine write_earr
    implicit none
end subroutine write_earr
subroutine escape_failed(parini,erat,erathopp)
    use mod_parini, only: typ_parini
    implicit none
    type(typ_parini), intent(in):: parini
    real(8), intent(in):: erat, erathopp
end subroutine escape_failed
subroutine local_minimum_accepted(atoms_hopp,atoms_curr,atoms_locmin)
    use mod_atoms, only: typ_atoms, typ_atoms_arr
    implicit none
    type(typ_atoms), intent(in):: atoms_hopp
    type(typ_atoms), intent(inout):: atoms_curr
    type(typ_atoms_arr), intent(inout):: atoms_locmin
end subroutine local_minimum_accepted
subroutine local_minimum_rejected(atoms_hopp)
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_atoms), intent(in):: atoms_hopp
end subroutine local_minimum_rejected
subroutine already_visited_minimum(parini)
    use mod_parini, only: typ_parini
    implicit none
    type(typ_parini), intent(in):: parini
end subroutine already_visited_minimum
subroutine new_minimum(atoms_hopp)
    use mod_atoms, only: typ_atoms, typ_file_info
    implicit none
    type(typ_atoms), intent(in):: atoms_hopp
end subroutine new_minimum
subroutine print_final_statistics
    implicit none
end subroutine print_final_statistics
subroutine MPI_atom_arr_copy(nat,atoms_arr)
    use mod_atoms, only: typ_atoms, typ_atoms_arr
    implicit none
    integer,intent(in)::nat
    type(typ_atoms_arr),intent(inout):: atoms_arr
end subroutine mpi_atom_arr_copy
! ./src/minhopp_pot.F90 :
subroutine setpot_init(parini,atoms_curr,paropt,paropt_prec)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_opt, only: typ_paropt
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms_curr
    type(typ_paropt), intent(inout):: paropt, paropt_prec
end subroutine setpot_init
subroutine setpot_final(parini,atoms_curr)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(in):: atoms_curr
end subroutine setpot_final
subroutine setpot_mdescape
    implicit none
end subroutine setpot_mdescape
subroutine setpot_soften
    implicit none
end subroutine setpot_soften
subroutine setpot_geopt_prec
    implicit none
end subroutine setpot_geopt_prec
subroutine setpot_geopt
    implicit none
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
    implicit none
    integer, intent(in):: iproc, nr, nwork
    real(8):: x(nr), f(nr), epot, work(nwork)
    type(typ_paropt), intent(inout):: paropt
end subroutine mybfgs
subroutine init_mybfgs(paropt,epot,fmax)
    use mod_opt, only: typ_paropt
    implicit none
    type(typ_paropt), intent(inout):: paropt
    real(8), intent(in):: epot, fmax
end subroutine init_mybfgs
! ./src/optimizer_cg.F90 :
subroutine cgminimum(iproc,n,nr,x,f,epot,paropt,nwork,work)
    use mod_opt, only: typ_paropt, frmt_base
    implicit none
    type(typ_paropt):: paropt
    integer, intent(in):: iproc, n, nr, nwork
    real(8):: x(n), f(n), epot, work(nwork)
end subroutine cgminimum
subroutine init_cgminimum(paropt,n,nr,f,nwork,work,epot,fnrm)
    use mod_opt, only: typ_paropt
    implicit none
    type(typ_paropt), intent(inout):: paropt
    integer, intent(in):: n, nr, nwork
    real(8), intent(in):: f(n), epot, fnrm
    real(8), intent(inout):: work(nwork)
end subroutine init_cgminimum
! ./src/optimizer_dfp.F90 :
subroutine mydfp(nr,x,epot,f,nwork,work,paropt)
    use mod_opt, only: typ_paropt
    implicit none
    integer::nr,nwork,mf,my,ms,nrsq,iw1,iw2,iw3,info,i,j,l,mx
    real(8)::x(nr),f(nr),epot,work(nwork)
    type(typ_paropt)::paropt
end subroutine mydfp
! ./src/optimizer_diis.F90 :
subroutine diisminimum(n,nr,x,epot,f,paropt,nwork,work)
    use mod_opt, only: typ_paropt
    implicit none
    integer:: n, nr, nwork, info, id, jd
    real(8):: x(n), f(n), epot, work(nwork), fnrm, dnrm2, ddot, fmax
    type(typ_paropt):: paropt
end subroutine diisminimum
! ./src/optimizer_drivers.F90 :
subroutine minimize(parini,iproc,atoms,paropt)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, typ_file_info
    use mod_opt, only: typ_paropt
    implicit none
    type(typ_parini), intent(in):: parini
    integer, intent(in):: iproc
    type(typ_atoms), intent(inout):: atoms
    type(typ_paropt), intent(inout):: paropt
end subroutine minimize
subroutine test_convergence(n,f,paropt)
    use mod_opt, only: typ_paropt
    implicit none
    integer, intent(in):: n
    real(8), intent(in):: f(n)
    type(typ_paropt), intent(inout):: paropt
end subroutine test_convergence
subroutine x_to_xr(n,x,f,bemoved,nr,xr,fr)
    implicit none
    integer, intent(in):: n, nr
    real(8), intent(in):: x(n), f(n)
    logical, intent(in):: bemoved(3,n/3)
    real(8), intent(inout):: xr(nr), fr(nr)
end subroutine x_to_xr
subroutine xr_to_x(nr,xr,n,bemoved,x)
    implicit none
    integer, intent(in):: n, nr
    real(8), intent(in):: xr(nr)
    logical, intent(in):: bemoved(3,n/3)
    real(8), intent(inout):: x(n)
end subroutine xr_to_x
subroutine report_param(paropt)
    use mod_opt, only: typ_paropt
    implicit none
    type(typ_paropt), intent(inout):: paropt
end subroutine report_param
subroutine initminimize(paropt)
    use mod_opt, only: typ_paropt
    implicit none
    type(typ_paropt), intent(inout):: paropt
end subroutine initminimize
subroutine finalminimize(paropt)
    use mod_opt, only: typ_paropt
    implicit none
    type(typ_paropt), intent(inout):: paropt
end subroutine finalminimize
! ./src/optimizer_drivers_vc.F90 :
subroutine vc_minimize(parini,iproc,atoms,paropt)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, typ_file_info
    use mod_opt, only: typ_paropt
    implicit none
    type(typ_parini), intent(in):: parini
    integer, intent(in):: iproc !, n, nr
    type(typ_atoms), intent(inout):: atoms
    type(typ_paropt), intent(inout):: paropt
end subroutine vc_minimize
subroutine vc_test_convergence(n,f,paropt)
    use mod_opt, only: typ_paropt
    implicit none
    integer, intent(in):: n
    real(8), intent(in):: f(n)
    type(typ_paropt), intent(inout):: paropt
end subroutine vc_test_convergence
subroutine vc_x_to_xr(atoms,nr,xr,fr)
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_atoms), intent(in):: atoms
    integer, intent(in):: nr
    real(8), intent(inout):: xr(nr), fr(nr)
end subroutine vc_x_to_xr
subroutine vc_xr_to_x(nr,xr,atoms)
    use mod_atoms, only: typ_atoms
    implicit none
    integer, intent(in):: nr
    real(8), intent(in):: xr(nr)
    type(typ_atoms), intent(inout):: atoms
end subroutine vc_xr_to_x
subroutine vc_report_param(paropt)
    use mod_opt, only: typ_paropt
    implicit none
    type(typ_paropt), intent(inout):: paropt
end subroutine vc_report_param
! ./src/optimizer_fire.F90 :
subroutine fire(parini,iproc,n,x,epot,f,work,paropt)
    use mod_parini, only: typ_parini
    use mod_opt, only: typ_paropt, frmt_base
    implicit none
    type(typ_parini), intent(in):: parini
    integer, intent(in):: n, iproc
    real(8), intent(inout):: x(n), epot, f(n)
    real(8), intent(inout):: work(3*n) !1:n velocities, n+1:2*n previous force
    type(typ_paropt), intent(inout):: paropt
end subroutine fire
subroutine init_fire(n,f,epot,work,paropt)
    use mod_opt, only: typ_paropt
    implicit none
    integer, intent(in):: n
    real(8), intent(in):: epot, f(n)
    real(8), intent(inout):: work(3*n)
    type(typ_paropt), intent(inout):: paropt
end subroutine init_fire
! ./src/optimizer_gmdfire.F90 :
subroutine gmdfire(nr,x,epot,f,work,paropt)
    use mod_opt, only: typ_paropt
    implicit none
    integer:: nr
    real(8):: x(nr), epot, f(nr), de, DDOT, fnrm, fmax, vnrm, dt, p
    real(8):: work(5*nr) !1:nr velocities, nr+1:2*nr previous force
    type(typ_paropt)::paropt
end subroutine gmdfire
! ./src/optimizer_nlbfgs.F90 :
! ./src/optimizer_sd.F90 :
subroutine sdminimum(parini,iproc,nr,x,f,epot,paropt,nwork,work)
    use mod_parini, only: typ_parini
    use mod_opt, only: typ_paropt, frmt_base
    implicit none
    type(typ_parini), intent(in):: parini
    integer, intent(in):: iproc, nr, nwork
    real(8), intent(inout):: x(nr), f(nr), work(nwork)
    real(8), intent(in):: epot
    type(typ_paropt), intent(inout):: paropt
end subroutine sdminimum
subroutine init_sdminimum(paropt,nr,x,nwork,work)
    use mod_opt, only: typ_paropt
    implicit none
    type(typ_paropt), intent(inout):: paropt
    integer, intent(in):: nr, nwork
    real(8), intent(in):: x(nr)
    real(8), intent(inout):: work(nwork)
end subroutine init_sdminimum
subroutine what_is_condition_of_feedback(paropt,de1,df1,feedbackcondition)
    use mod_opt, only: typ_paropt
    implicit none
    type(typ_paropt), intent(in):: paropt
    real(8), intent(in):: de1, df1
    logical, intent(out):: feedbackcondition
end subroutine what_is_condition_of_feedback
subroutine test_saturation(paropt,de1,de2,df2,fnrm)
    use mod_opt, only: typ_paropt
    implicit none
    type(typ_paropt), intent(inout):: paropt
    real(8), intent(in):: de1, de2, df2, fnrm !, fnrmitm1
end subroutine test_saturation
subroutine final_sdminimum(paropt)
    use mod_opt, only: typ_paropt
    implicit none
    type(typ_paropt), intent(inout):: paropt
end subroutine final_sdminimum
! ./src/optimizer_sqnm.F90 :
subroutine sqnm(parini,atoms,paropt,count_sqnm,fail)
   use mod_parini, only: typ_parini
   use mod_atoms, only: typ_atoms
   use mod_opt, only: typ_paropt, frmt_base
   implicit none
   type(typ_parini), intent(in):: parini
   type(typ_atoms), intent(inout):: atoms
   type(typ_paropt), intent(inout):: paropt
   real(8), intent(inout):: count_sqnm
   logical, intent(out):: fail
   integer :: nat    !< number of atoms
   real(8) :: trustr !< a single atoms is not allowed to be dsiplaced more than by trustr
end subroutine sqnm
subroutine minenergyandforces_alborz(parini,iproc,nproc,eeval,imode,atoms,nat,rat,fat,fstretch,fxyzraw,epot,iconnect,nbond_,wold,alpha_stretch0,alpha_stretch,infocode)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    implicit none
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
    implicit none
    integer, intent(in) :: iproc
    type(typ_atoms), intent(in) :: atoms
    real(8), intent(out) :: rcov(atoms%nat)
end subroutine give_rcov_sqnm
subroutine getSubSpaceEvecEval(label,iproc,verbosity,nat,nhist,nhistx,ndim,cutoffratio,lwork,work,idx,rxyz,fxyz,aa,rr,ff,rrr,fff,eval,res,success)
    implicit none
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
    implicit none
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
    use mod_atoms, only: typ_atoms
    implicit none
    integer, intent(in) :: iproc,verbosity
    character(len=*), intent(in) :: label
    type(typ_atoms), intent(in) :: atoms
    real(8), intent(in) :: rcov(atoms%nat)
    integer, intent(out) :: nbond
    integer, intent(out) :: iconnect(2,1000)
end subroutine findbonds
subroutine projectbond(nat,nbond,rat,fat,fstretch,iconnect,wold,alpha_stretch0,alpha_stretch)
    implicit none
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
! ./src/parser_all.F90 :
subroutine get_main_parameters(file_ini,parini)
    use mod_task, only: typ_file_ini
    use mod_parini, only: typ_parini
    implicit none
    type(typ_file_ini), intent(inout):: file_ini
    type(typ_parini), intent(inout):: parini
end subroutine get_main_parameters
subroutine get_minhopp_parameters(file_ini,parini)
    use mod_parini, only: typ_parini
    use mod_task, only: typ_file_ini
    implicit none
    type(typ_file_ini), intent(inout):: file_ini
    type(typ_parini), intent(inout):: parini
end subroutine get_minhopp_parameters
subroutine get_opt_param(file_ini,paropt)
    use mod_task, only: typ_file_ini
    use mod_opt, only: typ_paropt
    implicit none
    type(typ_file_ini), intent(inout):: file_ini
    type(typ_paropt), intent(inout):: paropt
end subroutine get_opt_param
subroutine get_geopt_parameters(file_ini,parini)
    use mod_parini, only: typ_parini
    use mod_task, only: typ_file_ini
    implicit none
    type(typ_file_ini), intent(inout):: file_ini
    type(typ_parini), intent(inout):: parini
end subroutine get_geopt_parameters
subroutine get_geopt_prec_parameters(file_ini,parini)
    use mod_parini, only: typ_parini
    use mod_task, only: typ_file_ini
    implicit none
    type(typ_file_ini), intent(inout):: file_ini
    type(typ_parini), intent(inout):: parini
end subroutine get_geopt_prec_parameters
subroutine get_saddle_1s_opt_parameters(file_ini,parini)
    use mod_parini, only: typ_parini
    use mod_task, only: typ_file_ini
    implicit none
    type(typ_file_ini), intent(inout):: file_ini
    type(typ_parini), intent(inout):: parini
end subroutine get_saddle_1s_opt_parameters
subroutine get_saddle_1s_parameters(file_ini,parini)
    use mod_parini, only: typ_parini
    use mod_task, only: typ_file_ini
    implicit none
    type(typ_file_ini), intent(inout):: file_ini
    type(typ_parini), intent(inout):: parini
end subroutine get_saddle_1s_parameters
subroutine get_potential_parameters(file_ini,parini)
    use mod_parini, only: typ_parini
    use mod_task, only: typ_file_ini
    implicit none
    type(typ_file_ini), intent(inout):: file_ini
    type(typ_parini), intent(inout):: parini
end subroutine get_potential_parameters
subroutine get_ann_parameters(file_ini,parini)
    use mod_parini, only: typ_parini
    use mod_task, only: typ_file_ini
    implicit none
    type(typ_file_ini), intent(inout):: file_ini
    type(typ_parini), intent(inout):: parini
end subroutine get_ann_parameters
subroutine get_dynamics_parameters(file_ini,parini)
    use mod_parini, only: typ_parini
    use mod_task, only: typ_file_ini
    implicit none
    type(typ_file_ini), intent(inout):: file_ini
    type(typ_parini), intent(inout):: parini
end subroutine get_dynamics_parameters
subroutine get_genconf_parameters(file_ini,parini)
    use mod_parini, only: typ_parini
    use mod_task, only: typ_file_ini
    implicit none
    type(typ_file_ini), intent(inout):: file_ini
    type(typ_parini), intent(inout):: parini
end subroutine get_genconf_parameters
subroutine get_conf_comp_parameters(file_ini,parini)
    use mod_parini, only: typ_parini
    use mod_task, only: typ_file_ini
    implicit none
    type(typ_file_ini), intent(inout):: file_ini
    type(typ_parini), intent(inout):: parini
end subroutine get_conf_comp_parameters
subroutine get_testforces_parameters(file_ini,parini)
    use mod_parini, only: typ_parini
    use mod_task, only: typ_file_ini
    implicit none
    type(typ_file_ini), intent(inout):: file_ini
    type(typ_parini), intent(inout):: parini
end subroutine get_testforces_parameters
subroutine get_single_point_parameters(file_ini,parini)
    use mod_parini, only: typ_parini
    use mod_task, only: typ_file_ini
    implicit none
    type(typ_file_ini), intent(inout):: file_ini
    type(typ_parini), intent(inout):: parini
end subroutine get_single_point_parameters
subroutine get_ewald_parameters(file_ini,parini)
    use mod_parini, only: typ_parini
    use mod_task, only: typ_file_ini
    implicit none
    type(typ_file_ini), intent(inout):: file_ini
    type(typ_parini), intent(inout):: parini
end subroutine get_ewald_parameters
subroutine get_misc_parameters(file_ini,parini)
    use mod_parini, only: typ_parini
    use mod_task, only: typ_file_ini
    implicit none
    type(typ_file_ini), intent(inout):: file_ini
    type(typ_parini), intent(inout):: parini
end subroutine get_misc_parameters
! ./src/parser_core.F90 :
subroutine read_file_input(file_ini)
    use mod_task, only: typ_file_ini
    implicit none
    type(typ_file_ini), intent(inout):: file_ini
end subroutine read_file_input
subroutine get_header_location(file_ini,str_header)
    use mod_task, only: typ_file_ini
    implicit none
    type(typ_file_ini), intent(inout):: file_ini
    character(*), intent(in):: str_header
end subroutine get_header_location
subroutine split_line(file_ini)
    use mod_task, only: typ_file_ini
    implicit none
    type(typ_file_ini), intent(inout):: file_ini
end subroutine split_line
subroutine get_one_param(file_ini,var_name,int_var,real_var,char_var,char_line_var,log_var)
    use mod_task, only: typ_file_ini
    implicit none
    type(typ_file_ini), intent(inout):: file_ini
    character(*), intent(in):: var_name
    integer, optional, intent(out):: int_var
    real(8), optional, intent(out):: real_var
    character(*), optional, intent(out):: char_var
    character(*), optional, intent(out):: char_line_var
    logical, optional, intent(out):: log_var
end subroutine get_one_param
! ./src/phonon.F90 :
subroutine cal_hessian_4p(parini)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_parini), intent(in):: parini
end subroutine cal_hessian_4p
! ./src/plain_ewald.F90 :
subroutine plain_ewald(atoms,en)
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_atoms), intent(inout):: atoms
    real(8):: pi, q(1:atoms%nat)
    real(8):: k, k2, alpha, sum_en, en
end subroutine plain_ewald
subroutine structur_factor_comput(atoms,kx,ky,kz,s1,s2,sfactor_norm)
    use mod_atoms, only: typ_atoms
    implicit none
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
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
end subroutine init_potential_ann
subroutine cal_potential_ann(parini,atoms)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_ann, only: typ_symfunc
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
end subroutine cal_potential_ann
subroutine final_potential_ann
    implicit none
end subroutine final_potential_ann
subroutine add_repulsive_potential(parini,atoms)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_linked_lists, only: typ_linked_lists
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
end subroutine add_repulsive_potential
! ./src/potential_BigDFT.F90 :
subroutine init_potential_forces_bigdft(atoms_t)
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_atoms), intent(in):: atoms_t
end subroutine init_potential_forces_bigdft
subroutine final_potential_forces_bigdft
    implicit none
end subroutine final_potential_forces_bigdft
subroutine cal_potential_forces_bigdft(atoms)
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_atoms), intent(inout):: atoms
end subroutine cal_potential_forces_bigdft
subroutine writexyz_bigdft(filename,nat,rat,comment)
    implicit none
    integer:: nat,iat
    real(8)::rat(3,nat),x,y,z,cellx,celly,cellz
    character(*)::filename,comment
end subroutine writexyz_bigdft
subroutine get_output_bigdft(iproc,filename,nat,fat,epot,success)
    implicit none
    integer, intent(in):: iproc, nat
    character(*):: filename
    real(8):: fat(3,nat),epot
    logical:: success
end subroutine get_output_bigdft
! ./src/potential_BLJ_vc.F90 :
subroutine init_lennardjones_vc(nat,sat)
    implicit none
    integer, intent(in):: nat
    character(5), intent(in):: sat(nat)
end subroutine init_lennardjones_vc
subroutine lennardjones_vc(iproc,nat,xred0,latvec,pressure,fxyz,celldv,stress,etot,enth)
    implicit none
    integer, intent(in):: iproc, nat
    real(8):: xred(3,nat),fxyz(3,nat),xred0(3,nat),dxyz(3),r1red(3),r2red(3),rcut2(2,2)
    real(8):: etot,latvec_x(3,3),rec_nec(3),enth,stressvol(3,3),vol
    real(8):: latvec(3,3),latvecinv(3,3),celldv(3,3),pressure,gradvol(3,3),stress(3,3),tmplat(3,3)
end subroutine lennardjones_vc
subroutine stress_volume_alborz(latvec,vol,pressure,stressvol)
    implicit none
    real(8):: latvec(3,3),vol,stressvol(3,3),inv_latvec(3,3),pressure
end subroutine stress_volume_alborz
subroutine cell_vol(nat,latvec,vol)
    implicit none
    integer:: nat
    real(8):: latvec(3,3),vol,a(3,3)
end subroutine cell_vol
! ./src/potential_DFTB.F90 :
subroutine init_potential_forces_dftb(atoms)
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_atoms), intent(in):: atoms
end subroutine init_potential_forces_dftb
subroutine final_potential_forces_dftb
    implicit none
end subroutine final_potential_forces_dftb
subroutine cal_potential_forces_dftb(atoms)
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_atoms), intent(inout):: atoms
end subroutine cal_potential_forces_dftb
subroutine writexyz_dftb(filename,atoms)
    use mod_atoms, only: typ_atoms
    implicit none
    character(*):: filename !,comment
    type(typ_atoms), intent(in):: atoms
end subroutine writexyz_dftb
subroutine get_output_dftb(filename,atoms,success)
    use mod_atoms, only: typ_atoms
    implicit none
    character(*), intent(in):: filename
    type(typ_atoms), intent(inout):: atoms
    logical, intent(out):: success
end subroutine get_output_dftb
! ./src/potential_FF.F90 :
subroutine init_potential_forces_ff(parini,atoms)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
end subroutine init_potential_forces_ff
subroutine cal_potential_forces_ff(parini,atoms)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
end subroutine cal_potential_forces_ff
subroutine final_potential_forces_ff(parini,atoms)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(in):: atoms
end subroutine final_potential_forces_ff
! ./src/potential_LJ.F90 :
subroutine init_lennardjones
    implicit none
end subroutine init_lennardjones
subroutine lennardjones(atoms)
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_atoms), intent(inout):: atoms
end subroutine lennardjones
! ./src/potential_LTB.F90 :
subroutine init_lenosky_tb(atoms_t)
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_atoms), intent(in):: atoms_t
end subroutine init_lenosky_tb
subroutine lenosky_tb(atoms)
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_atoms), intent(inout):: atoms
end subroutine lenosky_tb
! ./src/potential_main.F90 :
subroutine init_potential_forces(parini,atoms)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
end subroutine init_potential_forces
subroutine cal_potential_forces(parini,atoms)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
end subroutine cal_potential_forces
subroutine final_potential_forces(parini,atoms)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(in):: atoms
end subroutine final_potential_forces
subroutine remove_drift(atoms)
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_atoms), intent(inout):: atoms
end subroutine remove_drift
! ./src/potential_main_vc.F90 :
subroutine vc_init_potential_forces(atoms)
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_atoms), intent(inout):: atoms
end subroutine vc_init_potential_forces
subroutine cal_potential_forces_vc(iproc,nat,rat,cellvec,pressure,fat,celldv,stress,epot,enth)
    implicit none
    integer, intent(in):: iproc, nat
    real(8), intent(in):: rat(3,nat), cellvec(3,3), pressure
    real(8), intent(inout):: fat(3,nat), celldv(3,3), stress(3,3), epot, enth
end subroutine cal_potential_forces_vc
! ./src/potential_MPMD.F90 :
subroutine init_potential_forces_mpmd(atoms_t)
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_atoms), intent(in):: atoms_t
end subroutine init_potential_forces_mpmd
subroutine cal_potential_forces_mpmd(atoms)
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_atoms), intent(inout):: atoms
end subroutine cal_potential_forces_mpmd
subroutine mpmd_init
    implicit none
end subroutine mpmd_init
! ./src/potential_NetSock.F90 :
  subroutine cal_potential_forces_netsock(atoms)
  use mod_potential, only: sock_socket, sock_inet, sock_port,sock_host,MSGLEN,sock_extra_string,sock_ecutwf,reset
  use mod_atoms, only: typ_atoms
  implicit none
  type(typ_atoms), intent(inout):: atoms
  real(8):: xred(3,atoms%nat)
  real(8):: fcart(3,atoms%nat),energy,strten(6)
end subroutine cal_potential_forces_netsock
  subroutine init_netsock(parini)
  use mod_parini, only: typ_parini
  implicit none
  type(typ_parini), intent(in):: parini
end subroutine init_netsock
  subroutine send_data(pos,latvec,nat,repid,msg,nmsg,latvec_rot)
  use mod_potential, only: sock_socket, sock_inet, sock_port,sock_host,MSGLEN,sock_extra_string,sock_ecutwf
  implicit none
  integer,intent(in):: nat, nmsg, repid
  real(8),intent(in):: latvec(3,3),pos(3,nat)
  real(8),intent(out)::latvec_rot(3,3)
  real(8):: latvec_inv(3,3),pos_back(3,nat),pos_cart(3,nat),dist_ang(6)
  character*1024:: msg
end subroutine send_data
       subroutine latvec2dist_ang(dist_ang,latvec,pi)
       implicit none
       real(8):: dist_ang(6),latvec(3,3),pi,convang
end subroutine latvec2dist_ang
       subroutine dist_ang2latvec(dist_ang,latvec,pi)
       implicit none
       real(8):: dist_ang(6),latvec(3,3),pi,convang
end subroutine dist_ang2latvec
  subroutine get_data(etot,fcart,strten,latvec,latvec_rot,nat)
  use mod_potential, only: sock_socket, sock_inet, sock_port,sock_host,MSGLEN,sock_extra_string,sock_ecutwf
  implicit none
  integer,intent(in) :: nat
  real(8),intent(in) :: latvec(3,3),latvec_rot(3,3)
  real(8),intent(out):: fcart(3,nat),etot,strten(6)
end subroutine get_data
       subroutine rotmat_fcart_stress(latvec_init,latvec_trans,rotmat)
       implicit none
       real(8):: latvec_init(3,3),latvec_trans(3,3),latvec_trans_inv(3,3),rotmat(3,3)
end subroutine rotmat_fcart_stress
       subroutine rotate_stresstensor(strten,rotmat)
       implicit none
       real(8):: strten(6),rotmat(3,3),stress(3,3)
end subroutine rotate_stresstensor
  subroutine final_netsock()
  implicit none
  character*1024:: host
end subroutine final_netsock
! ./src/potential_PLATO.F90 :
subroutine init_potential_forces_plato(atoms_t)
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_atoms), intent(inout):: atoms_t
end subroutine init_potential_forces_plato
subroutine cal_potential_forces_plato(iproc,n,rat,fat,epot)
    implicit none
    integer, intent(in):: iproc, n
    real(8), intent(inout):: rat(3,n/3), fat(3,n/3), epot
end subroutine cal_potential_forces_plato
subroutine final_potential_forces_plato
    implicit none
end subroutine final_potential_forces_plato
! ./src/potential_QSC.F90 :
subroutine init_potential_forces_qsc(atoms_t)
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_atoms), intent(in):: atoms_t
end subroutine init_potential_forces_qsc
subroutine cal_potential_forces_qsc(atoms)
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_atoms), intent(inout):: atoms
end subroutine cal_potential_forces_qsc
subroutine final_potential_forces_qsc
    implicit none
end subroutine final_potential_forces_qsc
! ./src/potential_sec_main.F90 :
subroutine init_potential_forces_sec(parini,atoms)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
end subroutine init_potential_forces_sec
subroutine cal_potential_forces_sec(parini,atoms)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
end subroutine cal_potential_forces_sec
subroutine final_potential_forces_sec(atoms)
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_atoms), intent(in):: atoms
end subroutine final_potential_forces_sec
! ./src/potential_SIESTA.F90 :
subroutine init_cal_potential_forces_siesta(atoms_t)
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_atoms), intent(in):: atoms_t
end subroutine init_cal_potential_forces_siesta
subroutine cal_potential_forces_siesta(atoms)
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_atoms), intent(inout):: atoms
end subroutine cal_potential_forces_siesta
subroutine final_potential_forces_siesta
    implicit none
end subroutine final_potential_forces_siesta
! ./src/potential_VASP.F90 :
subroutine init_potential_forces_vasp(atoms_t)
        use mod_atoms, only: typ_atoms
    implicit none
    type(typ_atoms), intent(inout):: atoms_t
end subroutine init_potential_forces_vasp
subroutine cal_potential_forces_vasp(atoms)
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_atoms), intent(inout):: atoms
end subroutine cal_potential_forces_vasp
subroutine final_potential_forces_vasp
    implicit none
end subroutine final_potential_forces_vasp
subroutine add_repulsive_wall(iproc,nat,rat,cellvec,fat,epot)
    implicit none
    integer, intent(in):: iproc, nat
    real(8), intent(in):: rat(3,nat), cellvec(3,3)
    real(8), intent(inout):: fat(3,nat), epot
end subroutine add_repulsive_wall
subroutine get_output_vasp_geopt(filename1,filename2,nat,latvec,xred,fcart,energy,strten,success)
    implicit none
    character(*):: filename1
    character(*):: filename2
    integer:: nat
    real(8):: fcart(3,nat),energy,strten(6),value,latvec(3,3),xred(3,nat),str_matrix(3,3),vol,a(3,3),scaling
    logical:: success
end subroutine get_output_vasp_geopt
! ./src/processors.F90 :
subroutine initprocessors
    implicit none
end subroutine initprocessors
subroutine finalizeprocessors
    implicit none
end subroutine finalizeprocessors
! ./src/saddle_1s_dimer.F90 :
subroutine dimmethimproved(parini,iproc,atoms_s,nat,ndof,rat,epot,fat,curv,uvn,paropt)
    use mod_parini, only: typ_parini
    use mod_opt, only: typ_paropt
    use mod_atoms, only: typ_atoms, typ_file_info
    implicit none
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
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_parini), intent(in):: parini
    integer, intent(in):: iproc, nat, ndof, maxitlc, nw
    type(typ_atoms), intent(inout):: atoms_s
    real(8), intent(in):: rat(3,nat), fat(3,nat), angletol
    real(8), intent(inout):: uvn(3,nat), curv0, curv
end subroutine lowestcurvature
subroutine rotatedimer(parini,iproc,atoms_s,nat,ndof,rat,uvn,fat,curv0,curv,fnrm)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    implicit none
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
    use mod_atoms, only: typ_atoms, typ_file_info
    use mod_opt, only: typ_paropt
    implicit none
    type(typ_parini), intent(in):: parini
end subroutine surface_walking
subroutine read_input(atoms_s) !,paropt)
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_atoms), intent(inout):: atoms_s
end subroutine read_input
subroutine random_move_atoms(nat,atom_motion,cellvec,rat)
    implicit none
    integer, intent(in):: nat
    logical, intent(in):: atom_motion(3,nat)
    real(8), intent(in):: cellvec(3,3)
    real(8), intent(inout):: rat(3,nat)
end subroutine random_move_atoms
subroutine find_minima(parini,iproc,atoms_s,paropt_m,paropt_m_prec,uvn,curv,epot_m0)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, typ_file_info
    use mod_opt, only: typ_paropt
    implicit none
    type(typ_parini), intent(in):: parini
    integer, intent(in):: iproc !, nat, ndof
    type(typ_atoms), intent(in):: atoms_s
    type(typ_paropt), intent(inout):: paropt_m, paropt_m_prec
    real(8), intent(in):: uvn(3,atoms_s%nat), curv, epot_m0
end subroutine find_minima
subroutine alongnegcurvature(iproc,atoms,uvn,c)
    use mod_atoms, only: typ_atoms
    implicit none
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
    implicit none
    type(typ_parini), intent(in):: parini
    integer, intent(in):: iproc, n, nr
    type(typ_atoms), intent(inout):: atoms_s
    real(8), intent(inout):: x(n), f(n), uvn(n), epot
    type(typ_paropt), intent(inout):: paropt
end subroutine optimizer_saddle
subroutine cal_potential_forces_modified(parini,iproc,atoms_s,n,x,f,epot,nr,uvn,feff,curv0,curv,fold,paropt)
    use mod_parini, only: typ_parini
    use mod_opt, only: typ_paropt
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_parini), intent(in):: parini
    integer, intent(in):: iproc, n, nr
    type(typ_atoms), intent(inout):: atoms_s
    real(8), intent(inout):: x(3,atoms_s%nat)
    real(8), intent(inout):: f(3,atoms_s%nat), epot, uvn(3,atoms_s%nat), feff(3,atoms_s%nat), curv0, curv, fold(n)
    type(typ_paropt), intent(inout):: paropt
end subroutine cal_potential_forces_modified
subroutine test_convergence_saddle(n,f,curv,paropt)
    use mod_opt, only: typ_paropt
    implicit none
    integer, intent(in):: n
    real(8), intent(in):: f(n), curv
    type(typ_paropt), intent(inout):: paropt
end subroutine test_convergence_saddle
! ./src/saddle_1s_pot.F90 :
subroutine pot_initialize(parini,atoms,paropt,paropt_m)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_opt, only: typ_paropt
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
    type(typ_paropt), intent(inout):: paropt, paropt_m
end subroutine pot_initialize
! ./src/shortrange.F90 :
subroutine shortrange_init(atoms,shortrange,linked_lists,spline)
    use mod_atoms, only: typ_atoms
    use mod_shortrange, only: typ_shortrange
    use mod_linked_lists, only: typ_linked_lists
    use mod_spline, only: typ_spline
    implicit none
    type(typ_atoms), intent(in):: atoms
    type(typ_shortrange), intent(inout):: shortrange
    type(typ_linked_lists), intent(in):: linked_lists
    type(typ_spline), intent(inout):: spline
end subroutine shortrange_init
subroutine shortrange_final(linked_lists,spline)
    use mod_linked_lists, only: typ_linked_lists
    use mod_spline, only: typ_spline
    implicit none
    type(typ_linked_lists), intent(in):: linked_lists
    type(typ_spline), intent(inout):: spline
end subroutine shortrange_final
subroutine set_interaction(atoms,shortrange)
    use mod_atoms, only: typ_atoms
    use mod_shortrange, only: typ_shortrange
    implicit none
    type(typ_atoms), intent(inout):: atoms
    type(typ_shortrange), intent(inout):: shortrange
end subroutine set_interaction
subroutine cal_shortenergy(parini,shortrange,atoms,linked_lists,spline,alpha,cell,epot_short)
    use mod_parini, only: typ_parini
    use mod_shortrange, only: typ_shortrange
    use mod_linked_lists, only: typ_linked_lists
    use mod_spline, only: typ_spline
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_shortrange), intent(in):: shortrange
    type(typ_atoms), intent(inout):: atoms
    type(typ_linked_lists), intent(inout):: linked_lists
    type(typ_spline), intent(in):: spline
    real(8), intent(in):: alpha
    real(8), intent(out):: cell(3)
    real(8), intent(out):: epot_short !short range electrostatic energy
end subroutine cal_shortenergy
! ./src/solve_poisson_cube.F90 :
subroutine solve_poisson(parini)
    use mod_parini, only: typ_parini
    use mod_electrostatics, only: typ_poisson_p3d, typ_ewald_p3d
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_parini), intent(in):: parini
end subroutine solve_poisson
! ./src/task_ann.F90 :
subroutine task_ann(parini)
    use mod_parini, only: typ_parini
    implicit none
    type(typ_parini), intent(in):: parini
end subroutine task_ann
! ./src/task_confcomp.F90 :
subroutine conf_comp(parini)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms_all
    implicit none
    type(typ_parini), intent(in):: parini
end subroutine conf_comp
subroutine set_fpall_ann(atoms_all)
    use mod_atoms, only: typ_atoms_all
    implicit none
    type(typ_atoms_all), intent(inout):: atoms_all
end subroutine set_fpall_ann
subroutine set_fpall_angle(atoms_all)
    use mod_atoms, only: typ_atoms_all
    implicit none
    type(typ_atoms_all), intent(inout):: atoms_all
end subroutine set_fpall_angle
subroutine set_fpall_distance(atoms_all)
    use mod_atoms, only: typ_atoms_all
    implicit none
    type(typ_atoms_all), intent(inout):: atoms_all
end subroutine set_fpall_distance
subroutine build_images(atoms,natpmax,natp,ratp)
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_atoms), intent(in):: atoms
    integer, intent(in):: natpmax
    integer, intent(inout):: natp
    real(8), intent(inout):: ratp(3,natpmax)
end subroutine build_images
! ./src/task_genconf.F90 :
subroutine task_genconf(parini)
    use mod_genconf, only: typ_genconf
    use mod_parini, only: typ_parini
    implicit none
    type(typ_parini), intent(inout):: parini
end subroutine task_genconf
! ./src/task_geopt.F90 :
subroutine geopt(parini)
    use mod_parini, only: typ_parini
    use mod_opt, only: typ_paropt
    use mod_atoms, only: typ_atoms, typ_file_info
    implicit none
    type(typ_parini), intent(in):: parini
end subroutine geopt
subroutine init_geopt(parini,paropt,paropt_prec)
    use mod_parini, only: typ_parini
    use mod_opt, only: typ_paropt
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_paropt), intent(inout):: paropt, paropt_prec
end subroutine init_geopt
! ./src/task_linkedlist.F90 :
subroutine  linkedlist_test(parini)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, type_pairs, typ_file_info
    use mod_linked_lists, only: typ_linked_lists, typ_pia_arr
    implicit none
    type(typ_parini), intent(in):: parini
end subroutine linkedlist_test
subroutine callinkedlist(parini,atoms,rcut,posat1st,nim,conf)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, type_pairs
    use mod_linked_lists, only: typ_linked_lists
    implicit none
    type(typ_parini):: parini
    type(typ_atoms):: atoms 
    type(type_pairs):: posat1st(atoms%nat) 
    integer:: istat, nim(atoms%nat)
    integer:: iat, jat, niat,kat,kkz,conf
    real(8):: sclinv,cell(3) ,rcut
end subroutine callinkedlist
subroutine sort_alborz(i ,j ,k,conf)
    implicit none
    integer :: i, j, k
    integer ::conf
end subroutine sort_alborz
subroutine sort2_alborz(i ,j ,k,conf,num)
    implicit none
    integer :: i, j, k
    integer ::conf,num
end subroutine sort2_alborz
subroutine genrandomconf(atoms,numb,conf)
    use mod_atoms, only: typ_atoms
    use mod_atoms, only: typ_atoms, typ_file_info
    integer ::mat,conf
    character(2):: numb
    type(typ_atoms):: atoms 
end subroutine genrandomconf
! ./src/task_miscellaneous.F90 :
subroutine miscellaneous_task(parini)
    use mod_parini, only: typ_parini
    implicit none
    type(typ_parini), intent(in):: parini
end subroutine miscellaneous_task
! ./src/task_potential.F90 :
subroutine alborz_as_potential_init(nat,sat)
    implicit none
    integer, intent(in):: nat
    character(5), intent(in):: sat(nat)
end subroutine alborz_as_potential_init
subroutine alborz_as_potential_get(boundcond,nat,cellvec,rat,sat,fat,epot,stress)
    implicit none
    character(*), intent(in):: boundcond
    integer, intent(in):: nat
    real(8), intent(in):: rat(3,nat), cellvec(3,3)
    character(5), intent(in):: sat(nat)
    real(8), intent(out):: fat(3,nat), epot, stress(3,3)
end subroutine alborz_as_potential_get
subroutine alborz_as_potential_final
    implicit none
end subroutine alborz_as_potential_final
! ./src/task_single_point.F90 :
subroutine single_point_task(parini)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms_arr, typ_file_info
    implicit none
    type(typ_parini), intent(in):: parini
end subroutine single_point_task
! ./src/task_testforces.F90 :
subroutine task_testforces(parini)
    use mod_parini, only: typ_parini
    implicit none
    type(typ_parini), intent(in):: parini
end subroutine task_testforces
subroutine testforces_fd(parini)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_parini), intent(in):: parini
end subroutine testforces_fd
subroutine teststress_fd(parini)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_parini), intent(in):: parini
end subroutine teststress_fd
subroutine teststress_fd_cellvec(parini)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_parini), intent(in):: parini
end subroutine teststress_fd_cellvec
! ./src/tightbinding.F90 :
subroutine set_indorb(partb,atoms)
    use mod_atoms, only: typ_atoms
    use mod_tightbinding, only: typ_partb
    implicit none
    type(typ_partb), intent(inout):: partb
    type(typ_atoms), intent(inout):: atoms
end subroutine set_indorb
subroutine gammaenergy(partb,atoms,natsi,pplocal)
    use mod_atoms, only: typ_atoms
    use mod_potl, only: potl_typ
    use mod_tightbinding, only: typ_partb, lenosky
    implicit none
    type(typ_partb), intent(inout):: partb
    type(typ_atoms), intent(inout):: atoms
    integer, intent(in):: natsi
    type(potl_typ), intent(inout):: pplocal
end subroutine gammaenergy
subroutine gammamat(partb,atoms,natsi,flag2,pplocal)
    use mod_tightbinding, only: typ_partb
    use mod_atoms, only: typ_atoms
    use mod_potl, only: potl_typ
    implicit none
    type(typ_partb), intent(inout):: partb
    type(typ_atoms), intent(inout):: atoms
    integer, intent(in):: natsi, flag2
    type(potl_typ), intent(in):: pplocal
end subroutine gammamat
subroutine forcediagonalizeg(partb)
    use mod_tightbinding, only: typ_partb
    implicit none
    type(typ_partb), intent(inout):: partb
end subroutine forcediagonalizeg
subroutine gammacoupling(partb,atoms,flag2,iat,jat,atomtypei,atomtypej,pplocal,rem)
    use mod_tightbinding, only: typ_partb, lenosky
    use mod_atoms, only: typ_atoms
    use mod_potl, only: potl_typ
    implicit none
    type(typ_partb), intent(inout):: partb
    type(typ_atoms), intent(inout):: atoms
    integer, intent(in):: flag2, iat, jat, atomtypei, atomtypej
    type(potl_typ), intent(in):: pplocal
    real(8), intent(out):: rem(partb%nstride,partb%nstride)
end subroutine gammacoupling
subroutine slatercoupling(u,r,hgen,dhgen,flag2,mat)
    implicit none
    real(8), intent(in):: r, u(3)
    real(8), intent(in):: hgen(4), dhgen(4)
    integer, intent(in):: flag2
    real(8), intent(out):: mat(4,4)
end subroutine slatercoupling
subroutine yfdocclocal(partb)
    use mod_tightbinding, only: typ_partb
    implicit none
    type(typ_partb), intent(inout):: partb
end subroutine yfdocclocal
subroutine Hamiltonian_der(u,flag2,mat)
    implicit none
    real(8), intent(in):: u(3)
    integer, intent(in):: flag2
    real(8), intent(out):: mat(4,4)
end subroutine hamiltonian_der
! ./src/tosifumi.F90 :
subroutine set_tosifumi(atoms,tosifumi)
    use mod_atoms, only: typ_atoms
    use mod_shortrange, only: typ_tosifumi
    implicit none
    type(typ_atoms), intent(inout):: atoms
    type(typ_tosifumi), intent(inout):: tosifumi
end subroutine set_tosifumi
subroutine coulomb_free_direct(atoms)
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_atoms), intent(inout):: atoms
end subroutine coulomb_free_direct
subroutine calenergyforces(atoms,tosifumi)
    use mod_atoms, only: typ_atoms
    use mod_shortrange, only: typ_tosifumi
    implicit none
    type(typ_atoms), intent(inout):: atoms
    type(typ_tosifumi), intent(in):: tosifumi
end subroutine calenergyforces
subroutine tosifumi_parameters(s,p)
    implicit none
    character(6), intent(out):: s(10)
    real(8), intent(out):: p(5,10)
end subroutine tosifumi_parameters
end interface
end module mod_interface
!***************************************************************************************************
