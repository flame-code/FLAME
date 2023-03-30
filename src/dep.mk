./cell_oganov.o : ./cell_oganov.F90 
./inertia_tensor.o : ./inertia_tensor.F90 
./mpi_utilities.o : ./mpi_utilities.F90 ./processors_mod.o 
./gensymcrys_random_incell_p1.o : ./gensymcrys_random_incell_p1.F90 
./atoms_minhocao.o : ./atoms_minhocao.F90 ./minhocao_mod.o ./constants_minhocao_mod.o 
./io_utils.o : ./io_utils.F90 
./es_hartree_p3d.o : ./es_hartree_p3d.F90 fftw3.f.h ./parini_mod.o ./electrostatics_mod.o 
./optimizer_sqnm_minhocao.o : ./optimizer_sqnm_minhocao.F90 ./parini_mod.o ./optimizer_sqnm_minhocao_module.o ./minhocao_mod.o ./potential_main_minhocao.o ./constants_minhocao_mod.o ./minhocao_mod.o 
./bader_weight.o : ./bader_weight.F90 ./bader_mod.o ./parini_mod.o 
./ann_fit_bf_cent2.o : ./ann_fit_bf_cent2.F90 ./processors_mod.o ./ann_pot_cent2.o ./flm_futile_mod.o ./opt_mod.o ./cell_linkedlists.o ./cell_linkedlists.o ./ann_io_yaml.o ./parini_mod.o ./atoms_mod.o ./electrostatics_mod.o ./trial_energy_mod.o 
./envelope.o : ./envelope.F90 
./init_vel.o : ./init_vel.F90 ./parini_mod.o ./constants_minhocao_mod.o 
./shortrange_mod.o : ./shortrange_mod.F90 
./optimizer_sd_minhocao.o : ./optimizer_sd_minhocao.F90 ./parini_mod.o ./optimizer_sqnm_minhocao_module.o ./potential_main_minhocao.o ./constants_minhocao_mod.o ./minhocao_mod.o ./minhocao_mod.o 
./optimizer_sqnm.o : ./optimizer_sqnm.F90 ./potential_main.o ./processors_mod.o ./opt_mod.o ./atoms_mod.o ./parini_mod.o 
./potential_MOPAC.o : ./potential_MOPAC.F90 ./parini_mod.o ./constants_minhocao_mod.o ./minhocao_mod.o 
./save_low_conf.o : ./save_low_conf.F90 
./replace.o : ./replace.F90 
./optimizer_dfp.o : ./optimizer_dfp.F90 ./opt_mod.o 
./potential_abinit.o : ./potential_abinit.F90 ./parini_mod.o ./constants_minhocao_mod.o ./minhocao_mod.o 
./latticetools_minhocao.o : ./latticetools_minhocao.F90 ./parini_mod.o 
./task_mod.o : ./task_mod.F90 
./es_radpots_cent2.o : ./es_radpots_cent2.F90 ./electrostatics_mod.o 
./train_optimizer.o : ./train_optimizer.F90 ./basic_utilities.o ./atoms_mod.o ./ann_mod.o ./ann_symfunc_mod.o ./processors_mod.o ./ann_io_yaml.o ./parini_mod.o ./ann_mod.o 
./gensymcrys_random_atoms.o : ./gensymcrys_random_atoms.F90 
./parser_core.o : ./parser_core.F90 ./task_mod.o 
./gensymcrys_trans_uniprim.o : ./gensymcrys_trans_uniprim.F90 
./es_coulomb_spline.o : ./es_coulomb_spline.F90 ./definitions_mod.o ./spline_mod.o ./shortrange_mod.o 
./potential_sec_main.o : ./potential_sec_main.F90 ./processors_mod.o ./potential_main.o ./atoms_mod.o ./parini_mod.o 
./io_vasp_minhocao.o : ./io_vasp_minhocao.F90 ./constants_minhocao_mod.o ./minhocao_mod.o ./parini_mod.o 
./md_minhocao_rbmd.o : ./md_minhocao_rbmd.F90 ./parini_mod.o ./potential_main_minhocao.o ./constants_minhocao_mod.o ./minhocao_mod.o 
./optimizer_sd.o : ./optimizer_sd.F90 ./opt_mod.o ./parini_mod.o 
./ann_pot_cent3.o : ./ann_pot_cent3.F90 ./cell_linkedlists.o ./cell_linkedlists.o ./ann_symfunc_mod.o ./ann_mod.o ./atoms_mod.o ./parini_mod.o 
./ann_symfunc_atom_behler.o : ./ann_symfunc_atom_behler.F90 ./flm_futile_mod.o ./cell_linkedlists.o ./cell_linkedlists.o ./atoms_mod.o ./symfunc_data_mod.o ./ann_mod.o 
./genconf_mod.o : ./genconf_mod.F90 
./fingerprint_MOLGOM.o : ./fingerprint_MOLGOM.F90 ./constants_minhocao_mod.o ./parini_mod.o 
./trial_energy_mod.o : ./trial_energy_mod.F90 ./flm_futile_mod.o ./processors_mod.o ./atoms_mod.o ./electrostatics_mod.o ./parini_mod.o 
./constants_minhocao_mod.o : ./constants_minhocao_mod.F90 
./potential_ANN.o : ./potential_ANN.F90 ../src/act2_cell_linkedlist.inc ../src/act1_cell_linkedlist.inc ./cell_linkedlists.o ./cell_linkedlists.o ./train_optimizer.o ./potential_corerepulsion.o ./ann_symfunc_mod.o ./ann_io_yaml.o ./potential_main.o ./atoms_mod.o ./parini_mod.o 
./task_confcomp.o : ./task_confcomp.F90 ./io_acf.o ./atoms_mod.o ./parini_mod.o 
./task_saddle.o : ./task_saddle.F90 ./parini_mod.o 
./es_hartree_main.o : ./es_hartree_main.F90 act2_cell_linkedlist.inc act1_cell_linkedlist.inc fftw3.f.h ./cell_linkedlists.o ./cell_linkedlists.o ./electrostatics_mod.o ./atoms_mod.o ./parini_mod.o 
./task_ann.o : ./task_ann.F90 ./ann_fit_bf_cent2.o ./ann_train.o ./parini_mod.o 
./ann_symfunc_pair_behler.o : ./ann_symfunc_pair_behler.F90 ./flm_futile_mod.o ./cell_linkedlists.o ./cell_linkedlists.o ./atoms_mod.o ./symfunc_data_mod.o ./ann_mod.o 
./splinedsaddle.o : ./splinedsaddle.F90 ./splinedsaddle_mod.o ./splinedsaddle_mod.o ./io_yaml_conf.o ./potential_main.o ./atoms_mod.o ./processors_mod.o ./parini_mod.o 
./acceleration.o : ./acceleration.F90 
./buckingham.o : ./buckingham.F90 ./shortrange_mod.o ./atoms_mod.o 
./grid_basic.o : ./grid_basic.F90 ./parini_mod.o ./atoms_mod.o ./electrostatics_mod.o 
./gensymcrys_random_incell_p6.o : ./gensymcrys_random_incell_p6.F90 
./basic_minhocao.o : ./basic_minhocao.F90 
./splinedsaddle_mod.o : ./splinedsaddle_mod.F90 
./electrostatics_mod.o : ./electrostatics_mod.F90 ./greenf_kspace_mod.o ./spline_mod.o ./cell_linkedlists.o 
./flame_as_potential_mod.o : ./flame_as_potential_mod.F90 ./atoms_mod.o ./task_mod.o ./parini_mod.o 
./optimizer_bfgs.o : ./optimizer_bfgs.F90 ./opt_mod.o 
./md_util.o : ./md_util.F90 ./atoms_mod.o 
./optimizer_fire.o : ./optimizer_fire.F90 ./opt_mod.o ./parini_mod.o 
./potential_TERSOFF.o : ./potential_TERSOFF.F90 ./parini_mod.o ./minhocao_mod.o ./constants_minhocao_mod.o ./minhocao_mod.o 
./ann_train.o : ./ann_train.F90 ./cell_linkedlists.o ./cell_linkedlists.o ./basic_utilities.o ./processors_mod.o ./ann_io_yaml.o ./train_optimizer.o ./ann_symfunc_mod.o ./ann_mod.o ./parini_mod.o ./atoms_mod.o 
./constants_mod.o : ./constants_mod.F90 ./definitions_mod.o 
./gensymcrys_random_lattice.o : ./gensymcrys_random_lattice.F90 
./spglib_int.o : ./spglib_int.F90 
./cell_linkedlists.o : ./cell_linkedlists.F90 ./constants_mod.o ./flm_futile_mod.o ./atoms_mod.o 
./saddle_1s.o : ./saddle_1s.F90 ./basic_utilities.o ./constants_mod.o ./processors_mod.o ./io_yaml_conf.o ./opt_mod.o ./saddle_mod.o ./potential_main.o ./atoms_mod.o ./parini_mod.o 
./es_hartree_fourier.o : ./es_hartree_fourier.F90 ./electrostatics_mod.o ./atoms_mod.o ./parini_mod.o 
./bader_ongrid.o : ./bader_ongrid.F90 ./bader_mod.o ./parini_mod.o 
./compare_lammps.o : ./compare_lammps.F90 ./parini_mod.o ./constants_minhocao_mod.o ./potential_main_minhocao.o ./minhocao_mod.o 
./ann_mod.o : ./ann_mod.F90 ./electrostatics_mod.o ./cell_linkedlists.o 
./opt_mod.o : ./opt_mod.F90 
./optimizer_diis.o : ./optimizer_diis.F90 ./opt_mod.o 
./optimizer_drivers_vc.o : ./optimizer_drivers_vc.F90 ./potential_main.o ./opt_mod.o ./io_acf.o ./atoms_mod.o ./parini_mod.o 
./potential_MLJ.o : ./potential_MLJ.F90 ./parini_mod.o ./minhocao_mod.o ./constants_minhocao_mod.o ./minhocao_mod.o 
./saddle_1s_pot.o : ./saddle_1s_pot.F90 ./potential_main.o ./processors_mod.o ./opt_mod.o ./atoms_mod.o ./parini_mod.o 
./io_acf.o : ./io_acf.F90 ./parini_mod.o ./constants_mod.o ./atoms_mod.o 
./ann_pot_centt.o : ./ann_pot_centt.F90 ./electrostatics_mod.o ./cell_linkedlists.o ./cell_linkedlists.o ./ann_symfunc_mod.o ./ann_mod.o ./atoms_mod.o ./parini_mod.o 
./potential_main.o : ./potential_main.F90 ./processors_mod.o ./atoms_mod.o ./parini_mod.o ./shortrange_mod.o ./electrostatics_mod.o ./ann_mod.o 
./grid_rp2gto_s_ortho.o : ./grid_rp2gto_s_ortho.F90 ./parini_mod.o 
./task_bader.o : ./task_bader.F90 ./parini_mod.o 
./es_hartree_bps.o : ./es_hartree_bps.F90 ./parini_mod.o ./flm_futile_mod.o ./electrostatics_mod.o ./atoms_mod.o 
./io_xyz.o : ./io_xyz.F90 
./forcefield.o : ./forcefield.F90 ./electrostatics_mod.o ./potential_main.o ./atoms_mod.o ./parini_mod.o 
./tightbinding.o : ./tightbinding.F90 ./constants_mod.o ./tightbinding_mod.o ./cell_linkedlists.o ./tightbinding_mod.o ./atoms_mod.o 
./definitions_mod.o : ./definitions_mod.F90 
./task_lammps.o : ./task_lammps.F90 ./lammps_mod.o ./io_acf.o ./potential_LAMMPS_interface.o ./potential_main.o ./io_yaml_conf.o ./atoms_mod.o ./parini_mod.o 
./ann_symfunc_mod.o : ./ann_symfunc_mod.F90 ./ann_symfunc_pair_behler.o ./atoms_mod.o ./ann_mod.o ./flm_futile_mod.o ./symfunc_data_mod.o ./cell_linkedlists.o 
./bader_mod.o : ./bader_mod.F90 
./potential_MPMD.o : ./potential_MPMD.F90 ./processors_mod.o ./constants_mod.o ./atoms_mod.o ./potential_main.o 
./propagate.o : ./propagate.F90 ./parini_mod.o 
./saddle_1s_optimizer.o : ./saddle_1s_optimizer.F90 ./saddle_mod.o ./atoms_mod.o ./potential_main.o ./opt_mod.o ./parini_mod.o 
./optimizer_drivers.o : ./optimizer_drivers.F90 ./potential_main.o ./io_bin.o ./opt_mod.o ./atoms_mod.o ./parini_mod.o 
./splines_mod.o : ./splines_mod.F90 
./es_hartree_kspace.o : ./es_hartree_kspace.F90 fftw3.f.h ./processors_mod.o ./greenf_kspace_mod.o 
./flame.o : ./flame.F90 ./flame_as_potential_mod.o ./parini_mod.o ./task_mod.o 
./ann_pot_cent_common.o : ./ann_pot_cent_common.F90 ../src/act2_cell_linkedlist.inc ../src/act1_cell_linkedlist.inc ./cell_linkedlists.o ./cell_linkedlists.o ./atoms_mod.o ./ann_symfunc_mod.o ./ann_mod.o ./parini_mod.o 
./gensymcrys_random_incell_p2.o : ./gensymcrys_random_incell_p2.F90 
./phonon.o : ./phonon.F90 ./io_yaml_conf.o ./potential_main.o ./processors_mod.o ./atoms_mod.o ./parini_mod.o 
./minhocao_enthalpyrelax.o : ./minhocao_enthalpyrelax.F90 ./parini_mod.o ./potential_main_minhocao.o ./constants_minhocao_mod.o ./minhocao_mod.o 
./genconf_trimer.o : ./genconf_trimer.F90 ./io_acf.o ./potential_main.o ./processors_mod.o ./genconf_mod.o ./atoms_mod.o ./parini_mod.o 
./minhocao_rotate_like_crazy.o : ./minhocao_rotate_like_crazy.F90 ./parini_mod.o ./potential_main_minhocao.o ./constants_minhocao_mod.o ./minhocao_mod.o 
./grid_gto_s_ortho.o : ./grid_gto_s_ortho.F90 ./parini_mod.o 
./io_bin.o : ./io_bin.F90 ./atoms_mod.o ./parini_mod.o 
./fingerprint_gaussmol.o : ./fingerprint_gaussmol.F90 
./fingerprint_GOM.o : ./fingerprint_GOM.F90 
./potential_LJ_voids.o : ./potential_LJ_voids.F90 ./parini_mod.o ./minhocao_mod.o ./constants_minhocao_mod.o ./minhocao_mod.o 
./bader_neargrid.o : ./bader_neargrid.F90 ./bader_mod.o ./parini_mod.o 
./hung.o : ./hung.F90 
./md.o : ./md.F90 ./md_util.o ./processors_mod.o ./dynamics_mod.o ./atoms_mod.o ./potential_main.o ./io_acf.o ./parini_mod.o 
./gensymcrys_random_incell_p5.o : ./gensymcrys_random_incell_p5.F90 
./potential_LenoskyTB_minhocao.o : ./potential_LenoskyTB_minhocao.F90 ./parini_mod.o ./constants_minhocao_mod.o ./minhocao_mod.o 
./potential_BigDFT.o : ./potential_BigDFT.F90 ./processors_mod.o ./atoms_mod.o ./potential_main.o 
./grid_rp4gto_sym.o : ./grid_rp4gto_sym.F90 ./parini_mod.o ./atoms_mod.o 
./optimizer_cg.o : ./optimizer_cg.F90 ./opt_mod.o 
./task_minhocao.o : ./task_minhocao.F90 ./parini_mod.o ./minhocao_mod.o ./potential_main_minhocao.o ./constants_minhocao_mod.o ./minhocao_mod.o 
./io_vasp.o : ./io_vasp.F90 
./spher_harm_mathematica.o : ./spher_harm_mathematica.F90 
./minhocao_plot_fp_grid.o : ./minhocao_plot_fp_grid.F90 ./parini_mod.o 
./spline_mod.o : ./spline_mod.F90 
./minhocao_mod.o : ./minhocao_mod.F90 
./es_coulomb_p3d_bias.o : ./es_coulomb_p3d_bias.F90 fftw3.f.h ./cell_linkedlists.o ./parini_mod.o ./potential_main.o ./atoms_mod.o ./electrostatics_mod.o 
./potential_SIESTA_minhocao.o : ./potential_SIESTA_minhocao.F90 ./parini_mod.o ./constants_minhocao_mod.o 
./ann_pot_cent2.o : ./ann_pot_cent2.F90 ./ann_io_yaml.o ./trial_energy_mod.o ./flm_futile_mod.o ./processors_mod.o ./cell_linkedlists.o ./cell_linkedlists.o ./electrostatics_mod.o ./ann_symfunc_mod.o ./ann_mod.o ./atoms_mod.o ./parini_mod.o 
./soften.o : ./soften.F90 ./parini_mod.o ./minhocao_mod.o ./potential_main_minhocao.o ./constants_minhocao_mod.o ./minhocao_mod.o 
./potential_main_minhocao.o : ./potential_main_minhocao.F90 ./parini_mod.o ./potential_LJ_voids.o ./potential_MSOCK.o ./potential_IPI.o ./potential_EDIP.o ./potential_TERSOFF.o ./potential_flame.o ./potential_LenoskyTB_LJ_minhocao.o ./potential_PWSCF.o ./potential_MLJ.o ./potential_BLJ_minhocao.o ./potential_LenoskyMEAM.o ./potential_LenoskyTB_minhocao.o ./potential_DFTB_minhocao.o ./potential_VASP_minhocao.o ./potential_SIESTA_minhocao.o ./potential_MOPAC.o ./potential_CP2K.o ./potential_abinit.o ./constants_minhocao_mod.o ./potential_corerepulsion.o ./minhocao_mod.o ./minhocao_mod.o 
./gensymcrys_sg_ops.o : ./gensymcrys_sg_ops.F90 
./greenf_kspace_mod.o : ./greenf_kspace_mod.F90 fftw3.f.h ./splines_mod.o 
./ann_io_yaml.o : ./ann_io_yaml.F90 ./io_yaml_conf.o ./io_bin.o ./atoms_mod.o ./processors_mod.o ./ann_mod.o ./parini_mod.o 
./processors_mod.o : ./processors_mod.F90 ./flm_futile_mod.o 
./parini_mod.o : ./parini_mod.F90 ./opt_mod.o 
./cell_niggli.o : ./cell_niggli.F90 
./potential_LenoskyTB_LJ_minhocao.o : ./potential_LenoskyTB_LJ_minhocao.F90 ./parini_mod.o ./minhocao_mod.o ./constants_minhocao_mod.o ./minhocao_mod.o 
./ann_pot_atom.o : ./ann_pot_atom.F90 ./cell_linkedlists.o ./cell_linkedlists.o ./ann_symfunc_mod.o ./ann_mod.o ./atoms_mod.o ./parini_mod.o 
./flm_futile_mod.o : ./flm_futile_mod.F90 
./task_miscellaneous.o : ./task_miscellaneous.F90 ./get_qat_target.o ./parini_mod.o 
./optimizer_simplex.o : ./optimizer_simplex.F90 
./optimizer_fire_minhocao.o : ./optimizer_fire_minhocao.F90 ./parini_mod.o ./minhocao_mod.o ./potential_main_minhocao.o ./minhocao_mod.o ./constants_minhocao_mod.o ./minhocao_mod.o 
./ann_best_symfunc.o : ./ann_best_symfunc.F90 ./processors_mod.o ./ann_io_yaml.o ./atoms_mod.o ./ann_symfunc_mod.o ./ann_mod.o ./parini_mod.o 
./fingerprint_oganov_cont.o : ./fingerprint_oganov_cont.F90 
./fingerprint.o : ./fingerprint.F90 ./constants_minhocao_mod.o ./minhocao_mod.o ./minhocao_mod.o ./parini_mod.o 
./es_coulomb_p3d.o : ./es_coulomb_p3d.F90 ./parini_mod.o ./atoms_mod.o ./electrostatics_mod.o 
./find_symmetry.o : ./find_symmetry.F90 ./minhocao_mod.o ./minhocao_mod.o ./parini_mod.o 
./io_ascii.o : ./io_ascii.F90 ./constants_minhocao_mod.o ./minhocao_mod.o ./parini_mod.o 
./unitsconversion_mod.o : ./unitsconversion_mod.F90 
./grid_gto_p_ortho.o : ./grid_gto_p_ortho.F90 ./parini_mod.o 
./dynamics_md_fixlat.o : ./dynamics_md_fixlat.F90 ./potential_main_minhocao.o ./constants_minhocao_mod.o ./minhocao_mod.o ./parini_mod.o 
./logo_minhocao.o : ./logo_minhocao.F90 
./fingerprint_XYZ2SM.o : ./fingerprint_XYZ2SM.F90 
./ann_evaluate.o : ./ann_evaluate.F90 ./atoms_mod.o ./ann_symfunc_mod.o ./ann_mod.o ./parini_mod.o 
./io_yaml_conf.o : ./io_yaml_conf.F90 ./io_acf.o ./constants_mod.o ./atoms_mod.o ./parini_mod.o 
./optimizer_bfgs_minhocao.o : ./optimizer_bfgs_minhocao.F90 ./minhocao_mod.o ./potential_main_minhocao.o ./constants_minhocao_mod.o ./minhocao_mod.o ./parini_mod.o ./minhocao_mod.o 
./minhopp.o : ./minhopp.F90 ./io_bin.o ./basic_utilities.o ./potential_main.o ./io_yaml_conf.o ./opt_mod.o ./atoms_mod.o ./processors_mod.o ./minhopp_mod.o ./task_mod.o ./parini_mod.o 
./potential_PWSCF.o : ./potential_PWSCF.F90 ./parini_mod.o ./constants_minhocao_mod.o ./minhocao_mod.o 
./parser_core_minhocao.o : ./parser_core_minhocao.F90 ./minhocao_mod.o 
./potential_LTB.o : ./potential_LTB.F90 ./parini_mod.o ./tightbinding_mod.o ./atoms_mod.o ./potential_main.o 
./parser_minhocao.o : ./parser_minhocao.F90 ./minhocao_mod.o ./minhocao_mod.o ./minhocao_mod.o ./minhocao_mod.o ./potential_MSOCK.o ./potential_IPI.o ./constants_minhocao_mod.o ./minhocao_mod.o ./parini_mod.o 
./gensymcrys_random_incell_main.o : ./gensymcrys_random_incell_main.F90 
./shortrange.o : ./shortrange.F90 act2_cell_linkedlist.inc act1_cell_linkedlist.inc ./cell_linkedlists.o ./parini_mod.o ./spline_mod.o ./cell_linkedlists.o ./shortrange_mod.o ./atoms_mod.o 
./potential_BLJ_vc.o : ./potential_BLJ_vc.F90 ./potential_main.o 
./minhopp_mod.o : ./minhopp_mod.F90 
./write_restart.o : ./write_restart.F90 ./parini_mod.o 
./ann_io.o : ./ann_io.F90 ./io_acf.o ./atoms_mod.o ./processors_mod.o ./ann_mod.o ./parini_mod.o 
./grid_gto_sym.o : ./grid_gto_sym.F90 ./parini_mod.o ./atoms_mod.o 
./potential_QSC.o : ./potential_QSC.F90 ./constants_mod.o ./atoms_mod.o ./potential_main.o 
./lenosky_tightbinding.o : ./lenosky_tightbinding.F90 ./constants_mod.o ./tightbinding_mod.o ./cell_linkedlists.o ./cell_linkedlists.o ./tightbinding_mod.o ./tightbinding_mod.o ./atoms_mod.o ./parini_mod.o 
./flame_init_fini.o : ./flame_init_fini.F90 ./atoms_mod.o ./flm_futile_mod.o ./parser_core.o ./parini_mod.o ./task_mod.o ./processors_mod.o 
./potential_MSOCK.o : ./potential_MSOCK.F90 ./parini_mod.o ./minhocao_mod.o ./fsockets.o ./constants_minhocao_mod.o ./minhocao_mod.o 
./optimizer_bfgs_qe.o : ./optimizer_bfgs_qe.F90 ./parini_mod.o ./minhocao_mod.o ./potential_main_minhocao.o ./constants_minhocao_mod.o ./minhocao_mod.o 
./lammps_mod.o : ./lammps_mod.F90 ./potential_LAMMPS_interface.o ./potential_main.o ./atoms_mod.o ./parini_mod.o 
./gensymcrys_random_incell_p3.o : ./gensymcrys_random_incell_p3.F90 
./ann_symfunc_atom_stefan.o : ./ann_symfunc_atom_stefan.F90 ./atoms_mod.o ./symfunc_data_mod.o ./ann_mod.o 
./optimizer_nlbfgs.o : ./optimizer_nlbfgs.F90 ./opt_mod.o 
./barsaddle.o : ./barsaddle.F90 ./basic_utilities.o ./io_yaml_conf.o ./potential_main.o ./atoms_mod.o ./processors_mod.o ./parini_mod.o 
./init_rotvels.o : ./init_rotvels.F90 ./constants_minhocao_mod.o ./minhocao_mod.o ./parini_mod.o 
./potential_main_vc.o : ./potential_main_vc.F90 ./potential_main.o ./atoms_mod.o 
./solve_poisson_cube.o : ./solve_poisson_cube.F90 ./atoms_mod.o ./electrostatics_mod.o ./parini_mod.o 
./fragments.o : ./fragments.F90 ./minhocao_mod.o ./constants_minhocao_mod.o ./parini_mod.o 
./md_minhocao.o : ./md_minhocao.F90 ./parini_mod.o ./minhocao_mod.o ./potential_main_minhocao.o ./constants_minhocao_mod.o ./minhocao_mod.o 
./lmder_modified.o : ./lmder_modified.F90 ./ann_mod.o 
./atoms_mod.o : ./atoms_mod.F90 ./constants_mod.o ./processors_mod.o 
./potential_EDIP.o : ./potential_EDIP.F90 ./parini_mod.o ./constants_minhocao_mod.o ./minhocao_mod.o 
./potential_corerepulsion.o : ./potential_corerepulsion.F90 ./parini_mod.o ./constants_minhocao_mod.o ./minhocao_mod.o 
./task_netsock.o : ./task_netsock.F90 ./minhocao_mod.o ./io_yaml_conf.o ./io_acf.o ./fsockets.o ./constants_mod.o ./processors_mod.o ./potential_main.o ./atoms_mod.o ./parini_mod.o 
./task_genconf.o : ./task_genconf.F90 ./parini_mod.o ./genconf_mod.o 
./potential_NetSock.o : ./potential_NetSock.F90 ./parini_mod.o ./atoms_mod.o ./potential_main.o ./fsockets.o 
./saddle_mod.o : ./saddle_mod.F90 
./potential_VASP_minhocao.o : ./potential_VASP_minhocao.F90 ./parini_mod.o ./constants_minhocao_mod.o ./minhocao_mod.o 
./md_minhocao_andersen.o : ./md_minhocao_andersen.F90 ./parini_mod.o ./potential_main_minhocao.o ./constants_minhocao_mod.o ./minhocao_mod.o 
./slab_stress.o : ./slab_stress.F90 
./optimizer_nlbfgs_minhocao.o : ./optimizer_nlbfgs_minhocao.F90 ./minhocao_mod.o 
./basic.o : ./basic.F90 
./minhocao_poslowrelax.o : ./minhocao_poslowrelax.F90 ./parini_mod.o ./potential_main_minhocao.o ./constants_minhocao_mod.o ./minhocao_mod.o 
./ann_pot_cent1.o : ./ann_pot_cent1.F90 ./es_hartree_main.o ./cell_linkedlists.o ./cell_linkedlists.o ./electrostatics_mod.o ./ann_symfunc_mod.o ./ann_mod.o ./atoms_mod.o ./parini_mod.o 
./ann_process.o : ./ann_process.F90 ./ann_mod.o 
./task_geopt.o : ./task_geopt.F90 ./constants_mod.o ./io_yaml_conf.o ./processors_mod.o ./io_acf.o ./potential_main.o ./atoms_mod.o ./opt_mod.o ./parini_mod.o 
./pbc_distance.o : ./pbc_distance.F90 
./ann_check_symmetry_function.o : ./ann_check_symmetry_function.F90 ./flm_futile_mod.o ./atoms_mod.o ./io_yaml_conf.o ./ann_io_yaml.o ./ann_symfunc_mod.o ./ann_mod.o ./parini_mod.o 
./tosifumi.o : ./tosifumi.F90 ./shortrange_mod.o ./atoms_mod.o 
./optimizer_gmdfire.o : ./optimizer_gmdfire.F90 ./opt_mod.o 
./potential_DFTB.o : ./potential_DFTB.F90 ./constants_mod.o ./potential_main.o ./processors_mod.o ./atoms_mod.o 
./potential_LenoskyMEAM.o : ./potential_LenoskyMEAM.F90 ./parini_mod.o ./constants_minhocao_mod.o ./minhocao_mod.o 
./genconf_rangrow.o : ./genconf_rangrow.F90 ./io_acf.o ./constants_mod.o ./genconf_mod.o ./atoms_mod.o ./parini_mod.o 
./ann_gen_symmetry_function.o : ./ann_gen_symmetry_function.F90 ./processors_mod.o ./ann_io_yaml.o ./atoms_mod.o ./ann_symfunc_mod.o ./ann_mod.o ./parini_mod.o 
./insert.o : ./insert.F90 
./identical.o : ./identical.F90 ./parini_mod.o 
./fp_distance.o : ./fp_distance.F90 ./constants_minhocao_mod.o ./minhocao_mod.o ./parini_mod.o 
./fingerprint_atorb.o : ./fingerprint_atorb.F90 
./quaternions.o : ./quaternions.F90 
./timing_mod.o : ./timing_mod.F90 
./minhopp_pot.o : ./minhopp_pot.F90 ./potential_main.o ./processors_mod.o ./opt_mod.o ./atoms_mod.o ./parini_mod.o 
./cell_utils.o : ./cell_utils.F90 
./parser_yaml.o : ./parser_yaml.F90 ./opt_mod.o ./constants_minhocao_mod.o ./parini_mod.o 
./processors.o : ./processors.F90 ./processors_mod.o 
./ann_pot_tb.o : ./ann_pot_tb.F90 ./ann_mod.o ./constants_mod.o ./cell_linkedlists.o ./train_optimizer.o ./cell_linkedlists.o ./ann_symfunc_mod.o ./ann_mod.o ./atoms_mod.o ./tightbinding_mod.o ./tightbinding_mod.o ./parini_mod.o 
./potential_FF.o : ./potential_FF.F90 ./atoms_mod.o ./parini_mod.o 
./test_free_bps.o : ./test_free_bps.F90 ./io_yaml_conf.o ./atoms_mod.o ./electrostatics_mod.o ./parini_mod.o 
./io_cube.o : ./io_cube.F90 ./processors_mod.o ./electrostatics_mod.o ./atoms_mod.o 
./potential_DFTB_minhocao.o : ./potential_DFTB_minhocao.F90 ./parini_mod.o ./constants_minhocao_mod.o ./minhocao_mod.o 
./gensymcrys_subs_cell.o : ./gensymcrys_subs_cell.F90 
./convcheck.o : ./convcheck.F90 ./constants_minhocao_mod.o ./parini_mod.o 
./enthalpy.o : ./enthalpy.F90 
./torque_cell.o : ./torque_cell.F90 
./task_single_point.o : ./task_single_point.F90 ./minhocao_mod.o ./constants_mod.o ./processors_mod.o ./io_yaml_conf.o ./io_acf.o ./potential_main.o ./atoms_mod.o ./parini_mod.o 
./potential_IPI.o : ./potential_IPI.F90 ./parini_mod.o ./minhocao_mod.o ./fsockets.o ./constants_minhocao_mod.o ./minhocao_mod.o 
./saddle_1s_dimer.o : ./saddle_1s_dimer.F90 ./basic_utilities.o ./potential_main.o ./io_yaml_conf.o ./atoms_mod.o ./opt_mod.o ./saddle_mod.o ./parini_mod.o 
./get_qat_target.o : ./get_qat_target.F90 ./es_radpots_cent2.o ./trial_energy_mod.o ./ann_io_yaml.o ./ann_mod.o ./electrostatics_mod.o ./atoms_mod.o ./parini_mod.o 
./gensymcrys_input_creator.o : ./gensymcrys_input_creator.F90 ./constants_mod.o ./io_yaml_conf.o ./atoms_mod.o ./parini_mod.o 
./gensymcrys_random_incell_p4.o : ./gensymcrys_random_incell_p4.F90 
./genconf_random.o : ./genconf_random.F90 ./io_acf.o ./genconf_mod.o ./atoms_mod.o ./parini_mod.o 
./gaussdist.o : ./gaussdist.F90 ./basic_utilities.o 
./minhocao_varvol.o : ./minhocao_varvol.F90 ./parini_mod.o ./potential_main_minhocao.o ./constants_minhocao_mod.o ./minhocao_mod.o 
./basic_utilities.o : ./basic_utilities.F90 ./parini_mod.o ./processors_mod.o ./task_mod.o ./definitions_mod.o 
./task_testforces.o : ./task_testforces.F90 ./io_yaml_conf.o ./io_acf.o ./constants_mod.o ./processors_mod.o ./potential_main.o ./atoms_mod.o ./parini_mod.o 
./potential_flame.o : ./potential_flame.F90 ./parini_mod.o 
./potential_VASP.o : ./potential_VASP.F90 ./processors_mod.o ./atoms_mod.o ./potential_main.o 
./specfun.o : ./specfun.F90 
./fingerprint_BCM.o : ./fingerprint_BCM.F90 ./parini_mod.o 
./md_NVT.o : ./md_NVT.F90 ./processors_mod.o ./dynamics_mod.o ./md_util.o ./io_acf.o ./atoms_mod.o ./potential_main.o ./parini_mod.o 
./ann_pot_main.o : ./ann_pot_main.F90 ./ann_pot_cent2.o ./ann_symfunc_mod.o ./tightbinding_mod.o ./ann_train.o ./atoms_mod.o ./train_optimizer.o ./ann_mod.o ./parini_mod.o 
./genconf_diatomic.o : ./genconf_diatomic.F90 ./io_acf.o ./constants_mod.o ./potential_main.o ./processors_mod.o ./genconf_mod.o ./atoms_mod.o ./parini_mod.o 
./fingerprint_oganov.o : ./fingerprint_oganov.F90 ./parini_mod.o 
./symfunc_data_mod.o : ./symfunc_data_mod.F90 ./cell_linkedlists.o 
./potential_CP2K.o : ./potential_CP2K.F90 ./parini_mod.o ./constants_minhocao_mod.o ./minhocao_mod.o 
./minhopp_allocation.o : ./minhopp_allocation.F90 ./minhopp_mod.o 
./potential_BLJ_minhocao.o : ./potential_BLJ_minhocao.F90 ./parini_mod.o ./minhocao_mod.o ./constants_minhocao_mod.o ./minhocao_mod.o 
./optimizer_sqnm_minhocao_module.o : ./optimizer_sqnm_minhocao_module.F90 
./parser_all.o : ./parser_all.F90 ./genconf_mod.o ./dynamics_mod.o ./saddle_mod.o ./opt_mod.o ./parser_core.o ./parini_mod.o ./task_mod.o 
./dynamics_mod.o : ./dynamics_mod.F90 
./es_coulomb_p3d_dielec.o : ./es_coulomb_p3d_dielec.F90 fftw3.f.h ./parini_mod.o ./potential_main.o ./atoms_mod.o ./electrostatics_mod.o 
./fpos_flat.o : ./fpos_flat.F90 ./parini_mod.o 
./potential_confinement.o : ./potential_confinement.F90 ./constants_minhocao_mod.o ./minhocao_mod.o ./parini_mod.o 
./fit_elecpot.o : ./fit_elecpot.F90 ./io_yaml_conf.o ./ann_mod.o ./atoms_mod.o ./electrostatics_mod.o ./parini_mod.o 
./task_linkedlist.o : ./task_linkedlist.F90 act2_cell_linkedlist.inc act1_cell_linkedlist.inc ./io_acf.o ./cell_linkedlists.o ./cell_linkedlists.o ./constants_mod.o ./atoms_mod.o ./parini_mod.o 
./task_potential.o : ./task_potential.F90 ./atoms_mod.o ./potential_main.o ./flame_as_potential_mod.o 
./minhocao_pathintegral.o : ./minhocao_pathintegral.F90 ./parini_mod.o ./potential_main_minhocao.o ./constants_minhocao_mod.o ./minhocao_mod.o 
./tightbinding_mod.o : ./tightbinding_mod.F90 
./correct_latvec.o : ./correct_latvec.F90 
./potential_LAMMPS_interface.o : ./potential_LAMMPS_interface.F90 
./potential_LJ.o : ./potential_LJ.F90 ./atoms_mod.o 
./fsockets.o : ./fsockets.F90 
