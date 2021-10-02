noinst_LIBRARIES= libsrc.a libclowlevel.a

SUBDIRS = lenosky_tb

all: libsrc.a basic_c_wrappers.o input_variables_definition.o

bin_PROGRAMS = flame

SRC_F = minpack/dpmpar.f minpack/enorm.f  minpack/lmder.f minpack/lmpar.f minpack/qrfac.f minpack/qrsolv.f

AM_FCFLAGS = -DMPI 
if HAVE_MKL
AM_FCFLAGS += -DHAVE_MKL
endif
AM_CPPFLAGS = -DQSC_STANDALONE

SRC_ = $(SRC_F)

SRC_ += definitions_mod.F90 task_mod.F90 unitsconversion_mod.F90 processors_mod.F90 \
		constants_mod.F90 \
		spline_mod.F90 \
		splines_mod.F90 \
		greenf_kspace_mod.F90 \
		shortrange_mod.F90 \
		opt_mod.F90 genconf_mod.F90 parini_mod.F90 \
		tightbinding_mod.F90 \
		bader_mod.F90 \
		minhocao_mod.F90 \
		minhopp_mod.F90 \
		saddle_mod.F90 dynamics_mod.F90 \
		fsockets.F90 \
		timing_mod.F90 \
		constants_minhocao_mod.F90 \
		atoms_mod.F90 linked_lists_mod.F90 flame_as_potential_mod.F90 \
		electrostatics_mod.F90 \
		ann_mod.F90 \
		basic_utilities.F90 \
		train_optimizer.F90 \
		splinedsaddle_mod.F90 \
		interface_mod.F90 \
		parser_core.F90 parser_all.F90 parser_yaml.F90 \
		flame_init_fini.F90 processors.F90 \
		mpi_utilities.F90 \
		lmder_modified.F90 \
		io_xyz.F90 basic.F90 io_acf.F90 io_vasp.F90 io_cube.F90 \
		io_yaml_conf.F90 io_bin.F90 \
		ann_process.F90 ann_io.F90 ann_io_yaml.F90 ann_best_symfunc.F90 \
		ann_evaluate.F90 ann_train.F90 hung.F90 ann_check_symmetry_function.F90 ann_gen_symmetry_function.F90 task_ann.F90 \
		symfunc_data_mod.F90 \
		ann_symfunc_mod.F90 ann_symfunc_atom_behler.F90 ann_symfunc_atom_stefan.F90 \
		ann_symfunc_pair_behler.F90 \
		ann_pot_main.F90 ann_pot_atom.F90 \
		ann_pot_cent1.F90 ann_pot_cent2.F90 ann_pot_centt.F90 ann_pot_cent3.F90 ann_pot_cent_common.F90 \
		ann_pot_tb.F90 \
		shortrange.F90 tosifumi.F90 buckingham.F90 forcefield.F90 \
		minhopp_allocation.F90 minhopp.F90 minhopp_pot.F90 \
		md_util.F90 md.F90 md_NVT.F90\
		phonon.F90 \
		es_coulomb_p3d.F90 es_coulomb_spline.F90 es_hartree_p3d.F90 es_hartree_kspace.F90 specfun.F90 \
		es_coulomb_p3d_bias.F90 \
		es_coulomb_p3d_dielec.F90 \
		grid_rp4gto_sym.F90  grid_gto_sym_ortho.F90 \
		grid_gto_sym.F90 grid_basic.F90 \
		es_hartree_fourier.F90 es_hartree_main.F90 es_hartree_bps.F90 \
		cell_linkedlists.F90 \
		optimizer_sd.F90 optimizer_cg.F90 optimizer_diis.F90 optimizer_fire.F90 \
		optimizer_gmdfire.F90 optimizer_bfgs.F90 optimizer_dfp.F90 optimizer_nlbfgs.F90 \
		optimizer_sqnm.F90 \
		optimizer_drivers.F90 optimizer_drivers_vc.F90 task_geopt.F90 task_single_point.F90 \
		task_netsock.F90 task_lammps.F90 \
		saddle_1s_pot.F90 saddle_1s.F90 saddle_1s_dimer.F90 saddle_1s_optimizer.F90 \
		task_saddle.F90 \
		io_utils.F90 \
		bader_neargrid.F90 bader_ongrid.F90 bader_weight.F90 \
		task_confcomp.F90 \
		task_bader.F90 \
		genconf_random.F90 genconf_rangrow.F90 genconf_diatomic.F90 genconf_trimer.F90 task_genconf.F90 \
		task_testforces.F90 \
		task_miscellaneous.F90 solve_poisson_cube.F90 fit_elecpot.F90 get_qat_target.F90 \
		test_free_bps.F90 \
		tightbinding.F90 lenosky_tightbinding.F90 \
		task_linkedlist.F90 \
		task_potential.F90 \
		potential_main.F90 potential_sec_main.F90 potential_main_vc.F90 potential_ANN.F90 potential_LJ.F90 \
		potential_MPMD.F90 potential_BLJ_vc.F90 potential_BigDFT.F90 potential_VASP.F90 \
		potential_QSC.F90 potential_LTB.F90 potential_FF.F90 potential_NetSock.F90 \
		potential_DFTB.F90 \
		optimizer_simplex.F90 \
		splinedsaddle.F90 barsaddle.F90

if HAVE_LAMMPS
SRC_ += potential_LAMMPS_interface.F90 potential_LAMMPS.F90 lammps_mod.F90
AM_FCFLAGS += -DHAVE_LAMMPS
AM_CPPFLAGS += -I$(LAMMPS_ROOT)/src -DOMPI_SKIP_MPICXX=1 -DMPICH_SKIP_MPICXX -DLAMMPS_LIB_MPI
LIBS += $(LAMMPS_ROOT)/src/liblammps_mpi.a
endif

minhocao_src = potential_flame.F90  potential_corerepulsion.F90 \
		potential_confinement.F90 potential_IPI.F90 potential_MSOCK.F90 \
		potential_abinit.F90 potential_LJ_voids.F90 potential_BLJ_minhocao.F90 \
		potential_MLJ.F90 potential_LenoskyTB_minhocao.F90 \
		potential_LenoskyTB_LJ_minhocao.F90 potential_LenoskyMEAM.F90 \
		potential_PWSCF.F90 potential_VASP_minhocao.F90 \
		potential_SIESTA_minhocao.F90 potential_CP2K.F90 \
		potential_TERSOFF.F90 potential_EDIP.F90 potential_DFTB_minhocao.F90 \
		potential_MOPAC.F90 parser_minhocao.F90 parser_core_minhocao.F90 \
		io_vasp_minhocao.F90 io_ascii.F90 potential_main_minhocao.F90 \
		spher_harm_mathematica.F90 dynamics_md_fixlat.F90 \
		find_symmetry.F90 optimizer_sqnm_minhocao_module.F90 \
		optimizer_sqnm_minhocao.F90 \
		optimizer_bfgs_qe.F90 optimizer_sd_minhocao.F90 \
		optimizer_bfgs_minhocao.F90 optimizer_nlbfgs_minhocao.F90 \
		quaternions.F90 envelope.F90 cell_niggli.F90 \
		fingerprint_GOM.F90 fingerprint_MOLGOM.F90 fingerprint_oganov.F90 \
		fingerprint_oganov_cont.F90 fingerprint_BCM.F90 \
		fingerprint_atorb.F90 fingerprint_gaussmol.F90 \
		fingerprint_XYZ2SM.F90 atoms_minhocao.F90 \
		acceleration.F90 basic_minhocao.F90 cell_oganov.F90 \
		compare_lammps.F90 convcheck.F90 correct_latvec.F90 enthalpy.F90 \
		fingerprint.F90 fp_distance.F90 fpos_flat.F90 fragments.F90 \
		gaussdist.F90 identical.F90 inertia_tensor.F90 init_rotvels.F90 \
		init_vel.F90 insert.F90 latticetools_minhocao.F90 \
		logo_minhocao.F90 md_minhocao.F90 md_minhocao_andersen.F90 \
		md_minhocao_rbmd.F90 optimizer_fire_minhocao.F90 pbc_distance.F90 \
		propagate.F90 replace.F90 save_low_conf.F90 slab_stress.F90 soften.F90 \
		task_minhocao.F90 minhocao_enthalpyrelax.F90 minhocao_pathintegral.F90 \
		minhocao_plot_fp_grid.F90 minhocao_poslowrelax.F90 \
		minhocao_rotate_like_crazy.F90 minhocao_varvol.F90 torque_cell.F90 \
		write_restart.F90

#dependencies
#potential_main_minhocao.o: potential_abinit.o
#interface_mod.o: ekf_mod.o
include $(srcdir)/dep.mk

libsrc_a_SOURCES = $(SRC_) $(minhocao_src)

flame_SOURCES = flame.F90
flame_LDADD = libsrc.a  lenosky_tb/liblenosky.a libclowlevel.a

libclowlevel_a_SOURCES = $(_SOURCES_CPP) $(_SOURCES_C)

#libminhocao_a_SOURCES = $(minhocao_src)

_SOURCES_C = sockets.c basic_c_wrappers.c $(DATABASE_SRC)

_SOURCES_CPP = potential_QSC_main.cpp potential_QSC_wrapper.cpp 
if HAVE_LAMMPS
    _SOURCES_CPP += potential_LAMMPS_wrapper.cpp potential_LAMMPS_wrapper2.cpp
potential_LAMMPS.o : potential_LAMMPS.F90 potential_LAMMPS_interface.o
endif

if HAVE_SPGLIB
libsrc_a_SOURCES += spglib_f08.f90 spglib_int.F90
spglib_int.o: spglib_f08.o
AM_FCFLAGS += -DSPGLIB 
LIBS += $(SPGLIB_ROOT)/src/.libs/libsymspg.a
spglib_f08.f90:
	ln -s $(SPGLIB_ROOT)/example/spglib_f08.f90 spglib_f08.f90
endif

if HAVE_BPS
override LIBS := -L$(BDIR)/install/lib/ -lPSolver-1 -latlab-1 $(LIBS)
AM_FCFLAGS += -I$(BDIR)/install/include
AM_FCFLAGS += -DHAVE_BPS
endif

EXTRA_DIST = 

YAML_DATABASE = input_variables_definition.yaml

CLEANFILES = *.mod *.MOD 

include $(srcdir)/import_database.mk

