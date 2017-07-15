include ../machines/thismachine.make
include ../def_pot.mk
LIB_MOD= ../modules/libmodules.a
include ../def_opt.mk

all: libsrc.a ofiles/basic_c_wrappers.o 

SRC_F= dpmpar.f enorm.f lmder.f lmpar.f qrfac.f qrsolv.f
SRC_= alborz_init_final.F90 processors.F90 \
		mpi_utilities.F90 \
		lmder_modified.F90 \
		ann_basic.F90 ann_process.F90 ann_weights_init.F90 ann_io.F90 ann_best_symfunc.F90 \
		ann_evaluate.F90 ann_train.F90 hung.F90 ann_check_symmetry_function.F90 ann_gen_symmetry_function.F90 task_ann.F90 \
		ann_symfunc_main.F90 ann_symfunc_atom_behler.F90 ann_symfunc_atom_stefan.F90 \
		ann_symfunc_pair_behler.F90 \
		ann_ekf_rivals.F90 ann_ekf_behler.F90 ann_lm.F90 \
		ann_pot_main.F90 ann_pot_atom.F90 ann_pot_cent1.F90 ann_pot_cent2.F90 ann_pot_cent_common.F90 \
		ann_pot_tb.F90 \
		shortrange.F90 tosifumi.F90 forcefield.F90 \
		minhopp_allocation.F90 minhopp.F90 minhopp_pot.F90 \
		parser_core.F90 parser_all.F90 \
		md.F90 md_NVT.F90\
		phonon.F90 \
		es_coulomb_p3d.F90 es_coulomb_spline.F90 es_hartee_p3d.F90 \
		es_coulomb_p3d_bias.F90 \
		es_hartee_fourier.F90 es_hartee_main.F90 es_hartee_bps.F90 \
		cell_linkedlists.F90 \
		optimizer_sd.F90 optimizer_cg.F90 optimizer_diis.F90 optimizer_fire.F90 \
		optimizer_gmdfire.F90 optimizer_bfgs.F90 optimizer_dfp.F90 optimizer_nlbfgs.F90 \
		optimizer_sqnm.F90 \
		optimizer_drivers.F90 optimizer_drivers_vc.F90 task_geopt.F90 task_single_point.F90 \
		task_netsock.F90 \
		saddle_1s_pot.F90 saddle_1s.F90 saddle_1s_dimer.F90 saddle_1s_optimizer.F90 \
		basic_atoms.F90 basic_utilities.F90 \
		io_xyz.F90 basic.F90 io_acf.F90 io_vasp.F90 io_cube.F90 \
		bader_neargrid.F90 bader_ongrid.F90 bader_weight.F90 \
		task_confcomp.F90 \
		task_bader.F90 \
		genconf_random.F90 genconf_rangrow.F90 genconf_diatomic.F90 genconf_trimer.F90 task_genconf.F90 \
		task_testforces.F90 \
		task_miscellaneous.F90 solve_poisson_cube.F90 best_charge_density.F90 \
		best_charge_density.F90 \
		tightbinding.F90 lenosky_tightbinding.F90 \
		task_linkedlist.F90 \
		task_potential.F90 \
		potential_main.F90 potential_sec_main.F90 potential_main_vc.F90 potential_ANN.F90 potential_LJ.F90 \
		potential_MPMD.F90 potential_BLJ_vc.F90 potential_BigDFT.F90 potential_VASP.F90 \
		potential_QSC.F90 potential_LTB.F90 potential_FF.F90 potential_NetSock.F90 \
		potential_DFTB.F90\
		potential_alborz.F90  potential_corerepulsion.F90\
		potential_confinement.F90 potential_IPI.F90 potential_MSOCK.F90\
		potential_abinit.o potential_LJ_voids.F90 potential_BLJ_minhocao.F90\
		potential_MLJ.F90 potential_LenoskyTB_minhocao.F90\
		potential_LenoskyTB_LJ_minhocao.F90 potential_LenoskyMEAM.F90\
		potential_PWSCF.F90 potential_VASP_minhocao.F90\
		potential_SIESTA_minhocao.F90 potential_CP2K.F90\
		potential_TERSOFF.F90 potential_EDIP.F90 potential_DFTB_minhocao.F90\
		potential_MOPAC.F90 parser_minhocao.F90 parser_core_minhocao.F90\
		io_vasp_minhocao.F90 io_ascii.F90 potential_main_minhocao.F90\
		spher_harm_mathematica.F90 dynamics_md_fixlat.F90\
		find_symmetry.F90 optimizer_sqnm_minhocao_module.F90\
		optimizer_sqnm_minhocao.F90\
		optimizer_bfgs_qe.F90 optimizer_sd_minhocao.F90\
		optimizer_bfgs_minhocao.F90 optimizer_nlbfgs_minhocao.F90\
		quaternions.F90 envelope.F90 cell_niggli.F90\
		fingerprint_GOM.F90 fingerprint_MOLGOM.F90 fingerprint_oganov.F90\
		fingerprint_oganov_cont.F90 fingerprint_BCM.F90\
		fingerprint_atorb.F90 fingerprint_gaussmol.F90\
		fingerprint_XYZ2SM.F90 atoms_minhocao.F90 minhocao.F90\
		 vasp_recompute_kpt.F90 \
		expand_poslows.F90\
		convex_hull.F90 envelope.F90\
		vasp_recompute_kpt_odd.F90 \
		vasp_recompute_cell.F90 \
		binaries.F90 \
		ascii2POSCAR.F90\
		POSCAR2ascii.F90 \
		recompute_kpt.F90 \
		PWSCF_restruct.F90 \
		ternaries.F90


ifdef SPGLIB
    SRC_ += spglib_int.F90
endif

ifdef LAMMPS
    SRC_ += potential_LAMMPS_interface.F90 potential_LAMMPS.F90
endif
#LAMMPS-wrapper.o : LAMMPS-wrapper.cpp LAMMPS-wrapper.h
#	$(CXX) $(CPPFLAGS) $(CXXFLAGS)  -c $< -I$(LAMMPS_SRC)    $(MPI_STUBS)




OBJ_MINHOCAO_ = potential_alborz.o potential_corerepulsion.o \
	potential_confinement.o potential_IPI.o potential_MSOCK.o \
	potential_abinit.o potential_LJ_voids.o potential_BLJ_minhocao.o \
	potential_MLJ.o potential_LenoskyTB_minhocao.o \
	potential_LenoskyTB_LJ_minhocao.o potential_LenoskyMEAM.o \
	potential_PWSCF.o potential_VASP_minhocao.o \
	potential_SIESTA_minhocao.o potential_CP2K.o \
	potential_TERSOFF.o potential_EDIP.o potential_DFTB_minhocao.o \
	potential_MOPAC.o parser_minhocao.o \
	io_vasp_minhocao.o io_ascii.o potential_main_minhocao.o \
	spher_harm_mathematica.o dynamics_md_fixlat.o \
	find_symmetry.o optimizer_sqnm_minhocao_module.o \
	optimizer_sqnm_minhocao.o \
	optimizer_bfgs_qe.o optimizer_sd_minhocao.o \
	optimizer_bfgs_minhocao.o optimizer_nlbfgs_minhocao.o \
	quaternions.o envelope.o cell_niggli.o \
	fingerprint_GOM.o fingerprint_MOLGOM.o fingerprint_oganov.o \
	fingerprint_oganov_cont.o fingerprint_BCM.o \
	fingerprint_atorb.o fingerprint_gaussmol.o \
	fingerprint_XYZ2SM.o atoms_minhocao.o minhocao.o 

OBJ_MINHOCAO =$(patsubst %,$(OBJDIR)/%,$(OBJ_MINHOCAO_))


#ifeq ($(findstring BPS,$(OPTIONS)),BPS)
#   The following file is renamed to es_hartee_bps.F90
#	SRC_ += es_coulomb_bps.F90
#endif

_SOURCES_C = sockets.c lenosky_tb/clsfunc.c lenosky_tb/nrutil.c\
	lenosky_tb/tblink.c lenosky_tb/minimize.c lenosky_tb/tbnorth.c



lenosky_tn%.o: %.c
	$(CC) $(CFLAGS) -c $<
sockets.o: sockets.c
	$(CC) $(FLAGS) -c -o sockets.o sockets.c
%.o: %.F90
	$(FC) $(FFLAGS) $(PRE_PROC) -c $<

_SOURCES_CPP = potential_QSC_main.cpp potential_QSC_wrapper.cpp 
ifdef LAMMPS
    _SOURCES_CPP += potential_LAMMPS_wrapper.cpp
endif

ifeq ($(findstring SIESTA,$(POTENTIALS)),SIESTA)
	_SOURCES += potential_SIESTA.F90
	INCLUDES += -I../../../packages/siesta/siesta-3.1/Obj-serial
endif

ifeq ($(findstring PLATO,$(POTENTIALS)),PLATO)
	_SOURCES += potential_PLATO.F90
	_SOURCES_C += potential_PLATO_wrapper.c
	INCLUDES_C += -I/home/ghasemi/Alborz/packages/plato/plato_Horsfield/src/TB2
endif


SRCDIR = .
SRC = $(patsubst %,$(SRCDIR)/%,$(SRC_))

_OBJECTS_C = $(_SOURCES_C:.c=.o)
_OBJECTS_CPP = $(_SOURCES_CPP:.cpp=.o)

OBJDIR = ofiles
OBJ_F_ = $(SRC_F:.f=.o)
OBJ_ = $(SRC_:.F90=.o)
OBJ_F = $(patsubst %,$(OBJDIR)/%,$(OBJ_F_))
OBJ = $(patsubst %,$(OBJDIR)/%,$(OBJ_))
OBJ_ALL = $(OBJ_F) $(OBJ) $(patsubst %,$(OBJDIR)/%,$(_OBJECTS_C)) \
		  $(patsubst %,$(OBJDIR)/%,$(_OBJECTS_CPP))

ifdef SPGLIB
	# Path to spglib
	LIB_SPGLIB =  $(SPGLIB_ROOT)/src/.libs/libsymspg.a
	SPGLIB_EX = $(SPGLIB_ROOT)/example/
	SPGLIB_OBJ= ofiles/spglib_f08.o
	LIBS+= $(LIB_SPGLIB)
	override PRE_PROC+= -DSPGLIB
endif

print_vars:
	@echo $(LIB_POT) $(LIB_MOD)

$(OBJDIR)/%.o: minpack/%.f
	$(F90) $(FFLAGS) -c $< -o $@ 

$(OBJDIR)/basic_c_wrappers.o: basic_c_wrappers.c
	$(CC) -c basic_c_wrappers.c -o $(OBJDIR)/basic_c_wrappers.o

$(OBJDIR)/alborz.o: $(SRCDIR)/alborz.F90
	$(F90) $(FFLAGS) $(PRE_PROC) $(INCLUDES) -c $< -o $@

$(OBJDIR)/%.o: $(SRCDIR)/%.F90 $(OBJDIR)/alborz.o $(LIB_MOD) $(SPGLIB_OBJ)
	$(F90) $(FFLAGS) $(PRE_PROC) $(INCLUDES) -I../build/install/include -c $< -o $@

$(OBJDIR)/%.o: $(SRCDIR)/%.c
	$(CC) $(FFLAGS) -c $(INCLUDES_C) $< -o $@

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp potential_QSC_parameters.h potential_QSC_class.h
	$(CC) $(FFLAGS) -DQSC_STANDALONE -c $< -o $@

#print-%  : ; @echo $* = $($*)

$(SPGLIB_OBJ): $(SPGLIB_EX)/spglib_f08.f90
	$(F90) $(FFLAGS) -c $< -o $@
spglib_f08.mod: $(SPGLIB_EX)/spglib_f08.f90 $(SPGLIB_OBJ)
	@true

MODULE: $(MODULES) $(SPGLIB_OBJ) spglib_f08.mod 
	@true


libsrc.a: $(OBJ_ALL)
	ar -scruv libsrc.a $(OBJ_ALL)

libminhocao.a: $(OBJ_MINHOCAO) $(_OBJECTS_C)
	ar -scru libminhocao.a $(OBJ_MINHOCAO)


clean:
	rm -f $(OBJ_ALL) libsrc.a $(OBJDIR)/basic_c_wrappers.o $(OBJDIR)/alborz.o ofiles/spglib_int.o ofiles/spglib_f08.o
	rm -rf libminhocao.a vasp_recompute_kpt.o *.mod *~

#*****************************************************************************************
