# Give here the list of existing tests
SHORT_TESTDIRS = testforces_FD \
				 SP_QSC \
				 geopt_FIRE \
				 geopt_BFGS \
				 geopt_SQNM \
				 minhopp_LJ \
				 saddle_1s_Si64Bulk \
				 saddle_1s_barsad \
				 saddle_1s_barsad_dir \
				 SP_CENT1_bulk \
				 SP_CENT1_free \
				 SP_CENT1_slab \
				 SP_PP_p3d \
				 SP_PP_p3dbias \
				 geopt_CENT1_slab_fixed_efield \
				 vibfreq_fourpoints \
				 train_atombased \
				 train_CENT1 \
				 train_CENTT \
				 saddle_2s_splsad \
				 saddle_2s_neb

if HAVE_SPGLIB
SHORT_TESTDIRS += minhocao_CENT1 \
				  minhocao_LTB \
				  minhocao_LTB_fixkpt
endif

if HAVE_LAMMPS
SHORT_TESTDIRS += lammps_NVT \
				  lammps_NPT
endif

LONG_TESTDIRS = $(SHORT_TESTDIRS)

EXTRA_DIST = 

pythondir = $(FUTILE)/install/lib/python$(PYTHON_VERSION)/site-packages

include $(srcdir)/check.mk
