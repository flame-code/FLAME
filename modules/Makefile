all: libmodules.a
.NOTPARALLEL:

SRC_ =  task_mod.F90 unitsconversion_mod.F90 processors_mod.F90 atoms_mod.F90 \
		linked_lists_mod.F90 \
		constants_mod.F90 \
		spline_mod.F90 \
		electrostatics_mod.F90 \
		ann_mod.F90 \
		shortrange_mod.F90 \
		potential_mod.F90 opt_mod.F90 genconf_mod.F90 parini_mod.F90 \
		tightbinding_mod.F90 \
		bader_mod.F90 \
		minhocao_mod.F90 \
		interface_mod.F90 minhopp_mod.F90 \
		saddle_mod.F90 dynamics_mod.F90 \
		fsockets.F90 \
		alborz_as_potential_mod.F90 \
		timing_mod.F90 \
		constants_minhocao_mod.F90

SRCDIR = .
SRC = $(patsubst %,$(SRCDIR)/%,$(SRC_))

OBJDIR = ofiles
OBJ_ = $(SRC_:.F90=.o)
OBJ = $(patsubst %,$(OBJDIR)/%,$(OBJ_))

$(OBJDIR)/%.o: $(SRCDIR)/%.F90
	$(F90) $(FFLAGS) $(PRE_PROC) $(INCLUDES) -c $< -o $@

libmodules.a: $(OBJ)
	ar -scruv libmodules.a $(OBJ)

clean:
	rm -f $(OBJ) *.mod libmodules.a
#*****************************************************************************************
