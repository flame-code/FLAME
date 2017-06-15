ifdef BPS
    OPTIONS += BPS
endif

ifdef SPGLIB
	# Path to spglib
	LIB_SPGLIB =  $(SPGLIB_ROOT)/src/.libs/libsymspg.a
	SPGLIB_EX = $(SPGLIB_ROOT)/example/
	SPGLIB_OBJ= spglib_f08.o
	OBJ+= $(SPGLIB_OBJ) spglib_int.o
	LIBS+= $(LIB_SPGLIB)
	override PRE_PROC+= -DSPGLIB
endif
