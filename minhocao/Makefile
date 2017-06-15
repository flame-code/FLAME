include ../machines/thismachine.make

LD     = mpif90 #-C -g -traceback ##-m64 -I$(MKLROOT)/include
FC     = mpif90 #-C -g -traceback  ##-m64 -I$(MKLROOT)/include

CC     = icc
LDFLAGS =  -shared-intel -mcmodel=large -I$(MKLROOT)/include -mkl=sequential
FFLAGS  =  -shared-intel -mcmodel=large # -fPIC
CFLAGS  =  -mcmodel=large
CLINK =  -cxxlib -lm
CXXFLAGS = -O0 #-fPIC
CPPFLAGS =
CXX = icpc   # replace with your C++ compiler
LDFLAGS = -L$(MKLROOT)/lib/intel64 -lmkl_gf_lp64 -lmkl_sequential -lmkl_core -lpthread -lm



OBJ_MINHOCAO = $(MODULES) $(PARSER)

#print-%  : ; @echo $* = $($*)








