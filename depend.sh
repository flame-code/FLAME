unfinished but the following two sections may be useful for future:

----------------------------------------------------------------------
silicon:~/tmp/qqq >cat cdepend.sh
gcc -MM -DQSC_STANDALONE `find . -iname '*.c' |grep -v PLATO ; find . -iname '*.cpp'` >cfiles.dep
silicon:~/tmp/qqq >cat depend.sh
makedepf90 -b ofiles -o client `find . -iname '*.f90' |grep -v wrappers |grep -v SIESTA |grep -v PLATO |grep -v modules |grep -v ANN |grep -v potentials/QSC/test.f90 |grep -v GenConf` >ffiles.dep

----------------------------------------------------------------------

silicon:~/tmp/qqq >cat Makefile
# FC = the compiler to use 
FC=/home/ghasemi/Programs/openmpi-1.4.4/bin/mpif90
CC=/home/ghasemi/Programs/openmpi-1.4.4/bin/mpicc

# Compiler options 
FFLAGS= -O2 

# List libraries used by the program here 
LIBS= modules/libmodules.a -L${MKLPATH} ${MKLPATH}/libmkl_solver_lp64_sequential.a \
    -Wl,--start-group -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -Wl,--end-group -lpthread

# Suffix-rules:  Begin by throwing away all old suffix- 
# rules, and then create new ones for compiling  
# *.f90-files. 
.SUFFIXES:
.SUFFIXES: .c .cpp .o

all: client libcfiles.a
   
#   libcfiles.a

#ofiles/%.o: 
#   $(FC) -c -I modules $(FFLAGS) $< -o $@

.c.o: 
    $(CC) -c $(CFLAGS) $< -o $@

clean:
    rm -f *.o *.mod client ofiles/*.o libcfiles.a


include ffiles.dep
include cfiles.dep

----------------------------------------------------------------------
