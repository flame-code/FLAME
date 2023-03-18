
export FPM_FC=mpiifort
export FPM_CC=mpiicc
export FPM_CXX=mpiicpc

#FFLAGS="-O0 -C -g -traceback -DMPI -DSPGLIB -DHAVE_BPS -DHAVE_LAMMPS -I$BIGDFTROOT/install/include"
FFLAGS="-O2 -DMPI -DSPGLIB -DHAVE_BPS -DHAVE_LAMMPS -I$BIGDFTROOT/install/include"
CFLAGS="-O2"
CPPFLAGS="-O2 -DQSC_STANDALONE -DMPICH_SKIP_MPICXX -DOMPI_SKIP_MPICXX=1 -DLAMMPS_LIB_MPI -I$LAMMPSROOT/src"

SRC="src"

LIBS="-L$SPGLIBROOT/src/.libs -L$LAMMPSROOT/src -L$FUTILE/install/lib -L$MKLROOT/lib/intel64 -L$BIGDFTROOT/install/lib -cxxlib -shared-intel -qmkl=sequential"

export FPM_FFLAGS="-DHAVE_MKL -I$FUTILE/install/include $FFLAGS"
export FPM_CFLAGS="-I$SRC $CFLAGS"
export FPM_CXXFLAGS="-I$SRC $CPPFLAGS"
export FPM_LDFLAGS="$LIBS"

