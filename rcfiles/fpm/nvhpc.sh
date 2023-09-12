CLEAN_YESNO=yes

SPGLIB_YESNO=yes
#FLAME needs two files from spglib, one from source,
#namely the Fortran interface "spglib_f08.f90"
#and the other from installation, namely the library presumably
#called "libsymspg.a".
#In old versions of spglib, "spglib_f08.f90" was in "example" directory,
#but now it is in "fortran" directory.
#Depending on how you install spglib, these two files may share the parent
#directory or be in completely different paths.
SPGLIB_F08=full_path_to_Fortran_interface/spglib_f08.f90
LIBSYMSPG=full_path_to_library_file/libsymspg.a

#Linking to BigDFT built with NVHPC compilers:
#tested by S. A. Ghasemi, it does not work!
BigDFT_YESNO=no
#In the root directory of a BigDFT installation,
#one can find a directory called "install", in there,
#two subdirectories called "include" and "lib"
BIGDFTROOT=full_path_to_bigdft_root

#Linking to LAMMPS built with NVHPC compilers: not tested yet
LAMMPS_YESNO=no
#FLAME needs the header files and the library file from LAMMPS.
#They usually share the same parent directory, however, it may
#depend on installation.
#header files are found in the "src" directory of LAMMPS source.
#The name of library file usually depends on the compiler,
#you can search for *liblammps*.a in the LAMMPS installation directory.
LAMMPS_HEADER=full_path_to_LAMMPS_header_files
LIBLAMMPS=full_path_to_LAMMPS_library_file/liblammps.a

FC=mpif90
CC=mpicc
CXX=mpicxx

FFLAGS="-Mbackslash -fast -nomp -Minfo=all"
CFLAGS="-fast -nomp -Minfo=all"
CXXFLAGS="-fast -nomp -Minfo=all"
#FFLAGS="-Mbackslash -C -g -gopt -Mbounds -Mchkptr -Mchkstk -Mcoff -Mdwarf2 -Mdwarf3 -Melf -Mnodwarf -traceback"
#CFLAGS="-g -gopt -Mbounds -Mchkstk -Mcoff -Mdwarf2 -Mdwarf3 -Melf -Mnodwarf -traceback"
#CXXFLAGS="-g -gopt -Mbounds -Mchkstk -Mcoff -Mdwarf2 -Mdwarf3 -Melf -Mnodwarf -traceback"
PATH_OPENMPI_LIB=path_to_openmpi_lib
LDFLAGS="-c++libs -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl -L$PATH_OPENMPI_LIB -lmpi_cxx -lmpi"
MKLPATH=$MKLROOT/lib/intel64
