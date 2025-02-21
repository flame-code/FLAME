CLEAN_YESNO=yes

#Linking to spglib built with Intel LLVM compilers: not tested yet
SPGLIB_YESNO=no
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

#Linking to BigDFT built with Intel LLVM compilers: not tested yet
BigDFT_YESNO=no
#In the root directory of a BigDFT installation,
#one can find a directory called "install", in there,
#two subdirectories called "include" and "lib"
BIGDFTROOT=full_path_to_bigdft_root

#Linking to LAMMPS built with Intel LLVM compilers: not tested yet
LAMMPS_YESNO=no
#FLAME needs the header files and the library file from LAMMPS.
#They usually share the same parent directory, however, it may
#depend on installation.
#header files are found in the "src" directory of LAMMPS source.
#The name of library file usually depends on the compiler,
#e.g. in the case of Intel compiler and Intel MPI, it is
#called "liblammps_intel_cpu_intelmpi.a"
LAMMPS_HEADER=full_path_to_LAMMPS_header_files
LIBLAMMPS=full_path_to_LAMMPS_library_file/liblammps_intel_cpu_intelmpi.a

FC=mpiifx
CC=mpiicx
CXX=mpiicpx

FFLAGS=-O2
CFLAGS=-O2
CXXFLAGS=-O2
LDFLAGS="-cxxlib -shared-intel -qmkl=sequential -liomp5 -lm -ldl -lrt"
MKLPATH=$MKLROOT/lib/intel64
