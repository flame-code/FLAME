cd src
FILES=$(find . -type f -iname '*.F90')
FILES=$(sed -e 's/\.\/potential_TINKER\.F90//' <<< $FILES)
FILES=$(sed -e 's/\.\/potential_LAMMPS\.F90//' <<< $FILES)
FILES=$(sed -e 's/\.\/binaries\.F90//' <<< $FILES)
FILES=$(sed -e 's/\.\/ascii2POSCAR\.F90//' <<< $FILES)
FILES=$(sed -e 's/\.\/convex_hull\.F90//' <<< $FILES)
FILES=$(sed -e 's/\.\/expand_poslows\.F90//' <<< $FILES)
FILES=$(sed -e 's/\.\/msock_slave_template\.F90//' <<< $FILES)
FILES=$(sed -e 's/\.\/POSCAR2ascii\.F90//' <<< $FILES)
FILES=$(sed -e 's/\.\/PWSCF_restruct\.F90//' <<< $FILES)
FILES=$(sed -e 's/\.\/recompute_kpt\.F90//' <<< $FILES)
FILES=$(sed -e 's/\.\/ternaries\.F90//' <<< $FILES)
FILES=$(sed -e 's/\.\/vasp_recompute_cell\.F90//' <<< $FILES)
FILES=$(sed -e 's/\.\/vasp_recompute_kpt\.F90//' <<< $FILES)
FILES=$(sed -e 's/\.\/vasp_recompute_kpt_odd\.F90//' <<< $FILES)
#FILES=$(find . -type f -iname '*.F90' |grep -v potential_TINKER |grep -v potential_LAMMPS.F90)
makedepf90 $FILES >dep.mk
sed -i '/^\.\/interface_mod\.o/d' dep.mk
cd ..

#gcc -MM -DQSC_STANDALONE `find . -iname '*.c' ; find . -iname '*.cpp'` >cfiles.dep
