cd src
makedepf90 `find . -type f -iname '*.F90' |grep -v potential_TINKER |grep -v potential_LAMMPS.F90` >dep.mk
sed -i '/^\.\/interface_mod\.o/d' dep.mk
cd ..

#gcc -MM -DQSC_STANDALONE `find . -iname '*.c' ; find . -iname '*.cpp'` >cfiles.dep
