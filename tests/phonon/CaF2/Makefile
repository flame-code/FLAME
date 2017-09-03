
all: run

run:
#	phonopy -d --tol=0.001 --amplitude=0.02 --dim="-3 3 3 3 -3 3 3 3 -3"
	phonopy -d --tol=0.001 --dim="-3 3 3 3 -3 3 3 3 -3"
	./runall.sh
	python get_forces.py DIR_POSCAR-00?/force_posout.acf
	phonopy -d --tol=0.001 -p band.conf

clean:
	rm -rf DIR_POSCAR-??? FORCE_SETS POSCAR-??? SPOSCAR band.yaml disp.yaml band.pdf

.PHONY: run clean
