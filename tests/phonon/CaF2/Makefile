
all: run

run:
#	phonopy -d --tol=0.001 --amplitude=0.02 --dim="-3 3 3 3 -3 3 3 3 -3"
#	phonopy -d --tol=0.001 --dim="-3 3 3 3 -3 3 3 3 -3"
	phonopy -d --tol=0.001 --dim="4 4 4"
	./runall.sh
	python get_forces.py DIR_POSCAR-00?/posout.yaml
	phonopy -d --tol=0.001 -s -p band.conf

clean:
	rm -rf DIR_POSCAR-??? FORCE_SETS POSCAR-??? SPOSCAR band.yaml disp.yaml band.pdf phonopy_disp.yaml  phonopy.yaml 

.PHONY: run clean
