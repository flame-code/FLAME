
all: run

run:
	./density_potential.py 1.0 1.0 16.0 16 16 90
	../../flame >o1
	./prep_pot_for_comp.py

clean:
	rm -f o1 master.dat rho.cube pot_analytic.cube \
		pot_p3d.cube pot_analytic.txt pot_p3d.txt \
		total_rho.cube time.yaml flame_log.yaml fort.??

.PHONY: run clean
