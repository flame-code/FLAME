
all: run

run:
	mkdir -p monminhopp
	cp -f input.minhopp.bak input.minhopp
	../../../flame >o1

clean:
	rm -f posout.acf earr.dat poslow.acf posbest.acf CPUlimit fort.?? fort.??? o1 master.dat \
		monminhopp/monitoring.??? input.minhopp \
		traj_proc???_step???????_???.acf force_traj_proc???_step???????_???.acf \
		flame_log.yaml time.yaml

.PHONY: run clean
