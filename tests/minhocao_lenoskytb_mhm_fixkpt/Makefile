
all: run

run:
	cp -f earr.dat.bak earr.dat
	cp -f ioput.bak ioput
	cp -f poscur.ascii.bak poscur.ascii
	../../flame >o1

clean:
	rm -f posgeopt.????.ascii posmd.????.ascii possoft_lat????.ascii \
		possoft_at????.ascii fcart.bin fp.bin latvec.bin poslow.bin strten.bin \
		earr.dat ioput poscur.ascii o1 poslocm_?????.ascii poslow?????.ascii \
		time.yaml posinit.ascii global.mon flame_log.yaml

.PHONY: run clean
