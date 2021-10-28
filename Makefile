nbody_bh: nbody_bh.c timer.c timer.h
	gcc -O3 -o nbody_bh nbody_bh.c timer.c -lm -I.

clean:
	\rm -f *.o nbody *~ *#
