# CFLAGS = -fsanitize=undefined -std=c++11 -O3 -Iinclude/
CFLAGS = -std=c++11 -O3 -Iinclude/

all: pocketcube lightsout pancake permnet rubikscube zerowalk

clean:
	rm pocketcube lightsout pancake permnet rubikscube zerowalk

pocketcube: example/pocketcube.cc
	g++  $(CFLAGS) $< -o $@

rubikscube: example/rubikscube.cc include/gdd.hh include/permutation.hh
	g++  $(CFLAGS) $< -o $@

pancake: example/pancake.cc
	g++  $(CFLAGS) $< -o $@

lightsout: example/lightsout.cc include/gdd.hh include/permutation.hh
	g++  $(CFLAGS) $< -o $@

permnet: example/permnet.cc 
	g++  $(CFLAGS) $< -o $@

zerowalk: example/zerowalk.cc include/gdd.hh include/permutation.hh
	g++  $(CFLAGS) $< -o $@

pocketcube.gdd.txt: pocketcube
	(time ./pocketcube) > pocketcube.gdd.txt
lightsout.gdd.txt: lightsout
	(time ./lightsout 3) > lightsout.gdd.txt 2>&1
	(time ./lightsout 4) >> lightsout.gdd.txt 2>&1
	(time ./lightsout 5) >> lightsout.gdd.txt 2>&1
	(time ./lightsout 6) >> lightsout.gdd.txt 2>&1
	(time ./lightsout 7) >> lightsout.gdd.txt 2>&1
pancake.gdd.txt: pancake
	echo 11 > pancake.gdd.txt
	(time ./pancake 11) >> pancake.gdd.txt 2>&1
