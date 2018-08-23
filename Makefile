# CFLAGS = -fsanitize=undefined -std=c++11 -O3 -Iinclude/
CFLAGS = -std=c++11 -O2 -Iinclude/

all: pocketcube lightsout pancake permnet rubikscube zerowalk

clean:
	rm lightsout pancake pocketcube

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
