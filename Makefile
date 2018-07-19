CFLAGS = -O3 -Iinclude/

all: pancake lightsout

pancake: example/pancake.cc
	g++  $(CFLAGS) $< -o $@

lightsout: example/lightsout.cc 
	g++  $(CFLAGS) $< -o $@


