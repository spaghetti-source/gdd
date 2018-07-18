CFLAGS = -O3 -Iinclude/
pancake: example/pancake.cc
	g++  $(CFLAGS) $< -o $@
