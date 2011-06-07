# Makefile for part-nd, an arbitrary-dimensional particle system

# the neccessary c files
CFILES = main.c findvel.c ndtree.c density.c setup.c writeout.c
HFILES = structs.h

# To make a normal version 
#CFLAGS = -lm -O1 -funroll-loops
#CFLAGS = -lm -O1 -funroll-loops -pg
CFLAGS = -lm -lpng -O2 -funroll-loops -ffast-math -malign-double
#CFLAGS = -lm -pg -g -ggdb
#CFLAGS = -I/Developer/SDKs/MacOSX10.5.sdk/usr/X11/include -L/Developer/SDKs/MacOSX10.5.sdk/usr/X11/lib -lm -lpng -O2 -funroll-loops -ffast-math -malign-double


#all: part-2d part-3d
all: grav-2d

part-1d: $(CFILES) $(HFILES) Makefile
	gcc -o part-1d -DONE $(CFLAGS) $(CFILES)
	@echo "part-1d made"

part-2d: $(CFILES) $(HFILES) Makefile
	gcc -o part-2d -DTWO $(CFLAGS) $(CFILES)
	@echo "part-2d made"

part-3d: $(CFILES) $(HFILES) Makefile
	gcc -o part-3d -DTHREE $(CFLAGS) $(CFILES)
	@echo "part-3d made"

part-4d: $(CFILES) $(HFILES) Makefile
	gcc -o part-4d -DFOUR $(CFLAGS) $(CFILES)
	@echo "part-4d made"

grav-1d: $(CFILES) $(HFILES) Makefile
	gcc -o grav-1d -DONE -DGRAV_ONLY $(CFLAGS) $(CFILES)
	@echo "grav-1d made"

grav-2d: $(CFILES) $(HFILES) Makefile
	gcc -o grav-2d -DTWO -DGRAV_ONLY $(CFLAGS) $(CFILES)
	@echo "grav-2d made"

grav-3d: $(CFILES) $(HFILES) Makefile
	gcc -o grav-3d -DTHREE -DGRAV_ONLY $(CFLAGS) $(CFILES)
	@echo "grav-3d made"

grav-4d: $(CFILES) $(HFILES) Makefile
	gcc -o grav-4d -DFOUR -DGRAV_ONLY $(CFLAGS) $(CFILES)
	@echo "grav-4d made"

lint:
	lint -abchp $(CFILES)
