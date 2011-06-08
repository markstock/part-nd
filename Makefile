# Makefile for part-nd, an arbitrary-dimensional particle system

# build options
#DEBUG=1
OPENMP=1
CC=gcc

# build targets
all: .part .grav
#all: part-2d part-3d
#all: grav-2d


# the neccessary c files
CFILES = main.c findvel.c ndtree.c density.c setup.c writeout.c
HFILES = structs.h

# To make a normal version 
#CFLAGS = -lm -O1 -funroll-loops
#CFLAGS = -lm -O1 -funroll-loops -pg
CFLAGS = -lm -lpng -O2 -funroll-loops -ffast-math -malign-double
#CFLAGS = -lm -pg -g -ggdb
#CFLAGS = -I/Developer/SDKs/MacOSX10.5.sdk/usr/X11/include -L/Developer/SDKs/MacOSX10.5.sdk/usr/X11/lib -lm -lpng -O2 -funroll-loops -ffast-math -malign-double

ifdef DEBUG
  CFLAGS+=-g -p -ggdb -fbounds-check
else
  CFLAGS+=-O2 -funroll-loops -ffast-math -fomit-frame-pointer
  CFLAGS+=-mtune=generic
endif
ifdef OPENMP
  CFLAGS+=-fopenmp
endif

ifeq ($(UNAME), Linux)
  LFLAGS=
endif
ifeq ($(UNAME), Darwin)
  LFLAGS=-L/Developer/SDKs/MacOSX10.5.sdk/usr/X11/lib
  CFLAGS=-I/Developer/SDKs/MacOSX10.5.sdk/usr/X11/include
endif
LFLAGS+=-lm -lgfortran -lpng

DIMS=`seq 1 4`

.part : $(CFILES) $(HFILES) Makefile
	for dim in $(DIMS); do \
	  echo "$(CC) $(CFLAGS) -DDIM=$$dim -o part-$${dim}d $(CFILES) $(LFLAGS)" ;\
	  $(CC) $(CFLAGS) -DDIM=$$dim -o part-$${dim}d $(CFILES) $(LFLAGS) ;\
	done ;\
	touch $@

.grav : $(CFILES) $(HFILES) Makefile
	for dim in $(DIMS); do \
	  echo "$(CC) $(CFLAGS) -DGRAV_ONLY -DDIM=$$dim -o grav-$${dim}d $(CFILES) $(LFLAGS)" ;\
	  $(CC) $(CFLAGS) -DGRAV_ONLY -DDIM=$$dim -o grav-$${dim}d $(CFILES) $(LFLAGS) ;\
	done ;\
	touch $@

lint:
	lint -abchp $(CFILES)

clean:
	rm -f *.o part-?d grav-?d .part .grav
