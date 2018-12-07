# list of executable files to produce
DIST = dist
TREE = tree
ESTI = esti
DRAW = draw

# .o files necessary to build the executables
OBJ_DIST = Model.o TreeExtras.o Uncertainty.o Utils.o Tree.o Fossil.o FossilInt.o Distribution.o ComputeSampledDistributionThread.o
OBJ_TREE = DrawDensity.o DrawFossilInt.o DrawTreeCairo.o DrawTreePSTricks.o DrawTreeTikz.o DrawTreeGeneric.o Model.o TreeExtras.o Uncertainty.o Utils.o Tree.o Fossil.o FossilInt.o Distribution.o ComputeSampledDistributionSingleTreeThread.o
OBJ_ESTI = MinimizeNLOpt.o Model.o TreeExtras.o Uncertainty.o Utils.o Tree.o Fossil.o FossilInt.o  EstimateSampledThread.o 
OBJ_DRAW = DrawFossilInt.o DrawTreeCairo.o DrawTreePSTricks.o DrawTreeTikz.o DrawTreeGeneric.o Utils.o Tree.o Fossil.o FossilInt.o DrawFossilTree.o

########### MODIFY ONLY TO CHANGE OPTIONS ############

all = DIST TREE ESTI DRAW

# compiler and its options
CC = gcc
CFLAGS = -Wall -Wno-char-subscripts -D_POSIX_SOURCE -std=c99 -Wall -pedantic -march=native -I/usr/include/cairo -D_REENTRANT

# linker and its options
LD = $(CC)

############ LIST OF EXECUTABLE TARGETS (MODIFY ONLY TO ADD AN EXECUTABLE) ##############

all: Makefile.d $(DIST) $(TREE) $(ESTI) $(DRAW)

# build the executable
$(DIST): $(OBJ_DIST)
	$(LD) $^ -o $@ -L/usr/local/lib -L/usr/lib/ -lm -lgsl -lgslcblas -D_REENTRANT -L/usr/local/lib -L/usr/lib/ -lpthread  -O3

# build the executable
$(TREE): $(OBJ_TREE)
	$(LD) $^ -o $@ -L/usr/local/lib -L/usr/lib/ -lm -lgsl -lgslcblas -lcairo -D_REENTRANT -L/usr/local/lib -L/usr/lib/ -lpthread  -O3

# build the executable
$(ESTI): $(OBJ_ESTI)
	$(LD) $^ -o $@  -L/usr/local/lib -L/usr/lib/ -lm -lgsl -lgslcblas -lnlopt  -D_REENTRANT -L/usr/local/lib -L/usr/lib/ -lpthread  -O3

# build the executable
$(DRAW): $(OBJ_DRAW)
	$(LD) $^ -o $@ -L/usr/local/lib -L/usr/lib/ -lm -lcairo -O3

############ DO NOT MODIFY ANYTHING BELOW THIS LINE ##############

# create .o from .c
.c.o:
	$(CC) $(CFLAGS) -c $<

# remove non essential files
clean:
	$(RM) *.o *~ *.log Makefile.d

# clean everything but sources
distclean: clean
	$(RM) $(EXE)

# dependencies
Makefile.d:
	$(CC) -MM $(CFLAGS) *.c > Makefile.d

# only real files can be non phony targets
.PHONY: all clean distclean debug release
