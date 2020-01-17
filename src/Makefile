# list of executable files to produce
DIST = dist
DIAG = diag
TREE = tree
ESTI = esti
DRAW = draw

# .o files necessary to build the executables
OBJ_DIST = MyError.o Model.o MCMCImportanceSampling.o TreeExtras.o Uncertainty.o Utils.o Tree.o Fossil.o FossilInt.o Distribution.o ComputeDistribution.o
OBJ_DIAG = MyError.o Model.o MCMCImportanceSampling.o TreeExtras.o Uncertainty.o Utils.o Tree.o Fossil.o FossilInt.o Distribution.o ComputeDiagnosticMCMC.o
OBJ_TREE = MyError.o Model.o MCMCImportanceSampling.o DrawDensity.o DrawFossilInt.o DrawTreeCairo.o DrawTreePSTricks.o DrawTreeTikz.o DrawTreeGeneric.o TreeExtras.o Uncertainty.o Utils.o Tree.o Fossil.o FossilInt.o Distribution.o ComputeDistributionSingleTree.o
OBJ_ESTI = MyError.o MinimizeNLOpt.o Model.o TreeExtras.o Uncertainty.o Utils.o Tree.o Fossil.o FossilInt.o  EstimateRates.o 
OBJ_DRAW = DrawFossilInt.o DrawTreeCairo.o DrawTreePSTricks.o DrawTreeTikz.o DrawTreeGeneric.o Utils.o Tree.o Fossil.o FossilInt.o DrawFossilTree.o

########### MODIFY ONLY TO CHANGE OPTIONS ############

all = DIST DIAG TREE ESTI DRAW

# compiler and its options
CC = gcc
CFLAGS = -Wall -Wno-char-subscripts -D_POSIX_SOURCE -std=c99 -Wall -pedantic -march=native -I/usr/include/cairo -D_REENTRANT

# linker and its options
LD = $(CC)

############ LIST OF EXECUTABLE TARGETS (MODIFY ONLY TO ADD AN EXECUTABLE) ##############

all: Makefile.d $(DIST) $(DIAG) $(TREE) $(ESTI) $(DRAW)

# build the executable
$(DIST): $(OBJ_DIST)
	$(LD) $^ -o $@ -L/usr/local/lib -L/usr/lib/ -lm -lgsl -lgslcblas -D_REENTRANT -L/usr/local/lib -L/usr/lib/ -lpthread  -O3

# build the executable
$(DIAG): $(OBJ_DIAG)
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
