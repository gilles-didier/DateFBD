# list of executable files to produce
DSPE = dspe
TSPE = tspe
DEXT = dext
TEXT = text
QEXT = qext
CEXT = cext
DIAG = diag
ESTI = esti
DRAW = draw
SAMP = samp

# .o files necessary to build the executables
OBJ_DSPE = MyR.o Model.o MCMCImportanceSampling.o TreeExtras.o Uncertainty.o Utils.o Tree.o Fossil.o FossilInt.o Distribution.o ComputeSpeciationDistribution.o
OBJ_TSPE = MyR.o Model.o MCMCImportanceSampling.o DrawDensity.o DrawFossilInt.o DrawTreeCairo.o DrawTreePSTricks.o DrawTreeTikz.o DrawTreeGeneric.o TreeExtras.o Uncertainty.o Utils.o Tree.o Fossil.o FossilInt.o Distribution.o ComputeSpeciationDistributionSingleTree.o
OBJ_DEXT = MyR.o Model.o MCMCImportanceSampling.o TreeExtras.o Uncertainty.o Utils.o Tree.o Fossil.o FossilInt.o Distribution.o ComputeExtinctionDistribution.o
OBJ_TEXT = MyR.o Model.o MCMCImportanceSampling.o DrawTimePosition.o DrawDensity.o DrawFossilInt.o DrawTreeCairo.o DrawTreePSTricks.o DrawTreeTikz.o DrawTreeGeneric.o TreeExtras.o Uncertainty.o Utils.o Tree.o Fossil.o FossilInt.o Distribution.o ComputeExtinctionDistributionSingleTree.o
OBJ_QEXT = MyR.o Model.o MCMCImportanceSampling.o TreeExtras.o Uncertainty.o Utils.o Tree.o Fossil.o FossilInt.o Distribution.o ComputeExtinctionUpperBound.o
OBJ_CEXT = MyR.o Model.o MCMCImportanceSampling.o TreeExtras.o Uncertainty.o Utils.o Tree.o Fossil.o FossilInt.o Distribution.o ComputeExtinctionComparison.o
OBJ_DIAG = MyR.o Model.o MCMCImportanceSampling.o TreeExtras.o Uncertainty.o Utils.o Tree.o Fossil.o FossilInt.o Distribution.o ComputeDiagnosticMCMC.o
OBJ_ESTI = MyR.o MinimizeNLOpt.o Model.o TreeExtras.o Uncertainty.o Utils.o Tree.o Fossil.o FossilInt.o  EstimateRates.o 
OBJ_DRAW = MyR.o DrawFossilInt.o DrawTreeCairo.o DrawTreePSTricks.o DrawTreeTikz.o DrawTimePosition.o DrawTreeGeneric.o Utils.o Tree.o Fossil.o FossilInt.o DrawFossilTree.o
OBJ_SAMP =  MyR.o DrawFossil.o DrawTreeCairo.o DrawTreePSTricks.o DrawTreeTikz.o DrawTreeGeneric.o SimulTree.o SimulFossil.o Utils.o Tree.o Fossil.o SampleTreeFossil.o 

########### MODIFY ONLY TO CHANGE OPTIONS ############

all = DSPE TSPE DEXT TEXT QEXT CEXT DIAG ESTI DRAW SAMP

# compiler and its options
CC = gcc
CFLAGS = -Wall -Wno-char-subscripts -D_POSIX_SOURCE -DBUILD_NATIVE=OFF -std=c99 -Wall -pedantic -march=native -I/usr/include/cairo -D_REENTRANT -g
#CFLAGS = -Wall -Wno-char-subscripts -D_POSIX_SOURCE -DBUILD_NATIVE=OFF -std=c99 -Wall -pedantic -march=native -D_REENTRANT

# linker and its options
LD = $(CC)

############ LIST OF EXECUTABLE TARGETS (MODIFY ONLY TO ADD AN EXECUTABLE) ##############

all: Makefile.d $(DSPE) $(DEXT) $(QEXT) $(DIAG) $(TSPE) $(TEXT) $(ESTI) $(DRAW) $(SAMP)

# build the executable
$(DSPE): $(OBJ_DSPE)
	$(LD) $^ -o $@ -L/usr/local/lib -L/usr/lib/ -lm -lgsl -lgslcblas -D_REENTRANT -L/usr/local/lib -L/usr/lib/ -lpthread  -O3
# build the executable
$(DEXT): $(OBJ_DEXT)
	$(LD) $^ -o $@ -L/usr/local/lib -L/usr/lib/ -lm -lgsl -lgslcblas -D_REENTRANT -L/usr/local/lib -L/usr/lib/ -lpthread  -O3
# build the executable
$(QEXT): $(OBJ_QEXT)
	$(LD) $^ -o $@ -L/usr/local/lib -L/usr/lib/ -lm -lgsl -lgslcblas -D_REENTRANT -L/usr/local/lib -L/usr/lib/ -lpthread  -O3
# build the executable
$(CEXT): $(OBJ_CEXT)
	$(LD) $^ -o $@ -L/usr/local/lib -L/usr/lib/ -lm -lgsl -lgslcblas -D_REENTRANT -L/usr/local/lib -L/usr/lib/ -lpthread  -O3
# build the executable
$(DIAG): $(OBJ_DIAG)
	$(LD) $^ -o $@ -L/usr/local/lib -L/usr/lib/ -lm -lgsl -lgslcblas -D_REENTRANT -L/usr/local/lib -L/usr/lib/ -lpthread  -O3
# build the executable
$(TSPE): $(OBJ_TSPE)
	$(LD) $^ -o $@ -L/usr/local/lib -L/usr/lib/ -lm -lgsl -lgslcblas -lcairo -D_REENTRANT -L/usr/local/lib -L/usr/lib/ -lpthread  -O3
## build the executable
$(TEXT): $(OBJ_TEXT)
	$(LD) $^ -o $@ -L/usr/local/lib -L/usr/lib/ -lm -lgsl -lgslcblas  -lcairo -D_REENTRANT -L/usr/local/lib -L/usr/lib/ -lpthread  -g
## build the executable
$(ESTI): $(OBJ_ESTI)
	$(LD) $^ -o $@  -L/usr/local/lib -L/usr/lib/ -lm -lgsl -lgslcblas -lnlopt  -D_REENTRANT -L/usr/local/lib -L/usr/lib/ -lpthread  -O3
# build the executable
$(DRAW): $(OBJ_DRAW)
	$(LD) $^ -o $@ -L/usr/local/lib -L/usr/lib/ -lm -lgsl -lgslcblas -lcairo -g
# build the executable
$(SAMP): $(OBJ_SAMP)
	$(LD) $^ -o $@ -lm -lgsl -lgslcblas -lnlopt -lcairo -g

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
