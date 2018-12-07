#ifndef SimulFossilF
#define SimulFossilF

#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "Fossil.h"

/*adds fossil finds with rate "fossil" on the tree "tree"*/
TypeFossilFeature *addFossils(gsl_rng *rgen, double fossil, TypeTree *tree);
/*simulate a random tree with specified birth and death rates and fossil finds on this tree with rate "fossil"*/
TypeTree *simulFossilTree(gsl_rng *rgen, double birth, double death, double fossil, double maxTime);

#endif
