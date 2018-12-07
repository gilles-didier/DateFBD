#ifndef SimulTreeF
#define SimulTreeF

#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "Tree.h"


/*simulate a random tree with specified birth and death rates and fossil finds on this tree with rate "fossil"*/
TypeTree *simulTree(gsl_rng *rgen, double birth, double death, double maxTime);

#endif
