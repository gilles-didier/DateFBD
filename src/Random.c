#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "Random.h"



gsl_rng *rgen;


void initRandom() {
	rgen = gsl_rng_alloc(gsl_rng_random_glibc2);
}

void endRandom() {
	gsl_rng_free(rgen);
}

/*return an uniform  random value in [0,1]*/
double getUniformStd() {
	return gsl_rng_uniform(rgen);
}
/*return an uniform  random value in [0, m]*/
double getUniformCont(double m) {
	return m*gsl_rng_uniform(rgen);
}

/*return an uniform discrete random value in {0, 1, m}*/
double getUniformDisc(int n) {
	return gsl_rng_uniform_int(rgen, n);
}
/*return an exponential distributed random value in {0, 1, m}*/
double getExponential(double l) {
    return gsl_ran_exponential(rgen, l);
}
