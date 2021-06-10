#ifndef ExtinctionCIF
#define ExtinctionCIF

#include "Tree.h"
#include "Fossil.h"

#ifdef __cplusplus
extern "C" {
#endif

double maxIntervalStraussSadler(double risk, int n, TypeTree *tree, TypeFossilFeature *fos);
double maxLikelihoodStraussSadler(int n, TypeTree *tree, TypeFossilFeature *fos);

#ifdef __cplusplus
}
#endif


#endif
