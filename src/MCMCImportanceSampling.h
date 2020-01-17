#ifndef MCMCImportanceSamplingF
#define MCMCImportanceSamplingF

#include <stdlib.h>
#include <stdio.h>

#include "Tree.h"
#include "Fossil.h"
#include "FossilInt.h"
#include "Model.h"
#include "MinimizeNLOpt.h"

typedef struct MCMC_PARAM {
	double al;
	int burn, gap, iter;
} TypeMCMCParam;
	
#ifdef __cplusplus
extern "C" {
#endif

TypeDistribution MCMCSamplingDist(FILE *fout, FILE *find, TypeTree **tree, int nTree, int *node, TypeFossilIntFeature **fi, double step, double al, int burn, int gap, int iter, double prop, TypeModelParam *windSize, TypeModelParam *init, double probSpe, double probExt);
void MCMCSamplingDiag(FILE *fout, FILE *find, TypeTree **tree, int nTree, TypeFossilIntFeature **fi, double al, int burn, int gap, int iter, double prop, TypeModelParam *windSize, TypeModelParam *init, double probSpe, double probExt);
TypeDistribution *MCMCSamplingMultipleDistSingleTree(int *node, int sizeNode, TypeTree *tree, TypeFossilIntFeature *fi, double step, double al, int burn, int gap, int iter, double prop, TypeModelParam *windSize, TypeModelParam *init, double probSpe, double probExt);

#ifdef __cplusplus
}
#endif



#endif
