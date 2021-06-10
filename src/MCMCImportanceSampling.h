#ifndef MCMCImportanceSamplingF
#define MCMCImportanceSamplingF

#include <stdlib.h>
#include <stdio.h>

#include "Tree.h"
#include "Fossil.h"
#include "FossilInt.h"
#include "Model.h"




#ifdef __cplusplus
extern "C" {
#endif

TypeDistribution MCMCSamplingDist(FILE *fout, FILE *find, int *node, TypeTree **tree, TypeFossilIntFeature *fi, double step, double al, int burn, int gap, int iter, double prop, TypeModelParam *windSize, TypeModelParam *init, double probSpe, double probExt);
void MCMCSamplingDiag(FILE *fout, FILE *find, TypeTree **tree, TypeFossilIntFeature *fi, double al, int burn, int gap, int iter, double prop, TypeModelParam *windSize, TypeModelParam *init, double probSpe, double probExt);
TypeDistribution *MCMCSamplingMultipleDistSingleTree(int *node, int sizeNode, TypeTree *tree, TypeFossilIntFeature *fi, double step, double al, int burn, int gap, int iter, double prop, TypeModelParam *windSize, TypeModelParam *init, double probSpe, double probExt);
TypeDistribution *MCMCSamplingExtinctionMultipleDistSingleTree(int *node, int sizeNode, double maxDisplayTime, TypeTree *tree, TypeFossilIntFeature *fi, double step, double al, int burn, int gap, int iter, double prop, TypeModelParam *windSize, TypeModelParam *init, double probSpe, double probExt);
TypeDistribution *MCMCSamplingExtinctionDist(FILE *fout, FILE *find, int **taxa, TypeTree **tree, TypeFossilIntFeature *fi, double maxDisplayTime, double step, double al, int burn, int gap, int iter, double prop, TypeModelParam *windSize, TypeModelParam *init, double probSpe, double probExt);
TypeDistribution MCMCSamplingDistExtAll(FILE *fout, FILE *find, TypeTree **tree, TypeFossilIntFeature *fi, double maxDisplayTime, double step, double al, int burn, int gap, int iter, double prop, TypeModelParam *windSize, TypeModelParam *init, double probSpe, double probExt);
TypeDistribution *MCMCSamplingExtinctionMultipleDistSingleTreeFossilFixed(int *node, int sizeNode, double maxDisplayTime, TypeTree *tree, TypeFossilFeature *fp, double step, int burn, int gap, int iter, TypeModelParam *windSize, TypeModelParam *init, double probSpe, double probExt);
TypeDistribution *MCMCSamplingExtinctionMultipleDistSingleTreeFossilFixedX(FILE *fout, FILE *find, int *node, int sizeNode, double maxDisplayTime, TypeTree *tree, TypeFossilFeature *fp, double step, int burn, int gap, int iter, TypeModelParam *windSize, TypeModelParam *init, double probSpe, double probExt);
double *MCMCSamplingExtinctionMultipleDistSingleTreeFossilFixedQuantileMean(FILE *fout, FILE *find, double quant, int *node, int sizeNode, double maxDisplayTime, TypeTree *tree, TypeFossilFeature *fp, double step, int burn, int gap, int iter, TypeModelParam *windSize, TypeModelParam *init, double probSpe, double probExt);
double *MCMCSamplingExtinctionMultipleDistSingleTreeFossilFixedQuantileQuantile(FILE *fout, FILE *find, double order, int *node, int sizeNode, double maxDisplayTime, TypeTree *tree, TypeFossilFeature *fp, double step, int burn, int gap, int iter, TypeModelParam *windSize, TypeModelParam *init, double probSpe, double probExt);
double *MCMCSamplingExtinctionMultipleDistSingleTreeFossilFixedQuantileQuantileMean(int tmp, FILE *fout, FILE *find, double order, int *node, int sizeNode, TypeTree *tree, TypeFossilFeature *fp, int burn, int gap, int iter, TypeModelParam *windSize, TypeModelParam *init, double probSpe, double probExt);
double *MCMCSamplingExtinctionDistQuantMean(FILE *fout, FILE *find, double order, int **clade, TypeTree **tree, TypeFossilIntFeature *fi, double maxDisplayTime, double al, int burn, int gap, int iter, double prop, TypeModelParam *windSize, TypeModelParam *init, double probSpe, double probExt);
double MCMCSamplingExtinctionComp(FILE *fout, FILE *find, int *cladeA, int *cladeB, TypeTree **tree, TypeFossilIntFeature *fi, double al, int burn, int gap, int iter, double prop, TypeModelParam *windSize, TypeModelParam *init, double probSpe, double probExt);
#ifdef __cplusplus
}
#endif



#endif
