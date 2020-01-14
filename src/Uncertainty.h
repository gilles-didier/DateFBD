#ifndef UncertaintyF
#define UncertaintyF

#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_vector.h>
#include "Utils.h"
#include "Tree.h"
#include "Fossil.h"
#include "Model.h"
#include "Distribution.h"

typedef struct PARAM_UNCERTAINTY {
	TypeTree *tree;
	TypeFossilFeature *fos;
	double fossil;
} TypeParamUncertainty;

#ifdef __cplusplus
extern "C" {
#endif
void testUncertainty(int k, int l, double tstartA, double tstartB, double tend, TypeModelParam *model);

double getLogDensitySpecial(double time, int n, TypeTree *tree, TypeFossilFeature *fos, TypeModelParam *param);
void fillLogDistributionSpecial(TypeDistribution *logD, int n, TypeTree *tree, TypeFossilFeature *fos, TypeModelParam *param);
double getLogLikelihoodTreeFossilDebug(TypeTree *tree, TypeFossilFeature *fos, TypeModelParam *param);
void splitTreeFossil(TypeTree *tree, TypeFossilFeature *fos, TypeTree ***treeList, int *size);
double likelihood(double timeStart, double timeEnd, int nLeaves);
double getLogLikelihoodEvent(TypeTree *tree, TypeFossilFeature *fos, TypeModelParam *param);
double getLogLikelihoodTreeFossil(TypeTree *tree, TypeFossilFeature *fos, TypeModelParam *param);
double funcFossilUncertainGSL(const gsl_vector *x, void *params);
double funcFossilUncertainThreeParGSL(const gsl_vector *x, void *params);
double getItemNumber(TypeTree *tree, TypeFossilFeature *fos);
double getItemNumberSplitted(TypeTree *tree);
TypeDistribution getDistribution(int n, double step, TypeTree *tree, TypeFossilFeature *fos, TypeModelParam *param);
TypeDistribution getDistributionUn(int n, double val, TypeTree *tree, TypeFossilFeature *fos, TypeModelParam *param);
void fillDistribution(TypeDistribution *d, int n, TypeTree *tree, TypeFossilFeature *fos, TypeModelParam *param);
void fillLogDistribution(TypeDistribution *logD, double *logCond, int n, TypeTree *tree, TypeFossilFeature *fos, TypeModelParam *param);

#ifdef __cplusplus
}
#endif


#endif
