#ifndef TreeExtrasF
#define TreeExtrasF

#include <stdlib.h>
#include <stdio.h>
#include "Tree.h"

typedef struct EMPIRICAL_ITEM {
    double val, eff;
} TypeEmpiricalDistributionItem;

#ifdef __cplusplus
extern "C" {
#endif

double numberRankings(TypeTree *tree);
double probTreeShape(TypeTree *tree);
double probSubTreeShape(int n, TypeTree *tree);
double logNumberRankings(TypeTree *tree);
double logNumberOrders(TypeTree *tree);
double logSubNumberOrders(int n, TypeTree *tree);
double logProbTreeShape(TypeTree *tree);
double logProbSubTreeShape(int n, TypeTree *tree);
double numberTrees(int n);
double logNumberTrees(int n);


int compareEmpiricalItem(const void* a, const void* b);
TypeEmpiricalDistributionItem **rankingNumberDistribution(int n);
double *getThresholds(double prop, TypeEmpiricalDistributionItem **emp, int n);
double *getThresholdsBis(double prop, TypeEmpiricalDistributionItem **emp, int n);
void toDistribution(TypeEmpiricalDistributionItem **emp, int n);

#ifdef __cplusplus
}
#endif



#endif
