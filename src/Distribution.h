#ifndef DistributionF
#define DistributionF

#include <stdlib.h>
#include <stdio.h>

typedef struct DENSITY_ITEM {
    double val, dens;
} TypeDistributionItem;

typedef struct DENSITY {
    int size;
    TypeDistributionItem *item;
} TypeDistribution;

#ifdef __cplusplus
extern "C" {
#endif


TypeDistribution readDistribution(FILE *f);
int compareDistributionItem(const void* a, const void* b);
//TypeDistribution meanDistribution(TypeDistribution d1, double w1, TypeDistribution d2, double w2);
void fprintDistribution(FILE *f, TypeDistribution d);
double sumDistribution(TypeDistribution d);
double sumTrapezeDistribution(TypeDistribution d);
double sumExpTrapezeDistribution(TypeDistribution d);
double sumSimpsonDistribution(TypeDistribution d);
double sumExpSimpsonDistribution(TypeDistribution d);
void fprintDerive(FILE *f, TypeDistribution d);
double getMean(TypeDistribution d);
double getMode(TypeDistribution d);
double getMedian(TypeDistribution d);
double getMeanDens(TypeDistribution d);
double getMedianDens(TypeDistribution d);
double getQuantileInf(TypeDistribution d, double q);
double getQuantileSup(TypeDistribution d, double q);
void deriveDistribution(TypeDistribution *d);
void integrateDistribution(TypeDistribution *d);
TypeDistribution agregDistribution(TypeDistribution *d, int size, int def);
TypeDistribution resampleDistribution(TypeDistribution d, double min, double max, int def);
double getEmpiricalQuantile(TypeDistribution d, double q);

#ifdef __cplusplus
}
#endif



#endif
