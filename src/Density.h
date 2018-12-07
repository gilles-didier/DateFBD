#ifndef DensityF
#define DensityF

#include <stdlib.h>
#include <stdio.h>

typedef struct DENSITY_ITEM {
    double val, dens;
} TypeDensityItem;

typedef struct DENSITY {
    int size;
    TypeDensityItem *item;
} TypeDensity;

#ifdef __cplusplus
extern "C" {
#endif


int compareDensityItem(const void* a, const void* b);
TypeDensity meanDensity(TypeDensity d1, double w1, TypeDensity d2, double w2);
void fprintDensity(FILE *f, TypeDensity d);
double sumDensity(TypeDensity d);

#ifdef __cplusplus
}
#endif



#endif
