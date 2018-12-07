#ifndef EmpiricalDistributionF
#define EmpiricalDistributionF

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

#ifdef __cplusplus
}
#endif



#endif
