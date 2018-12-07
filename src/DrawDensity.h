#ifndef DrawDensityF
#define DrawDensityF
#include <stdio.h>
#include "DrawTreeGeneric.h"
#include "Distribution.h"

typedef struct DATA_DRAW_DENSITY {
	TypeRGB color;
	double alpha;
	void (*fillPolygon)(TypeRGB, double, double*, double*, int, TypeParamDrawTreeGeneric*);
	TypeDistribution *dens;
} TypeDataDrawDensity;


#ifdef __cplusplus
extern "C" {
#endif

void drawDensity(int n , double x, double y, TypeInfoDrawTreeGeneric *info, void *data);

#ifdef __cplusplus
}
#endif

#endif
