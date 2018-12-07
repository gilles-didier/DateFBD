#ifndef DrawFossilF
#define DrawFossilF
#include <stdio.h>
#include "DrawTreeGeneric.h"
#include "FossilInt.h"


typedef struct DATA_DRAW_FOSSIL_INT {
	TypeRGB color;
	double alpha, radius;
	void (*drawLineDot)(TypeRGB, double, double, double, double, double, double, TypeParamDrawTreeGeneric*);
	TypeFossilIntFeature *fos;
} TypeDataDrawFossilInt;


#ifdef __cplusplus
extern "C" {
#endif

void drawFossilInt(int n , double x, double y, TypeInfoDrawTreeGeneric *info, void *data);

#ifdef __cplusplus
}
#endif

#endif
