#ifndef DrawFossilF
#define DrawFossilF
#include <stdio.h>
#include "DrawTreeGeneric.h"
#include "Fossil.h"

typedef struct DATA_DRAW_FOSSIL {
	TypeRGB color;
	double alpha, radius;
	void (*drawDot)(TypeRGB, double, double, double, double, TypeParamDrawTreeGeneric*);
	TypeFossilFeature *fos;
} TypeDataDrawFossil;


#ifdef __cplusplus
extern "C" {
#endif

void drawFossil(int n , double x, double y, TypeInfoDrawTreeGeneric *info, void *data);

#ifdef __cplusplus
}
#endif

#endif
