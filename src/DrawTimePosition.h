#ifndef DrawTimePositionF
#define DrawTimePositionF
#include <stdio.h>
#include "DrawTreeGeneric.h"
#include "Distribution.h"

typedef struct TIME_POSITION {
	double start, end;
} TypeTimePosition;

typedef struct DATA_DRAW_TIME_POSITION {
	TypeRGB color;
	double alpha;
	void (*fillPolygon)(TypeRGB, double, double*, double*, int, TypeParamDrawTreeGeneric*);
	TypeTimePosition pos;
} TypeDataDrawTimePosition;


#ifdef __cplusplus
extern "C" {
#endif

void drawTimePosition(TypeInfoDrawTreeGeneric *info, void *data);

#ifdef __cplusplus
}
#endif

#endif
