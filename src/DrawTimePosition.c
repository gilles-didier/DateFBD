#include "DrawTimePosition.h"

#define PER 0.05
#define MIN 2

void drawTimePosition(TypeInfoDrawTreeGeneric *info, void *data) {
	if(data) {
		double xtab[4], ytab[4];
		ytab[0] = info->param.leafSep*(info->param.nleaves+1)+info->param.labelSep;
		ytab[1] = 0.;
		xtab[0] = (((TypeDataDrawTimePosition*)data)->pos.start-info->param.tmin)*info->param.scale+info->param.xoffset;
		xtab[2] = (((TypeDataDrawTimePosition*)data)->pos.end-info->param.tmin)*info->param.scale+info->param.xoffset;
		if(xtab[2]-xtab[0]<MIN) {
			info->funct.drawLineColorAlpha(((TypeDataDrawTimePosition*)data)->color, ((TypeDataDrawTimePosition*)data)->alpha, (xtab[2]+xtab[0])/2., ytab[0], (xtab[2]+xtab[0])/2., ytab[1], &(info->param));
		} else {
			xtab[1] = xtab[0];
			xtab[3] = xtab[2];
			ytab[2] = 0.;
			ytab[3] = ytab[0];
			((TypeDataDrawTimePosition*)data)->fillPolygon(((TypeDataDrawTimePosition*)data)->color, ((TypeDataDrawTimePosition*)data)->alpha, xtab, ytab, 4, &(info->param));
		}
	}
}

