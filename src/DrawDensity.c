#include "DrawDensity.h"

#define PER 0.05

void drawDensity(int n , double x, double y, TypeInfoDrawTreeGeneric *info, void *data) {
	if(data && ((TypeDataDrawDensity*)data)->dens && ((TypeDataDrawDensity*)data)->dens[n].size>0) {
		double xscale, yscale, xoffset, yoffset, *xtab, *ytab, max;
		int i, ind;
		yoffset = y+(PER-0.5)*info->param.leafSep;
		xscale = info->param.scale;
		xoffset = info->param.xoffset;
		max = ((TypeDataDrawDensity*)data)->dens[n].item[0].dens;
		for(i=1; i<((TypeDataDrawDensity*)data)->dens[n].size; i++)
			if(((TypeDataDrawDensity*)data)->dens[n].item[i].dens>max)
				max = ((TypeDataDrawDensity*)data)->dens[n].item[i].dens;
		yscale = (1-2.*PER)*info->param.leafSep/max;
		xtab = (double*) malloc((((TypeDataDrawDensity*)data)->dens[n].size+2)*sizeof(double));
		ytab = (double*) malloc((((TypeDataDrawDensity*)data)->dens[n].size+2)*sizeof(double));
		ind = 0;
		for(i=0; i<((TypeDataDrawDensity*)data)->dens[n].size && ((TypeDataDrawDensity*)data)->dens[n].item[i].dens*yscale<0.001; i++)
			;
		xtab[ind] = (((TypeDataDrawDensity*)data)->dens[n].item[i].val-info->param.tmin)*xscale+xoffset;
		ytab[ind] = (1-2.*PER)*info->param.leafSep+yoffset;
		ind++;
		for(; i<((TypeDataDrawDensity*)data)->dens[n].size; i++) {
			xtab[ind] = (((TypeDataDrawDensity*)data)->dens[n].item[i].val-info->param.tmin)*xscale+xoffset;
			ytab[ind] = (1-2.*PER)*info->param.leafSep-((TypeDataDrawDensity*)data)->dens[n].item[i].dens*yscale+yoffset;
			ind++;
		}
		xtab[ind] = (((TypeDataDrawDensity*)data)->dens[n].item[((TypeDataDrawDensity*)data)->dens[n].size-1].val-info->param.tmin)*xscale+xoffset;
		ytab[ind] = (1-2.*PER)*info->param.leafSep+yoffset;
		ind++;
		((TypeDataDrawDensity*)data)->fillPolygon(((TypeDataDrawDensity*)data)->color, ((TypeDataDrawDensity*)data)->alpha, xtab, ytab, ind, &(info->param));
		free((void*)xtab);
		free((void*)ytab);
	}
}

