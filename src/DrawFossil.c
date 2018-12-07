#include "DrawFossil.h"

void drawFossil(int n , double x, double y, TypeInfoDrawTreeGeneric *info, void *data) {
	if(data && ((TypeDataDrawFossil*)data)->fos) {
		int f;
//printf("y %.2lf\n", y);
        for(f=((TypeDataDrawFossil*)data)->fos->fossil[n]; f>=0; f=((TypeDataDrawFossil*)data)->fos->fossilList[f].prec)
            ((TypeDataDrawFossil*)data)->drawDot(((TypeDataDrawFossil*)data)->color, ((TypeDataDrawFossil*)data)->alpha, ((TypeDataDrawFossil*)data)->radius, (((TypeDataDrawFossil*)data)->fos->fossilList[f].time-info->param.tmin)*info->param.scale+info->param.xoffset, info->param.yoffset+y, &(info->param));
	}
}


