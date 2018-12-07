#include "DrawFossilInt.h"

void drawFossilInt(int n , double x, double y, TypeInfoDrawTreeGeneric *info, void *data) {
	if(data && ((TypeDataDrawFossilInt*)data)->fos) {
		int f;
        for(f=((TypeDataDrawFossilInt*)data)->fos->fossilInt[n]; f>=0; f=((TypeDataDrawFossilInt*)data)->fos->fossilIntList[f].prec)
            ((TypeDataDrawFossilInt*)data)->drawLineDot(
            ((TypeDataDrawFossilInt*)data)->color, 
            ((TypeDataDrawFossilInt*)data)->alpha, 
            ((TypeDataDrawFossilInt*)data)->radius, 
            (((TypeDataDrawFossilInt*)data)->fos->fossilIntList[f].fossilInt.inf-info->param.tmin)*info->param.scale+info->param.xoffset, 
            info->param.yoffset+y, 
            (((TypeDataDrawFossilInt*)data)->fos->fossilIntList[f].fossilInt.sup-info->param.tmin)*info->param.scale+info->param.xoffset, 
            info->param.yoffset+y, 
            &(info->param));
	}
}

