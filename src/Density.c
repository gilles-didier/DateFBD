#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "Utils.h"
#include "Density.h"



int compareDensityItem(const void* a, const void* b) {
    if(((TypeDensityItem*)a)->val>((TypeDensityItem*)b)->val)
        return 1;
    if(((TypeDensityItem*)a)->val<((TypeDensityItem*)b)->val)
        return -1;
    return 0;
}

TypeDensity meanDensity(TypeDensity d1, double w1, TypeDensity d2, double w2) {
	int i;
	TypeDensity d;
	if(d1.size != d2.size) {
		d.size = 0;
		d.item = NULL;
		return d;
	}
	return d;
}

double sumDensity(TypeDensity d) {
	int i;
	double sum;
	if(d.size == 0)
		return 0;
	if(d.size == 1)
		return d.item[0].dens;
	sum = (d.item[1].val-d.item[0].val)*d.item[0].dens;
	for(i=1; i<d.size-1; i++)
		sum += ((d.item[i].val-d.item[i-1].val)/2.+(d.item[i+1].val-d.item[i].val)/2.)*d.item[i].dens;
	sum += (d.item[d.size-1].val-d.item[d.size-2].val)*d.item[d.size-1].dens;
	return sum;
}

void fprintDensity(FILE *f, TypeDensity d) {
	int i;
	for(i=0; i<d.size; i++)
		fprintf(f, "%.2lf\t%le\n", d.item[i].val, d.item[i].dens);
}
