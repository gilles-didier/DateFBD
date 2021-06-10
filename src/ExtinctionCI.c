#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include "Utils.h"
#include "Tree.h"
#include "ExtinctionCI.h"


double maxIntervalStraussSadler(double risk, int n, TypeTree *tree, TypeFossilFeature *fos) {
	double min, max;
	int f, tot;
min = fos->fossilList[fos->fossil[n]].time;
max = fos->fossilList[fos->fossil[n]].time;
	if(fos->fossil[n] == NOSUCH || fos->fossilList[fos->fossil[n]].prec == NOSUCH) {
//		printf("%d\t%lf\t%lf\tnan\t", 1, min, max);
		return MY_NAN;
	}
	min = fos->fossilList[fos->fossil[n]].time;
	max = fos->fossilList[fos->fossil[n]].time;
	tot=1;
	for(f=fos->fossilList[fos->fossil[n]].prec; f!=NOSUCH; f=fos->fossilList[f].prec) {
		tot++;
		if(fos->fossilList[f].time<min)
			min = fos->fossilList[f].time;
		if(fos->fossilList[f].time>max)
			max = fos->fossilList[f].time;
	}
//printf("%d\t%lf\t%lf\t%lf\t", tot, min, max, (pow(risk, -1./((double)tot-1.))-1.)*(max-min));
	return max+(pow(risk, -1./((double)tot-1.))-1.)*(max-min);
//	return (((double)tot)*max-min)/((double)tot-1.);
}

double maxLikelihoodStraussSadler(int n, TypeTree *tree, TypeFossilFeature *fos) {
	double min, max;
	int f, tot;
	if(fos->fossil[n] == NOSUCH || fos->fossilList[fos->fossil[n]].prec == NOSUCH)
		return MY_NAN;
	min = fos->fossilList[fos->fossil[n]].time;
	max = fos->fossilList[fos->fossil[n]].time;
	tot=1;
	for(f=fos->fossilList[fos->fossil[n]].prec; f!=NOSUCH; f=fos->fossilList[f].prec) {
		tot++;
		if(fos->fossilList[f].time<min)
			min = fos->fossilList[f].time;
		if(fos->fossilList[f].time>max)
			max = fos->fossilList[f].time;
	}
	return (((double)tot)*max-min)/((double)tot-1.);
}
