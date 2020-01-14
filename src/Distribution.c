#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "Utils.h"
#include "Distribution.h"



int compareDistributionItem(const void* a, const void* b) {
    if(((TypeDistributionItem*)a)->val>((TypeDistributionItem*)b)->val)
        return 1;
    if(((TypeDistributionItem*)a)->val<((TypeDistributionItem*)b)->val)
        return -1;
    return 0;
}

//TypeDistribution meanDistribution(TypeDistribution d1, double w1, TypeDistribution d2, double w2) {
	//TypeDistribution d;
	//if(d1.size != d2.size) {
		//d.size = 0;
		//d.item = NULL;
		//return d;
	//}
	//return d;
//}

double sumDistribution(TypeDistribution d) {
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

TypeDistribution resampleDistribution(TypeDistribution d, double min, double max, int def) {
	int i, cur;
	double step;
	TypeDistribution res;
	res.size = def;
	if(res.size == 0) {
		res.item = NULL;
		return res;
	}
	res.item = (TypeDistributionItem*)malloc(res.size*sizeof(TypeDistributionItem));
	step = (max-min)/((double)res.size-1.);
	cur = 0;
	for(i=0; i<res.size; i++) {
		res.item[i].val = min+((double)i)*step;
		if(cur<d.size-1) {
			while(d.item[cur+1].val<=res.item[i].val)
				cur++;
			if(res.item[i].val>=d.item[cur].val && res.item[i].val<d.item[cur+1].val)
				res.item[i].dens = d.item[cur].dens+((res.item[i].val-d.item[cur].val)*(d.item[cur+1].dens-d.item[cur].dens))/(d.item[cur+1].val-d.item[cur].val);
			else
				res.item[i].dens = 0.;	
		} else
			res.item[i].dens = 1.;
	}
	return res;
}

TypeDistribution agregDistribution(TypeDistribution *d, int size, int def) {
	int i, *cur;
	double min, max, step;
	TypeDistribution res;
	if(size ==0) {
		res.size = 0;
		res.item = NULL;
		return res;
	}
	min = d[0].item[0].val;
	max = d[0].item[d[0].size-1].val;
	for(i=1; i<size; i++) {
		if(d[i].item[0].val<min)
			min = d[i].item[0].val;
		if(d[i].item[d[i].size-1].val>max)
			max = d[i].item[d[i].size-1].val;
	}
	res.size = def;
	res.item = (TypeDistributionItem*)malloc(res.size*sizeof(TypeDistributionItem));
	step = (max-min)/((double)res.size-1.);
	cur = (int*) malloc(size*sizeof(int));
	for(i=0; i<size; i++)
		cur[i] = 0;
	for(i=0; i<res.size; i++) {
		int j;
		res.item[i].val = min+((double)i)*step;
		res.item[i].dens = 0.;
		for(j=0; j<size; j++) {
			if(cur[j]<d[j].size-1) {
				while(d[j].item[cur[j]+1].val<=res.item[i].val)
					cur[j]++;
				if(res.item[i].val>=d[j].item[cur[j]].val && res.item[i].val<d[j].item[cur[j]+1].val)
					res.item[i].dens += d[j].item[cur[j]].dens+((res.item[i].val-d[j].item[cur[j]].val)*(d[j].item[cur[j]+1].dens-d[j].item[cur[j]].dens))/(d[j].item[cur[j]+1].val-d[j].item[cur[j]].val);
			} else
				res.item[i].dens += 1.;
		}
	}
	free((void*)cur);
	for(i=0; i<res.size; i++)
		res.item[i].dens /= (double) size;
	return res;
}

#define MAX_SIZE_TMP 50
#define INC_BUFFER 50
#define IS_SEP(c) (c == ' ' || c == '\t' || c == '\r' || c == '\n' || c == ';')

TypeDistribution readDistribution(FILE *f) {
	char c, tmpA[MAX_SIZE_TMP+1], tmpB[MAX_SIZE_TMP+1];
	int sizeBuffer;
	TypeDistribution d;
	sizeBuffer = INC_BUFFER;
	d.item = (TypeDistributionItem*) malloc(sizeBuffer*sizeof(TypeDistributionItem));
	d.size = 0;
	do {
		c = getc(f);
	} while(c!=EOF && IS_SEP(c)); 
	while(c != EOF) {
		int i;
		i = 0;
		while(i<MAX_SIZE_TMP && c !=EOF && !IS_SEP(c)) {
			tmpA[i] = c;
			c = getc(f);
			i++;
		}
		tmpA[i++] = '\0';
		if(i == MAX_SIZE_TMP) {
			fprintf(stderr, "Ident too long (%s) ...", tmpA);
			exit(1);
		}
		if(i <= 1) {
			fprintf(stderr, "Problem empty (%s) ...", tmpA);
			exit(1);
		}
		while(c!=EOF && IS_SEP(c))
			c=getc(f);
		i = 0;
		while(i<MAX_SIZE_TMP && c !=EOF && !IS_SEP(c)) {
			tmpB[i] = c;
			c = getc(f);
			i++;
		}
		tmpB[i++] = '\0';
		if(i == MAX_SIZE_TMP) {
			fprintf(stderr, "Ident too long (%s) ...", tmpB);
			exit(1);
		}
		if(i <= 1) {
			fprintf(stderr, "Problem empty (%s) ...", tmpB);
			exit(1);
		}
		if(d.size >= sizeBuffer) {
			sizeBuffer += INC_BUFFER;
			d.item = (TypeDistributionItem*) realloc((void *) d.item, sizeBuffer*sizeof(TypeDistributionItem));
		}
		d.item[d.size].val = atof(tmpA);
		d.item[d.size].dens = atof(tmpB);
		d.size++;
		while(c!=EOF && IS_SEP(c))
			c=getc(f);
	}
	d.item = (TypeDistributionItem*) realloc((void *) d.item, d.size*sizeof(TypeDistributionItem));
	return d;
}

void fprintDistribution(FILE *f, TypeDistribution d) {
	int i;
	for(i=0; i<d.size; i++)
		fprintf(f, "%lf\t%le\n", d.item[i].val, d.item[i].dens);
}

void fprintDerive(FILE *f, TypeDistribution d) {
	int i;
	for(i=0; i<d.size-1; i++)
		fprintf(f, "%lf\t%le\n", d.item[i].val, d.item[i+1].dens- d.item[i].dens);
}

//void deriveDistribution(TypeDistribution *d) {
	//if(d->size>0) {
		//int i;
		//for(i=0; i<d->size-1; i++)
			//d->item[i].dens = d->item[i+1].dens-d->item[i].dens;
		//d->size--;
	//}
//}

//void deriveDistribution(TypeDistribution *d) {
	//if(d->size>0) {
		//int i;
		//for(i=0; i<d->size-1; i++) {
			//d->item[i].dens = (d->item[i+1].dens-d->item[i].dens)/(d->item[i+1].val-d->item[i].val);
			//d->item[i].val = (d->item[i].val+d->item[i+1].val)/2.;
		//}
		//d->item[d->size-1].dens = 0.;
		//d->item[d->size-1].val = d->item[d->size-1].val+(d->item[d->size-1].val-d->item[d->size-2].val)/2.;
	//}
//}

//void deriveDistribution(TypeDistribution *d) {
	//if(d->size>0) {
		//int i;
		//for(i=0; i<d->size-1; i++) {
			//d->item[i].dens = (d->item[i+1].dens-d->item[i].dens)/(d->item[i+1].val-d->item[i].val);
			//d->item[i].val = (d->item[i].val+d->item[i+1].val)/2.;
		//}
		//d->item[d->size-1].dens = 2*(1.-d->item[d->size-1].dens)/(d->item[d->size-1].val-d->item[d->size-2].val);
		//d->item[d->size-1].val = d->item[d->size-1].val;
	//}
//}

void deriveDistribution(TypeDistribution *d) {
	if(d->size>0) {
		int i;
		double tmp0 = d->item[0].dens;
		d->item[0].dens = (d->item[0].dens)/(d->item[1].val-d->item[0].val);
		for(i=1; i<d->size; i++) {
			double tmp1 = d->item[i].dens;
			d->item[i].dens = (tmp1-tmp0)/(d->item[i].val-d->item[i-1].val);
			tmp0 = tmp1;
		}
	}
}

void integrateDistribution(TypeDistribution *d) {
	if(d->size>0) {
		int i;
		d->item[0].dens = (d->item[0].dens)*(d->item[1].val-d->item[0].val);
		for(i=1; i<d->size; i++)
			d->item[i].dens = d->item[i-1].dens+(d->item[i].dens)*(d->item[i].val-d->item[i-1].val);
	}
}

double getMean(TypeDistribution d) {
	int i;
	double sum=0.;
	for(i=0; i<d.size-1; i++)
		sum += d.item[i].val*(d.item[i+1].dens-d.item[i].dens);
	return sum;
}

double getMode(TypeDistribution d) {
	int i, imax=0;
	double max=NEG_INFTY;
	for(i=0; i<d.size-1; i++)
		if((d.item[i+1].dens-d.item[i].dens)>max) {
			max = d.item[i+1].dens-d.item[i].dens;
			imax = i;
		}
	return d.item[imax].val;
}

double getMedianDens(TypeDistribution d) {
	int i;
	for(i=0; i<d.size && d.item[i].dens <= 0.5; i++)
		;
	return d.item[i-1].val;
}

double getMeanDens(TypeDistribution d) {
	int i;
	double sum=0.;
	for(i=0; i<d.size; i++)
		sum += d.item[i].val*d.item[i].dens;
	return sum;
}

double getMedian(TypeDistribution d) {
	int i;
	if(d.size>1) {
		for(i=0; i<d.size && d.item[i].dens <= 0.5; i++)
			;
		return d.item[i-1].val;
	} else {
		if(d.size>0)
			return d.item[0].val;
		else
			return 0.;
	}
}

double getQuantileInf(TypeDistribution d, double q) {
	int i;
	if(d.size <= 1)
		return 0;
	for(i=0; i<d.size && d.item[i].dens<=q; i++)
		;
	return d.item[i-1].val;
}

double getQuantileSup(TypeDistribution d, double q) {
	int i;
	if(d.size <= 1)
		return 0;
	if(d.size == 1)
		return d.item[0].dens;
	for(i=d.size-1; i>=0 && d.item[i].dens>=(1.-q); i--)
		;
	return d.item[i+1].val;
}	
