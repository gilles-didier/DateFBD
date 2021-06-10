#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "Utils.h"
#include "TreeExtras.h"
#include "Uncertainty.h"

typedef enum UNCERTAINTY_LIKE_TYPE {
	UncertaintyTypeX=0,
	UncertaintyTypeA,
	UncertaintyTypeB,
	UncertaintyTypeC,
	UncertaintyTypeD,
	UncertaintyTypeBbase,
	UncertaintyTypeCbase,
	UncertaintyTypeDbase
} TypeUncertaintyLikeType;

typedef struct UNCERTAINTY_LR_DATA {
	int nLeaves;
	double logRank;
} TypeUncertaintyLogRankData;

typedef struct UNCERTAINTY_LIKE_DATA_A {
	int numberTaxa;
	double logRank;
} TypeUncertaintyLikeDataA;

typedef struct UNCERTAINTY_LIKE_DATA_BCD {
	int min, size, nLeafMin;
	double stopTime, *like;
} TypeUncertaintyLikeDataBCD;

typedef union UNCERTAINTY_LIKE_DATA_X {
	TypeUncertaintyLikeDataA A;
	TypeUncertaintyLikeDataBCD BCD;
} TypeUncertaintyLikeDataX;

typedef struct UNCERTAINTY_LIKE_DATA {
	TypeUncertaintyLikeType type;
	int nLeaves;
	TypeUncertaintyLikeDataX dataStd;
} TypeUncertaintyLikeData;

typedef struct UNCERTAINTY_LIKE_PARAMETERS {
    double birth, death, fossil, sampling;
    double alpha, beta, omega;
    double ab, bma, bmab, amab;
    double bs, as;
} TypeUncertaintyLikeParameters;

static double getLogLikelihoodSplitted(TypeTree *tree, TypeUncertaintyLikeParameters *param);
static int fillSplitTreeFossil(TypeTree *tcur, int n, TypeTree *tree, TypeFossilTab *ftab, TypeTree **treeList, int *size);
static double getLogLikelihoodSubtree(double t, int n, TypeTree *tree, TypeUncertaintyLikeData *data, TypeUncertaintyLikeParameters *param);
static double getLogLikelihoodSubtreeRoot(double a, double b, int n, TypeTree *tree, TypeUncertaintyLikeData *data, TypeUncertaintyLikeParameters *param);
static void fillIndexMin(int n, TypeTree *tree, int *indexMin);
static TypeUncertaintyLikeParameters getUncertaintyLikeParameters(TypeModelParam *param);
static double getLogIntProbTypeA(double a, double b, double maxTime, int k,  TypeUncertaintyLikeParameters *param);
static double getLogIntProbTypeC(double a, double b, double endTime, double maxTime, int k, int l, TypeUncertaintyLikeParameters *param);
static double getLogIntProbTypeD(double a, double b, double endTime, double maxTime, int k, TypeUncertaintyLikeParameters *param);
static double getLogIntProbTypeB(double a, double b, double endTime, double maxTime, int k, TypeUncertaintyLikeParameters *param);
static double getLogProbTypeA(double startTime, double maxTime, int k,  TypeUncertaintyLikeParameters *param);
static double getLogProbTypeB(double startTime, double endTime, double maxTime, int k, TypeUncertaintyLikeParameters *param);
static double getLogProbTypeC(double startTime, double endTime, double maxTime, int k, int l, TypeUncertaintyLikeParameters *param);
static double getLogProbTypeD(double startTime, double endTime, double maxTime, int k, TypeUncertaintyLikeParameters *param);
static double getLogProbNotObservableU(double t, double maxTime,  TypeUncertaintyLikeParameters *param);
static void fillTableTreeIter(int n, TypeTree *tree, int *indexMin, double *K, double *T);
static int fillBasicTree(TypeTree *basic, int n, TypeTree *tree, TypeFossilTab *ftab, int toCompute, int *newIndex);
static TypeTree *getBasicTreeDivergence(int n, TypeTree *tree, TypeFossilTab *ftab, int *newIndex, int *root);
static double getMinTimeTree(int n, TypeTree *tree, TypeFossilTab *ftab);
static void freeUncertaintyLikeData(TypeUncertaintyLikeData* data);
static void fillUncertaintyLikeData(int n, TypeTree *tree, int *indexMin, TypeUncertaintyLogRankData *lr, TypeUncertaintyLikeData *data, TypeUncertaintyLikeParameters *param);
static TypeUncertaintyLikeDataBCD getUncertaintyLikeDataBCD(int n, double val, TypeTree *tree, TypeUncertaintyLikeData *data, int *indexMin, TypeUncertaintyLikeParameters *param);
static TypeUncertaintyLikeDataBCD joinUncertaintyLikeDataBCD(double likeOne, TypeUncertaintyLikeDataBCD d0, TypeUncertaintyLikeDataBCD d1, TypeUncertaintyLikeParameters *param);
static void fillUncertaintyLogRankData(int node, TypeTree *tree, TypeUncertaintyLogRankData *lr);
static double logBinomial(unsigned int k, unsigned int n);
static double logFactorial(unsigned int n);
static double myLog(double a,double b);

double logFactorial(unsigned int n) {
	return lgammafn((double) n+1.);
}

double logBinomial(unsigned int k, unsigned int n) {
    return logFactorial(n)-(logFactorial(k)+logFactorial(n-k));
}

void fillUncertaintyLogRankData(int node, TypeTree *tree, TypeUncertaintyLogRankData *lr) {
	if(tree->node[node].child == NOSUCH) {
		lr[node].nLeaves = 1;
		lr[node].logRank = 0.;
	} else {
		fillUncertaintyLogRankData(tree->node[node].child, tree, lr);
		fillUncertaintyLogRankData(tree->node[tree->node[node].child].sibling, tree, lr);		
		lr[node].nLeaves = lr[tree->node[node].child].nLeaves+lr[tree->node[tree->node[node].child].sibling].nLeaves;
		lr[node].logRank = lr[tree->node[node].child].logRank+lr[tree->node[tree->node[node].child].sibling].logRank+log(2)+logBinomial((unsigned int) (lr[tree->node[node].child].nLeaves-1), (unsigned int) (lr[node].nLeaves-2));
	}
}

TypeUncertaintyLikeDataBCD joinUncertaintyLikeDataBCD(double likeOne, TypeUncertaintyLikeDataBCD d0, TypeUncertaintyLikeDataBCD d1, TypeUncertaintyLikeParameters *param) {
	TypeUncertaintyLikeDataBCD res;
	int i, j, inc;
	double *like, *offset;
	if(d0.stopTime != d1.stopTime) {
		error("Execution error: call of joinUncertaintyLikeDataBCD with different stop times (d0 %.2lf d1 %.2lf)\n", d0.stopTime, d1.stopTime);
	}
	if(likeOne != NEG_INFTY) {
		res.min = 1;
		res.size = d0.size+d0.min+d1.size+d1.min-2;
	} else {
		res.min = d0.min+d1.min;
		res.size = d0.size+d1.size-1;
	}
	res.stopTime = d0.stopTime;
	res.nLeafMin = d0.nLeafMin+d1.nLeafMin;
	res.like = (double*) malloc(res.size*sizeof(double));
	like = (double*) malloc(res.size*sizeof(double));
	offset = (double*) malloc(res.size*sizeof(double));
	for(i=0; i<res.size; i++) {
		offset[i] = NEG_INFTY;
		like[i] = 0.;
	}
	if(likeOne != NEG_INFTY)
		res.like[0] = likeOne;
	inc = d0.min+d1.min-res.min;
	for(i=0; i<d0.size; i++)
		if(d0.like[i] != NEG_INFTY) {
			for(j=0; j<d1.size; j++)
				if(d1.like[j] != NEG_INFTY) {
					double logLike = d0.like[i]+d1.like[j]+log(2.)+logBinomial((unsigned int) i+d0.min-1, (unsigned int) (i+j+d0.min+d1.min-2));
					if(offset[i+j+inc] == NEG_INFTY) {
						offset[i+j+inc] = logLike;
						like[i+j+inc] = 1.;
					} else {
						if(logLike>offset[i+j+inc]) { /*compute max in offset just to avoid numerical precision issues*/
							like[i+j+inc] *= exp(offset[i+j+inc]-logLike);
							offset[i+j+inc] = logLike;
							like[i+j+inc]++;
						} else
							like[i+j+inc] += exp(logLike-offset[i+j+inc]);
					}
				}
		}
	for(i=(likeOne != NEG_INFTY)?1:0; i<res.size; i++)
		if(like[i]>0.)
			res.like[i] = log(like[i])+offset[i];
		else
			res.like[i] = NEG_INFTY;
	free((void*)like);
	free((void*)offset);
	return res;
}

TypeUncertaintyLikeDataBCD getUncertaintyLikeDataBCD(int n, double val, TypeTree *tree, TypeUncertaintyLikeData *data, int *indexMin, TypeUncertaintyLikeParameters *param) {
	TypeUncertaintyLikeDataBCD res;
	if(tree->node[n].child == NOSUCH) {
		res.stopTime = val;
		res.min = 1;
		res.size = 1;
		res.like = (double*) malloc(sizeof(double));
		if(tree->time[n] == val) {
			res.nLeafMin = 1;
			res.like[0] = 0.;
		} else {
			res.nLeafMin = 0;		
			res.like[0] = getLogLikelihoodSubtree(val, n, tree, data, param);
		}
	} else {
		TypeUncertaintyLikeDataBCD d0, d1;
		d0 = getUncertaintyLikeDataBCD(tree->node[n].child, val, tree, data, indexMin, param);
		d1 = getUncertaintyLikeDataBCD(tree->node[tree->node[n].child].sibling, val, tree, data, indexMin, param);
		if(tree->time[indexMin[n]] > val)
			res = joinUncertaintyLikeDataBCD(getLogLikelihoodSubtree(val, n, tree, data, param), d0, d1, param);
		else
			res = joinUncertaintyLikeDataBCD(NEG_INFTY, d0, d1, param);		
		free((void*)d0.like);
		free((void*)d1.like);
	}
	return res;
}


void fillUncertaintyLikeData(int n, TypeTree *tree, int *indexMin, TypeUncertaintyLogRankData *lr, TypeUncertaintyLikeData *data, TypeUncertaintyLikeParameters *param) {
	if(tree->time[indexMin[n]] < tree->maxTime) { /*meaning it can happen something after tree->time[indexMin[n]], i.e. types B,C,D*/
 		if(tree->node[n].child != NOSUCH) {
			int l, child[2];
			TypeUncertaintyLikeDataBCD d[2];
			child[0] = tree->node[n].child;
			child[1] = tree->node[child[0]].sibling;
			for(l=0; l<2; l++)
				fillUncertaintyLikeData(child[l], tree, indexMin, lr, data, param);
			data[n].dataStd.BCD.stopTime = tree->time[indexMin[n]];
			switch(((TypeNodeStatus*)tree->info)[indexMin[n]]) {
				case extinctNodeStatus:
					data[n].type = UncertaintyTypeB;
					break;
				case divergenceNodeStatus:
					data[n].type = UncertaintyTypeD;
					break;
				default:
					data[n].type = UncertaintyTypeC;
			}
			for(l=0; l<2; l++)
				if(tree->time[indexMin[child[l]]] == tree->time[indexMin[n]])
					d[l] = data[child[l]].dataStd.BCD;
				else 
					d[l] = getUncertaintyLikeDataBCD(child[l], tree->time[indexMin[n]], tree, data, indexMin, param);
			data[n].dataStd.BCD = joinUncertaintyLikeDataBCD(NEG_INFTY, d[0], d[1], param);
			for(l=0; l<2; l++)
				if(tree->time[indexMin[child[l]]] > tree->time[indexMin[n]] && d[l].like != NULL)
					free((void*)d[l].like);
		} else {
			switch(((TypeNodeStatus*)tree->info)[indexMin[n]]) {
				case extinctNodeStatus:
					data[n].type = UncertaintyTypeBbase;
					break;
				case divergenceNodeStatus:
					data[n].type = UncertaintyTypeDbase;
					break;
				default:
					data[n].type = UncertaintyTypeCbase;
			}
			data[n].dataStd.BCD.stopTime = tree->time[indexMin[n]];
			data[n].dataStd.BCD.min = 1;
			data[n].dataStd.BCD.size = 1;
			data[n].dataStd.BCD.like = (double*) malloc(sizeof(double));
			data[n].dataStd.BCD.like[0] = 0.;
			data[n].dataStd.BCD.nLeafMin = 1;
		}
	} else {
 		if(tree->node[n].child != NOSUCH) {
			int l, child[2];
			child[0] = tree->node[n].child;
			child[1] = tree->node[child[0]].sibling;
			for(l=0; l<2; l++)
				fillUncertaintyLikeData(child[l], tree, indexMin, lr, data, param);
		}
		data[n].type = UncertaintyTypeA;
		data[n].dataStd.A.numberTaxa = lr[n].nLeaves;
		data[n].dataStd.A.logRank = lr[n].logRank;
	}
}

double getLogDensitySpecialX(double *tot, double time, int n, TypeTree *tree, TypeFossilFeature *fos, TypeModelParam *param) {
	int index, i, j, root, *indexMin;
	TypeTree *basicU;
	double maxOrigin, lltot, llpart, diff, res;
	TypeUncertaintyLikeParameters pu = getUncertaintyLikeParameters(param);
	TypeUncertaintyLikeData *data;
	TypeUncertaintyLogRankData *lr;
	lltot = getLogLikelihoodTreeFossil(tree, fos, param);
	TypeFossilTab *ftab = listToFossilTab(fos, tree->size);
	basicU = getBasicTreeDivergence(n, tree, ftab, &index, &root);
	for(i=0; i<basicU->size; i++)
		if(basicU->node[i].child == NOSUCH && basicU->time[i] == NO_TIME)
			warning("Execution warning %d %s (Uncertainty.c:197)\n", i, tree->name[i]);
	if(basicU->parent == NULL)
		setParent(basicU);
	if(basicU->minTimeInt.sup != NO_TIME)
		maxOrigin = basicU->minTimeInt.sup;
	else
		maxOrigin = getMinTimeTree(basicU->root, basicU, NULL);
	indexMin = (int*) malloc(basicU->size*sizeof(int));
	fillIndexMin(basicU->root, basicU, indexMin);
	lr = (TypeUncertaintyLogRankData*) malloc(basicU->size*sizeof(TypeUncertaintyLogRankData));
	fillUncertaintyLogRankData(basicU->root, basicU, lr);
	data = (TypeUncertaintyLikeData*) malloc(basicU->size*sizeof(TypeUncertaintyLikeData));
	fillUncertaintyLikeData(basicU->root, basicU, indexMin, lr, data, &pu);
	if(basicU->minTimeInt.inf<maxOrigin)
		llpart = getLogLikelihoodSubtreeRoot(basicU->minTimeInt.inf, maxOrigin, basicU->root, basicU, data, &pu);
	else
		llpart = getLogLikelihoodSubtree(basicU->minTimeInt.inf, basicU->root, basicU, data, &pu);
	diff = lltot-llpart;
	for(j=0; j<basicU->size; j++)
		freeUncertaintyLikeData(&(data[j]));
	((TypeNodeStatus*)basicU->info)[index] = divergenceNodeStatus;
	double end = utils_MIN(time, maxOrigin);
	basicU->time[index] = time;
	fillIndexMin(tree->root, basicU, indexMin);
	fillUncertaintyLikeData(basicU->root, basicU, indexMin, lr, data, &pu);
	res = getLogLikelihoodSubtreeRoot(basicU->minTimeInt.inf, end, basicU->root, basicU, data, &pu)+diff;
	for(j=0; j<basicU->size; j++)
		freeUncertaintyLikeData(&(data[j]));
	for(i=0; i<tree->size; i++)
		if(ftab[i].time != NULL)
			free((void*)ftab[i].time);
	free((void*)ftab);
	if(basicU->info!=NULL)
		free((void*)basicU->info);
	freeTree(basicU);
	free((void*)indexMin);
	free((void*)lr);
	free((void*)data);
	*tot = lltot;
	return res;
}

//#define PRECISION_LOG_LIKE 6*log(10.)
#define PRECISION_LOG_LIKE -5*log(10.)

void fillLogDistribution(TypeDistribution *logD, double *logCond, int n, TypeTree *tree, TypeFossilFeature *fos, TypeModelParam *param) {
	int index, i, j, root, *indexMin;
	TypeTree *basicU;
	double max, maxOrigin, lltot, llpart, diff;
	TypeUncertaintyLikeParameters pu = getUncertaintyLikeParameters(param);
	TypeUncertaintyLikeData *data;
	TypeUncertaintyLogRankData *lr;
	lltot = getLogLikelihoodTreeFossil(tree, fos, param);
	TypeFossilTab *ftab = listToFossilTab(fos, tree->size);
	basicU = getBasicTreeDivergence(n, tree, ftab, &index, &root);
	for(i=0; i<basicU->size; i++)
		if(basicU->node[i].child == NOSUCH && basicU->time[i] == NO_TIME)
			warning("Execution warning %d %s (Uncertainty.c:197)\n", i, tree->name[i]);
	if(basicU->parent == NULL)
		setParent(basicU);
	max = getMinTimeTree(index, basicU, NULL);
	if(basicU->minTimeInt.sup != NO_TIME)
		maxOrigin = basicU->minTimeInt.sup;
	else
		maxOrigin = getMinTimeTree(basicU->root, basicU, NULL);
	indexMin = (int*) malloc(basicU->size*sizeof(int));
	fillIndexMin(basicU->root, basicU, indexMin);
	lr = (TypeUncertaintyLogRankData*) malloc(basicU->size*sizeof(TypeUncertaintyLogRankData));
	fillUncertaintyLogRankData(basicU->root, basicU, lr);
	data = (TypeUncertaintyLikeData*) malloc(basicU->size*sizeof(TypeUncertaintyLikeData));
	fillUncertaintyLikeData(basicU->root, basicU, indexMin, lr, data, &pu);
	if(basicU->minTimeInt.inf<maxOrigin)
		llpart = getLogLikelihoodSubtreeRoot(basicU->minTimeInt.inf, maxOrigin, basicU->root, basicU, data, &pu);
	else
		llpart = getLogLikelihoodSubtree(basicU->minTimeInt.inf, basicU->root, basicU, data, &pu);
	diff = lltot-llpart;
	for(j=0; j<basicU->size; j++)
		freeUncertaintyLikeData(&(data[j]));
	((TypeNodeStatus*)basicU->info)[index] = divergenceNodeStatus;
	for(i=0; i<logD->size && logD->item[i].val<basicU->minTimeInt.inf; i++)
		logD->item[i].dens = NEG_INFTY;
	int a=i, b=logD->size-1, m;
	while(b>a+1) {
		m = (a+b)/2;
		if(logD->item[m].val<max)
			a = m;
		else
			b = m;
	}
	b=a; a=i;
	while(b>a+1) {
		double tmp, end;
		m = (a+b)/2;
		end = utils_MIN(logD->item[m].val, maxOrigin);
		basicU->time[index] = logD->item[m].val;
		fillIndexMin(tree->root, basicU, indexMin);
		fillUncertaintyLikeData(basicU->root, basicU, indexMin, lr, data, &pu);
		tmp = getLogLikelihoodSubtreeRoot(basicU->minTimeInt.inf, end, basicU->root, basicU, data, &pu);
		if(tmp<llpart+PRECISION_LOG_LIKE)
			a = m;
		else
			b = m;
		for(j=0; j<basicU->size; j++)
			freeUncertaintyLikeData(&(data[j]));
	}
	for(; i<a; i++)
		logD->item[i].dens = NEG_INFTY;
	for(; i<logD->size && logD->item[i].val<max; i++) {
		double end = utils_MIN(logD->item[i].val, maxOrigin);
		basicU->time[index] = logD->item[i].val;
		fillIndexMin(tree->root, basicU, indexMin);
		fillUncertaintyLikeData(basicU->root, basicU, indexMin, lr, data, &pu);
		logD->item[i].dens = getLogLikelihoodSubtreeRoot(basicU->minTimeInt.inf, end, basicU->root, basicU, data, &pu)+diff;
		for(j=0; j<basicU->size; j++)
			freeUncertaintyLikeData(&(data[j]));
	}
	for(; i<logD->size; i++)
		logD->item[i].dens = lltot;
	*logCond = lltot;
	for(i=0; i<tree->size; i++)
		if(ftab[i].time != NULL)
			free((void*)ftab[i].time);
	free((void*)ftab);
	if(basicU->info!=NULL)
		free((void*)basicU->info);
	freeTree(basicU);
	free((void*)indexMin);
	free((void*)lr);
	free((void*)data);
}


void fillDistribution(TypeDistribution *d, int n, TypeTree *tree, TypeFossilFeature *fos, TypeModelParam *param) {
	int index, i, j, root, *indexMin;
	TypeTree *basicU;
	double max, logDensA;;
	TypeUncertaintyLikeParameters pu = getUncertaintyLikeParameters(param);
	TypeUncertaintyLikeData *data;
	TypeUncertaintyLogRankData *lr;
	TypeFossilTab *ftab = listToFossilTab(fos, tree->size);
	basicU = getBasicTreeDivergence(n, tree, ftab, &index, &root);
	for(i=0; i<basicU->size; i++)
		if(basicU->node[i].child == NOSUCH && basicU->time[i] == NO_TIME)
			warning("Execution warning %d %s (Uncertainty.c:197)\n", i, tree->name[i]);
	if(basicU->parent == NULL)
		setParent(basicU);
	max = getMinTimeTree(index, basicU, NULL);
	indexMin = (int*) malloc(basicU->size*sizeof(int));
	fillIndexMin(basicU->root, basicU, indexMin);
	lr = (TypeUncertaintyLogRankData*) malloc(basicU->size*sizeof(TypeUncertaintyLogRankData));
	fillUncertaintyLogRankData(basicU->root, basicU, lr);
	data = (TypeUncertaintyLikeData*) malloc(basicU->size*sizeof(TypeUncertaintyLikeData));
	fillUncertaintyLikeData(basicU->root, basicU, indexMin, lr, data, &pu);	
	logDensA = getLogLikelihoodSubtreeRoot(basicU->minTimeInt.inf, basicU->minTimeInt.sup, basicU->root, basicU, data, &pu);
	for(j=0; j<basicU->size; j++)
		freeUncertaintyLikeData(&(data[j]));
	((TypeNodeStatus*)basicU->info)[index] = divergenceNodeStatus;
	
	for(i=0; i<d->size && d->item[i].val<basicU->minTimeInt.sup; i++)
		d->item[i].dens = 0.;
	for(; i<d->size && d->item[i].val<max; i++) {
		double tmp;
		basicU->time[index] = d->item[i].val;
		fillIndexMin(tree->root, basicU, indexMin);
		fillUncertaintyLikeData(basicU->root, basicU, indexMin, lr, data, &pu);
		tmp = getLogLikelihoodSubtreeRoot(basicU->minTimeInt.inf, basicU->minTimeInt.sup, basicU->root, basicU, data, &pu);
		d->item[i].dens = exp(tmp-logDensA);
		for(j=0; j<basicU->size; j++)
			freeUncertaintyLikeData(&(data[j]));
	}
	for(; i<d->size; i++)
		d->item[i].dens = 1.;
	for(i=0; i<tree->size; i++)
		if(ftab[i].time != NULL)
			free((void*)ftab[i].time);
	free((void*)ftab);
	if(basicU->info!=NULL)
		free((void*)basicU->info);
	freeTree(basicU);
	free((void*)indexMin);
	free((void*)lr);
	free((void*)data);
}

TypeDistribution getDistribution(int n, double step, TypeTree *tree, TypeFossilFeature *fos, TypeModelParam *param) {
	int index, i, j, root, *indexMin;
	TypeTree *basicU;
	double min, max, logDensA;;
	TypeUncertaintyLikeParameters pu = getUncertaintyLikeParameters(param);
	TypeUncertaintyLikeData *data;
	TypeUncertaintyLogRankData *lr;
	TypeFossilTab *ftab = listToFossilTab(fos, tree->size);
	TypeDistribution d;
	if(tree->time[n] != NO_TIME) {
		d.size = 1;
		d.item = (TypeDistributionItem*)malloc(sizeof(TypeDistributionItem));
		d.item[0].val = tree->time[n];
		d.item[0].dens = 1.;
		return d;
	}
	basicU = getBasicTreeDivergence(n, tree, ftab, &index, &root);
	for(i=0; i<basicU->size; i++)
		if(basicU->node[i].child == NOSUCH && basicU->time[i] == NO_TIME)
			warning("Execution warning %d %s (Uncertainty.c:197)\n", i, tree->name[i]);
	if(basicU->parent == NULL)
		setParent(basicU);
	min = floor(basicU->minTimeInt.inf/step)*step;
	max = ceil(getMinTimeTree(index, basicU, NULL)/step)*step;
	if(min>max) {
		error("Error: origin time and oldest fossil age overlap.\n");
	}
	d.size = 1+round((max-min)/step);
	d.item = (TypeDistributionItem*)malloc(d.size*sizeof(TypeDistributionItem));
	indexMin = (int*) malloc(basicU->size*sizeof(int));
	fillIndexMin(basicU->root, basicU, indexMin);
	lr = (TypeUncertaintyLogRankData*) malloc(basicU->size*sizeof(TypeUncertaintyLogRankData));
	fillUncertaintyLogRankData(basicU->root, basicU, lr);
	data = (TypeUncertaintyLikeData*) malloc(basicU->size*sizeof(TypeUncertaintyLikeData));
	fillUncertaintyLikeData(basicU->root, basicU, indexMin, lr, data, &pu);
	if(basicU->minTimeInt.inf<basicU->minTimeInt.sup)
		logDensA = getLogLikelihoodSubtreeRoot(basicU->minTimeInt.inf, basicU->minTimeInt.sup, basicU->root, basicU, data, &pu);
	else
		logDensA = getLogLikelihoodSubtree(basicU->minTimeInt.inf, basicU->root, basicU, data, &pu);
	for(j=0; j<basicU->size; j++)
		freeUncertaintyLikeData(&(data[j]));
	((TypeNodeStatus*)basicU->info)[index] = divergenceNodeStatus;
	for(i=0; i<d.size-1; i++) {
		double tmp;
		d.item[i].val = min+((double)i)*step;
		basicU->time[index] = d.item[i].val;
		fillIndexMin(tree->root, basicU, indexMin);
		fillUncertaintyLikeData(basicU->root, basicU, indexMin, lr, data, &pu);
		if(basicU->minTimeInt.inf<basicU->minTimeInt.sup)
			tmp = getLogLikelihoodSubtreeRoot(basicU->minTimeInt.inf, basicU->minTimeInt.sup, basicU->root, basicU, data, &pu);
		else
			tmp = getLogLikelihoodSubtree(basicU->minTimeInt.inf, basicU->root, basicU, data, &pu);
		d.item[i].dens = exp(tmp-logDensA);
		for(j=0; j<basicU->size; j++)
			freeUncertaintyLikeData(&(data[j]));
	}
	d.item[d.size-1].val = max;
	d.item[d.size-1].dens = 1.;
	for(i=0; i<tree->size; i++)
		if(ftab[i].time != NULL)
			free((void*)ftab[i].time);
	free((void*)ftab);
	if(basicU->info!=NULL)
		free((void*)basicU->info);
	freeTree(basicU);
	free((void*)indexMin);
	free((void*)lr);
	free((void*)data);
	return d;
}

/*get the whole likelihood*/
/* by convention if a leaf l is such that
 * - tree->time[l] == NO_TIME => it is extinct (it must have fossil(s)),
 * - tree->time[l] == tree->maxTime => it is contemporary,
 * - otherwise => it is unknown after but observable until tree->time[l]
 */
double getLogLikelihoodTreeFossil(TypeTree *tree, TypeFossilFeature *fos, TypeModelParam *param) {
    TypeTree **treeList;
    int l, size;
    double res = 0.;
    if(param->birth == 0. && tree->size>1)
		return NEG_INFTY;
	if(tree->minTimeInt.sup == NO_TIME)
		tree->minTimeInt.sup = getMinFossilTime(fos);
    TypeUncertaintyLikeParameters pu = getUncertaintyLikeParameters(param);
    splitTreeFossil(tree, fos, &treeList, &size);
    for(l=0; l<size; l++) {
		double tmp = getLogLikelihoodSplitted(treeList[l], &pu);
		if(isnan(tmp)) {
			error("Error getLogLikelihoodTreeFossil  of tree %d l 532\ntime %lf %lf %lf\n", l, treeList[l]->minTimeInt.inf, treeList[l]->minTimeInt.sup, treeList[l]->minTime);
		}
       res += tmp;
      
      if(treeList[l]->info!=NULL)
            free((void*)treeList[l]->info);
        freeTree(treeList[l]);
    }
    free((void*)treeList);
    return res;
}

double getLogLikelihoodTreeFossilDebug(TypeTree *tree, TypeFossilFeature *fos, TypeModelParam *param) {
    TypeTree **treeList;
    int l, size;
    double res = 0.;
    TypeUncertaintyLikeParameters pu = getUncertaintyLikeParameters(param);
    splitTreeFossil(tree, fos, &treeList, &size);


    for(l=0; l<size; l++) {

      res += getLogLikelihoodSplitted(treeList[l], &pu);

        if(treeList[l]->info!=NULL)
            free((void*)treeList[l]->info);
        freeTree(treeList[l]);
    }

    free((void*)treeList);
    return res;
}

/*get the likelihood of a tree ending at each fossil or contemporary time*/
double getLogLikelihoodSplitted(TypeTree *tree, TypeUncertaintyLikeParameters *param) {
    int  *indexMin, j;
	TypeUncertaintyLikeData *data;
	TypeUncertaintyLogRankData *lr;
    double res;
    if(tree == NULL || tree->size == 0) {
        if(param->death>0)
            return getLogProbNotObservableU(tree->minTime, tree->maxTime, param);
        else
            return NEG_INFTY;
    }
	indexMin = (int*) malloc(tree->size*sizeof(int));
	fillIndexMin(tree->root, tree, indexMin);
	lr = (TypeUncertaintyLogRankData*) malloc(tree->size*sizeof(TypeUncertaintyLogRankData));
	fillUncertaintyLogRankData(tree->root, tree, lr);
	data = (TypeUncertaintyLikeData*) malloc(tree->size*sizeof(TypeUncertaintyLikeData));
	fillUncertaintyLikeData(tree->root, tree, indexMin, lr, data, param);

	res = getLogLikelihoodSubtreeRoot(tree->minTimeInt.inf, tree->minTimeInt.sup, tree->root, tree, data, param);
	free((void*)indexMin);
	for(j=0; j<tree->size; j++)
		freeUncertaintyLikeData(&(data[j]));
	free((void*)lr);
	free((void*)data);
    return res;
}

double getLogLikelihoodSubtreeRoot(double a, double b, int n, TypeTree *tree, TypeUncertaintyLikeData *data, TypeUncertaintyLikeParameters *param) {
	int l;
	double like = 0., offset = NEG_INFTY;
	switch(data[n].type) {
		case UncertaintyTypeB:
			if(data[n].dataStd.BCD.nLeafMin>1) {
				warning("Execution error: two simultaneous fossils at time %.2lf (subtree %d)\n", data[n].dataStd.BCD.stopTime, n);
				return NEG_INFTY;
			}
			for(l=0; l<data[n].dataStd.BCD.size; l++) {
				double logLike;
				if(b>a)
					logLike = getLogIntProbTypeB(a, b, data[n].dataStd.BCD.stopTime, tree->maxTime, l+data[n].dataStd.BCD.min-1, param)+data[n].dataStd.BCD.like[l]-logFactorial(l+data[n].dataStd.BCD.min-1);
				else
					logLike = getLogProbTypeB(a, data[n].dataStd.BCD.stopTime, tree->maxTime, l+data[n].dataStd.BCD.min-1, param)+data[n].dataStd.BCD.like[l]-logFactorial(l+data[n].dataStd.BCD.min-1);
				if(isinf(logLike) || isnan(logLike)) {
					like += 0.;
				} else {
					if(offset == NEG_INFTY) {
						offset = logLike;
						like = 1.;
					} else {
						if(logLike>offset) { /*compute max in offset just to avoid numerical precision issues*/
							like *= exp(offset-logLike);
							offset = logLike;
						}
						like += exp(logLike-offset);
					}
				}
			}
			if(like>0.)
				return log(like)+offset;
			else
				return NEG_INFTY;
		case UncertaintyTypeC:
			for(l=0; l<data[n].dataStd.BCD.size; l++) {
				double logLike;
				if(b>a)
					logLike = getLogIntProbTypeC(a, b, data[n].dataStd.BCD.stopTime, tree->maxTime, l+data[n].dataStd.BCD.min, data[n].dataStd.BCD.nLeafMin, param)+data[n].dataStd.BCD.like[l]-logFactorial(l+data[n].dataStd.BCD.min-1);
				else
					logLike = getLogProbTypeC(a, data[n].dataStd.BCD.stopTime, tree->maxTime, l+data[n].dataStd.BCD.min, data[n].dataStd.BCD.nLeafMin, param)+data[n].dataStd.BCD.like[l]-logFactorial(l+data[n].dataStd.BCD.min-1);
				if(isinf(logLike) || isnan(logLike)) {
					like += 0.;
				} else {
					if(offset == NEG_INFTY) {
						offset = logLike;
						like = 1.;
					} else {
						if(logLike>offset) { /*compute max in offset just to avoid numerical precision issues*/
							like *= exp(offset-logLike);
							offset = logLike;
						}
						like += exp(logLike-offset);
					}
				}
			}
			if(like>0.)
				return log(like)+offset;
			else
				return NEG_INFTY;
		case UncertaintyTypeD:
			for(l=0; l<data[n].dataStd.BCD.size; l++) {
				double logLike;
				if(b>a)
					logLike = getLogIntProbTypeD(a, b, data[n].dataStd.BCD.stopTime, tree->maxTime, l+data[n].dataStd.BCD.min, param)+data[n].dataStd.BCD.like[l]-logFactorial(l+data[n].dataStd.BCD.min-1);
				else
					logLike = getLogProbTypeD(a, data[n].dataStd.BCD.stopTime, tree->maxTime, l+data[n].dataStd.BCD.min, param)+data[n].dataStd.BCD.like[l]-logFactorial(l+data[n].dataStd.BCD.min-1);
				if(isinf(logLike) || isnan(logLike)) {
					like += 0.;
				} else {
					if(offset == NEG_INFTY) {
						offset = logLike;
						like = 1.;
					} else {
						if(logLike>offset) { /*compute max in offset just to avoid numerical precision issues*/
							like *= exp(offset-logLike);
							offset = logLike;
						}
						like += exp(logLike-offset);
					}
				}
			}
			if(like>0.)
				return log(like)+offset;
			else
				return NEG_INFTY;
		case UncertaintyTypeA:
			if(b>a)
				return getLogIntProbTypeA(a, b, tree->maxTime, data[n].dataStd.A.numberTaxa, param)+data[n].dataStd.A.logRank-logFactorial(data[n].dataStd.A.numberTaxa-1);
			else
				return getLogProbTypeA(a, tree->maxTime, data[n].dataStd.A.numberTaxa, param)+data[n].dataStd.A.logRank-logFactorial(data[n].dataStd.A.numberTaxa-1);
		case UncertaintyTypeBbase:
			if(b>a)
				return getLogIntProbTypeB(a, b, data[n].dataStd.BCD.stopTime, tree->maxTime, 0, param);
			else
				return getLogProbTypeB(a, data[n].dataStd.BCD.stopTime, tree->maxTime, 0, param);
		case UncertaintyTypeCbase:
			if(b>a)
				return getLogIntProbTypeC(a, b, data[n].dataStd.BCD.stopTime, tree->maxTime, 1, 1, param);
			else
				return getLogProbTypeC(a, data[n].dataStd.BCD.stopTime, tree->maxTime, 1, 1, param);
		case UncertaintyTypeDbase:
			if(b>a)
				return getLogIntProbTypeD(a, b, data[n].dataStd.BCD.stopTime, tree->maxTime, 1, param);
			else
				return getLogProbTypeD(a, data[n].dataStd.BCD.stopTime, tree->maxTime, 1, param);
		default:
			error("\nExecution error in function getLogLikelihoodSubtreeRoot: No type while computing the likelihood of %d\n", n);
	}
	return 0.0;
}

double getLogLikelihoodSubtree(double t, int n, TypeTree *tree, TypeUncertaintyLikeData *data, TypeUncertaintyLikeParameters *param) {
	int l;
	double like = 0., offset = NEG_INFTY;
	switch(data[n].type) {
		case UncertaintyTypeB:
			if(data[n].dataStd.BCD.nLeafMin>1) {
				warning("Execution error: two simultaneous fossils at time %.2lf (subtree %d)\n", data[n].dataStd.BCD.stopTime, n);
				return NEG_INFTY;
			}
			for(l=0; l<data[n].dataStd.BCD.size; l++) {
				double logLike;
				logLike = getLogProbTypeB(t, data[n].dataStd.BCD.stopTime, tree->maxTime, l+data[n].dataStd.BCD.min-1, param)+data[n].dataStd.BCD.like[l]-logFactorial(l+data[n].dataStd.BCD.min-1);
				if(isinf(logLike) || isnan(logLike)) {
					like += 0.;
				} else {
					if(offset == NEG_INFTY) {
						offset = logLike;
						like = 1.;
					} else {
						if(logLike>offset) { /*compute max in offset just to avoid numerical precision issues*/
							like *= exp(offset-logLike);
							offset = logLike;
						}
						like += exp(logLike-offset);
					}
				}
			}
			if(like>0.)
				return log(like)+offset;
			else
				return NEG_INFTY;
		case UncertaintyTypeC:
			for(l=0; l<data[n].dataStd.BCD.size; l++) {
				double logLike;
				logLike = getLogProbTypeC(t, data[n].dataStd.BCD.stopTime, tree->maxTime, l+data[n].dataStd.BCD.min, data[n].dataStd.BCD.nLeafMin, param)+data[n].dataStd.BCD.like[l]-logFactorial(l+data[n].dataStd.BCD.min-1);
				if(isinf(logLike) || isnan(logLike)) {
					like += 0.;
				} else {
					if(offset == NEG_INFTY) {
						offset = logLike;
						like = 1.;
					} else {
						if(logLike>offset) { /*compute max in offset just to avoid numerical precision issues*/
							like *= exp(offset-logLike);
							offset = logLike;
						}
						like += exp(logLike-offset);
					}
				}
			}
			if(like>0.)
				return log(like)+offset;
			else
				return NEG_INFTY;
		case UncertaintyTypeD:
			for(l=0; l<data[n].dataStd.BCD.size; l++) {
				double logLike;
				logLike = getLogProbTypeD(t, data[n].dataStd.BCD.stopTime, tree->maxTime, l+data[n].dataStd.BCD.min, param)+data[n].dataStd.BCD.like[l]-logFactorial(l+data[n].dataStd.BCD.min-1);
				if(isinf(logLike) || isnan(logLike)) {
					like += 0.;
				} else {
					if(offset == NEG_INFTY) {
						offset = logLike;
						like = 1.;
					} else {
						if(logLike>offset) { /*compute max in offset just to avoid numerical precision issues*/
							like *= exp(offset-logLike);
							offset = logLike;
						}
						like += exp(logLike-offset);
					}
				}
			}
			if(like>0.)
				return log(like)+offset;
			else
				return NEG_INFTY;
		case UncertaintyTypeA:
			return getLogProbTypeA(t, tree->maxTime, data[n].dataStd.A.numberTaxa, param)+data[n].dataStd.A.logRank-logFactorial(data[n].dataStd.A.numberTaxa-1);
		case UncertaintyTypeBbase:
			return getLogProbTypeB(t, data[n].dataStd.BCD.stopTime, tree->maxTime, 0, param);
		case UncertaintyTypeCbase:
			return getLogProbTypeC(t, data[n].dataStd.BCD.stopTime, tree->maxTime, 1, 1, param);
		case UncertaintyTypeDbase:
			return getLogProbTypeD(t, data[n].dataStd.BCD.stopTime, tree->maxTime, 1, param);
		default:
			error("\nExecution error in function getLogLikelihoodSubtree: No type while computing the likelihood of %d\n", n);
	}
	return 0.0;
}

/*return p(k, startTime)*/
double getLogIntProbTypeA(double a, double b, double maxTime, int k,  TypeUncertaintyLikeParameters *param) {
	return
		((double)k)*log(param->sampling)
		-log(param->birth)
		-log((double)k)
		+log(pow((1.-exp(param->omega*(maxTime-a)))/(param->bs-param->as*exp(param->omega*(maxTime-a))), (double)k)-pow((1.-exp(param->omega*(maxTime-b)))/(param->bs-param->as*exp(param->omega*(maxTime-b))), (double)k));
}

double getLogIntProbTypeC(double a, double b, double endTime, double maxTime, int k, int l, TypeUncertaintyLikeParameters *param) {
	return
		((double)k-l)*myLog(param->bs,-param->as*exp(param->omega*(maxTime-endTime)))
		+((double)l)*myLog((1-param->alpha)*param->bs,-(1-param->beta)*param->as*exp(param->omega*(maxTime-endTime)))
		-log((double)k)
		-log(param->birth)
		+log(pow((1-exp(param->omega*(endTime-a)))/(param->bs-param->as*exp(param->omega*(maxTime-a))),k)-pow((1-exp(param->omega*(endTime-b)))/(param->bs-param->as*exp(param->omega*(maxTime-b))),k));
}


double getLogIntProbTypeD(double a, double b, double endTime, double maxTime, int k,  TypeUncertaintyLikeParameters *param) {
	return
		((double)k)*(log(param->bs-param->as*exp(param->omega*(maxTime-endTime)))-log(param->bma))
		-log((double)k)
		-log(param->birth)
		+log(pow((1-exp(param->omega*(endTime-a)))/(param->bs-param->as*exp(param->omega*(maxTime-a))),k)-pow((1-exp(param->omega*(endTime-b)))/(param->bs-param->as*exp(param->omega*(maxTime-b))),k));
}

/*return p(k, startTime)*/
double getLogIntProbTypeB(double a, double b, double endTime, double maxTime, int k,  TypeUncertaintyLikeParameters *param) {
	return
		+log(param->fossil)
		-log(param->birth)
		-log((double)k+1.)
		+((double)k+1.)*(log(param->bs-param->as*exp(param->omega*(maxTime-endTime)))-log(param->bma))
		+log(pow((1.-exp(param->omega*(endTime-a)))/(param->bs-param->as*exp(param->omega*(maxTime-a))),(double)k+1.)
			-pow((1.-exp(param->omega*(endTime-b)))/(param->bs-param->as*exp(param->omega*(maxTime-b))),(double)k+1.));
}

/*return p(k, startTime)*/
double getLogProbTypeA(double startTime, double maxTime, int k,  TypeUncertaintyLikeParameters *param) {
    double time = maxTime-startTime;
	return
		((double)k)*log(param->sampling)
		+2.*log(param->bma)
		+param->omega*time
		+((double)k-1)*log(1.-exp(param->omega*time))
		-((double)k+1)*log(param->bs-param->as*exp(param->omega*time));
}

/* return P_x(k, startTime, endTime)/(P^\star_o)^k*/
double getLogProbTypeB(double startTime, double endTime, double maxTime, int k,  TypeUncertaintyLikeParameters *param) {
	if(startTime == endTime) {
		if(k == 1)
			return log(param->fossil);
		else
			return NEG_INFTY;
	} else
		return        
			log(param->fossil)
			+param->omega*(endTime-startTime)
			+((double)k)*(log(1-exp(param->omega*(endTime-startTime)))-log(param->bma))
			+((double)k+2.)*(log(param->bs-param->as*exp(param->omega*(maxTime-endTime)))-log(param->bs-param->as*exp(param->omega*(maxTime-startTime))));
}

/* return P_y(k, startTime, endTime)/(P^\star_o)^(k-l)*/
double getLogProbTypeC(double startTime, double endTime, double maxTime, int k, int l, TypeUncertaintyLikeParameters *param) {
    return
        param->omega*(endTime-startTime)
        +((double)k-1.)*(myLog(1.,-exp(param->omega*(endTime-startTime)))-log(param->bma))
        -((double)k+1.)*myLog(param->bs,-param->as*exp(param->omega*(maxTime-startTime)))
        +((double)k-l+1.)*myLog(param->bs,-param->as*exp(param->omega*(maxTime-endTime)))
       +((double)l)*myLog((1-param->alpha)*param->bs,-(1-param->beta)*param->as*exp(param->omega*(maxTime-endTime)));
}

/* return P_y(k, startTime, endTime)/(P^\star_o)^(k)*/
double getLogProbTypeD(double startTime, double endTime, double maxTime, int k, TypeUncertaintyLikeParameters *param) {
    return
        param->omega*(endTime-startTime)
        +((double)k-1.)*(myLog(1.,-exp(param->omega*(endTime-startTime)))-log(param->bma))
        +((double)k+1.)*(myLog(param->bs,-param->as*exp(param->omega*(maxTime-endTime)))-myLog(param->bs,-param->as*exp(param->omega*(maxTime-startTime))));
}

double getLogProbNotObservableU(double t, double maxTime, TypeUncertaintyLikeParameters *param) {
    return log(param->alpha)+log(param->beta)+log((1.-exp(param->omega*(maxTime-t))))-log(param->beta-param->alpha*exp(param->omega*(maxTime-t)));
}


#define PRECISION_LOG 0.000001
//return log(a+b) by approximating by log(max(a,b))+min(a,b)/max(a,b)
double myLog(double a,double b) {
	double max, min, rat;
	if(a>b) {
		max = a;
		min = b;
	} else {
		min = a;
		max = b;
	}
	rat = min/max;
	if(fabs(min/max)<PRECISION_LOG)
		return log(max)+rat;
	else
		return log(a+b);
}


/*For all nodes n of tree, indexMin[n] = leaf with the smallest time in subtree n or descendent diff from n with the smallest time*/
void fillIndexMin(int n, TypeTree *tree, int *indexMin) {
	if(tree->node[n].child == NOSUCH) {
		indexMin[n] = n;
	} else {
		int c;
		for(c=tree->node[n].child; c != NOSUCH; c=tree->node[c].sibling)
			fillIndexMin(c, tree, indexMin);
		if(tree->time[n] != NO_TIME)
			indexMin[n] = n;
		else {
			indexMin[n] = indexMin[tree->node[n].child];
			for(c=tree->node[tree->node[n].child].sibling; c != NOSUCH; c=tree->node[c].sibling)
				if(tree->time[indexMin[c]]<tree->time[indexMin[n]])
					indexMin[n] = indexMin[c];
				else
					if(tree->time[indexMin[c]]==tree->time[indexMin[n]] && tree->node[indexMin[c]].child != NOSUCH)
						indexMin[n] = indexMin[c];
		}
	}
}

void fillTableTreeIter(int n, TypeTree *tree, int *indexMin, double *K, double *T) {
    if(tree->node[n].child != NOSUCH) {
        int child[2];
        child[0] = tree->node[n].child;
        child[1] = tree->node[child[0]].sibling;
        if(child[1] == NOSUCH) {
            error("Execution error: node %d with a single child\n", n);
        }
        if(tree->node[child[1]].sibling != NOSUCH) {
            error("Execution error: node %d with more than 2 children\n", n);
        }
        fillTableTreeIter(child[0], tree, indexMin, K, T);
        fillTableTreeIter(child[1], tree, indexMin, K, T);
        if(tree->time[indexMin[n]]<tree->maxTime) {
            double size[2];
            if(tree->time[indexMin[n]] == tree->time[indexMin[child[0]]])
                size[0] = K[child[0]];
            else
                size[0] = T[child[0]];
            if(tree->time[indexMin[n]] == tree->time[indexMin[child[1]]])
                size[1] = K[child[1]];
            else
                size[1] = T[child[1]];
            K[n] = size[0]*size[1];
        } else
            K[n] = 1;
        T[n] = T[child[0]]*T[child[1]]+1.;
    } else {
        if(tree->time[indexMin[n]]<tree->maxTime)
            K[n] = 1.;
        else
            K[n] = 0.;
        T[n] = 1.;
    }
}

double getItemNumber(TypeTree *tree, TypeFossilFeature *fos) {
    TypeTree **treeList;
    int l, size;
    double res = 0.;
    splitTreeFossil(tree, fos, &treeList, &size);
    for(l=0; l<size; l++) {
        res += getItemNumberSplitted(treeList[l]);
        if(treeList[l]->info!=NULL)
            free((void*)treeList[l]->info);
        freeTree(treeList[l]);
    }
    free((void*)treeList);
   return res;
}

double getItemNumberSplitted(TypeTree *tree) {
    int  *indexMin, n;
    double sum = 0., *K, *T;
    if(tree == NULL || tree->size == 0)
            return 1.;
    indexMin = (int*) malloc(tree->size*sizeof(int));
    fillIndexMin(tree->root, tree, indexMin);
    K = (double*) malloc(tree->size*sizeof(double));
    T = (double*) malloc(tree->size*sizeof(double));
    fillTableTreeIter(tree->root, tree, indexMin, K, T);
    for(n=0; n<tree->size; n++)
        sum += K[n]+T[n];
    free((void*)indexMin);
    free((void*)K);
    free((void*)T);
    return sum;
}

TypeUncertaintyLikeParameters getUncertaintyLikeParameters(TypeModelParam *param) {
	TypeUncertaintyLikeParameters res;
	res.birth = param->birth;
	res.death = param->death;
	res.fossil = param->fossil;
	res.sampling = param->sampling;
	res.alpha = (res.birth+res.death+res.fossil-sqrt(pow(res.birth+res.death+res.fossil, 2.)-4.*res.birth*res.death))/(2*res.birth);
	res.beta = (res.birth+res.death+res.fossil+sqrt(pow(res.birth+res.death+res.fossil, 2.)-4.*res.birth*res.death))/(2*res.birth);
	res.bma = res.beta-res.alpha;
	res.omega = -res.birth*res.bma;
	res.ab = res.alpha*res.beta;
	res.bmab = res.beta-res.ab;
	res.amab = res.alpha-res.ab;
	res.bs = res.beta-1+res.sampling;
	res.as = res.alpha-1+res.sampling;
	return res;
}


int fillBasicTree(TypeTree *basic, int n, TypeTree *tree, TypeFossilTab *ftab, int toCompute, int *newIndex) {
	if(n == toCompute)
		*newIndex = basic->size;
	int curInd = basic->size++;
	if(ftab[n].size>0) {
		basic->time[curInd] = ftab[n].time[0];
		basic->node[curInd].child = NOSUCH;
		((TypeNodeStatus*)basic->info)[curInd]=extinctNodeStatus;
	} else {
		basic->time[curInd] = tree->time[n];
		if(tree->node[n].child != NOSUCH) {
			if(tree->time[n] != NO_TIME)
				((TypeNodeStatus*)basic->info)[curInd] = divergenceNodeStatus;
			else
				((TypeNodeStatus*)basic->info)[curInd] = noneNodeStatus;
			basic->node[curInd].child = fillBasicTree(basic, tree->node[n].child, tree, ftab, toCompute, newIndex);
			basic->node[basic->node[curInd].child].sibling = fillBasicTree(basic, tree->node[tree->node[n].child].sibling, tree, ftab, toCompute, newIndex);
			basic->node[basic->node[basic->node[curInd].child].sibling].sibling = NOSUCH;
		} else {
			if(tree->time[n] == tree->maxTime)
				((TypeNodeStatus*)basic->info)[curInd]=contempNodeStatus;
			else
				((TypeNodeStatus*)basic->info)[curInd]=unknownNodeStatus;
			basic->node[curInd].child = NOSUCH;
		}
	}
	return curInd;
}

/*Get the basic tree containing the divergence n*/
TypeTree *getBasicTreeDivergence(int n, TypeTree *tree, TypeFossilTab *ftab, int *newIndex, int *root) {
    TypeTree *basic;
    if(tree->parent == NULL)
        setParent(tree);
    for(*root=n; *root!=tree->root && ftab[*root].size==0; *root=tree->parent[*root])
        ;
    basic = newTree(tree->size);
    basic->root = 0;
    if(ftab[*root].size>0) {
        basic->minTimeInt.inf = ftab[*root].time[ftab[*root].size-1];
		basic->minTimeInt.sup = ftab[*root].time[ftab[*root].size-1];
    } else {
        basic->minTimeInt.inf = tree->minTimeInt.inf;
		basic->minTimeInt.sup = tree->minTimeInt.sup;
	}
    basic->maxTime = tree->maxTime;
    basic->info = (TypeNodeStatus*) malloc(tree->size*sizeof(TypeNodeStatus));
    ftab[*root].size = 0;
    fillBasicTree(basic, *root, tree, ftab, n, newIndex);
    reallocTree(basic->size, basic);
    basic->info = (TypeNodeStatus*) realloc((void*)basic->info, basic->size*sizeof(TypeNodeStatus));
    return basic;
}

double getMinTimeTree(int n, TypeTree *tree, TypeFossilTab *ftab) {
	if(ftab!= NULL && ftab[n].size>0)
		return ftab[n].time[0];
	if(tree->time[n] != NO_TIME)
		return tree->time[n];
	if(tree->node[n].child == NOSUCH) {
		error("Execution error leaf %d (%s) with no time\n", n, (tree->name && tree->name[n])?tree->name[n]:"??");
	}
	return utils_MIN(getMinTimeTree(tree->node[n].child, tree, ftab), getMinTimeTree(tree->node[tree->node[n].child].sibling, tree, ftab));
}

int fillSplitTreeFossil(TypeTree *tcur, int n, TypeTree *tree, TypeFossilTab *ftab, TypeTree **treeList, int *size) {
    int curInd = tcur->size++;
    if(ftab[n].size>0) {
        tcur->time[curInd] = ftab[n].time[0];
        tcur->node[curInd].child = -1;
        ((TypeNodeStatus*)tcur->info)[curInd]=extinctNodeStatus;
        if(ftab[n].size == 1 && tree->node[n].child == NOSUCH && (tree->time[n] == ftab[n].time[0] || tree->time[n]==NO_TIME)) { /*empty tree*/
            treeList[(*size)] = newTree(0);
            treeList[(*size)]->root = 0;
            treeList[(*size)]->minTime = ftab[n].time[0];
            treeList[(*size)]->minTimeInt.inf = ftab[n].time[0];
            treeList[(*size)]->minTimeInt.sup = ftab[n].time[0];
            treeList[(*size)]->maxTime = tree->maxTime;
            (*size)++;
            ftab[n].size--;
            ftab[n].time++;
        } else  { /*tree not empty*/
            treeList[(*size)] = newTree(tree->size);
            treeList[(*size)]->root = 0;
            treeList[(*size)]->minTime = tcur->time[curInd];
            treeList[(*size)]->minTimeInt.inf = tcur->time[curInd];
            treeList[(*size)]->minTimeInt.sup = tcur->time[curInd];
            treeList[(*size)]->maxTime = tree->maxTime;
            treeList[(*size)]->info = (TypeNodeStatus*) malloc(tree->size*sizeof(TypeNodeStatus));
            ((TypeNodeStatus*)treeList[(*size)]->info)[0]=noneNodeStatus;
            (*size)++;
            ftab[n].size--;
            ftab[n].time++;
            fillSplitTreeFossil(treeList[(*size)-1], n, tree, ftab, treeList, size);
        }
    } else {
        tcur->time[curInd] = tree->time[n];
        if(tree->node[n].child!=NOSUCH) {
            int c, prec;
			if(tree->time[n] == NO_TIME)
				((TypeNodeStatus*)tcur->info)[curInd] = noneNodeStatus;
			else
				((TypeNodeStatus*)tcur->info)[curInd] = divergenceNodeStatus;
            tcur->node[curInd].child = fillSplitTreeFossil(tcur, tree->node[n].child, tree, ftab, treeList, size);
            prec = tcur->node[curInd].child;
            for(c=tree->node[tree->node[n].child].sibling; c!=NOSUCH; c=tree->node[c].sibling) {
                tcur->node[prec].sibling = fillSplitTreeFossil(tcur, c, tree, ftab, treeList, size);
                prec = tcur->node[prec].sibling;
            }
            tcur->node[prec].sibling = -1;
        } else {
            if(tree->time[n]==tree->maxTime)
                ((TypeNodeStatus*)tcur->info)[curInd]=contempNodeStatus;
            else
                ((TypeNodeStatus*)tcur->info)[curInd]=unknownNodeStatus;
            tcur->node[curInd].child = -1;
        }
    }
    return curInd;
}

/*split the tree "tree" each time a fossill occurs*/
void splitTreeFossil(TypeTree *tree, TypeFossilFeature *fos, TypeTree ***treeList, int *size) {
    int i;
    TypeFossilTab *ftab, *ftmp;
    ftab = listToFossilTab(fos, tree->size);
    ftmp = (TypeFossilTab*) malloc(tree->size*sizeof(TypeFossilTab));
    for(i=0; i<tree->size; i++) {
        ftmp[i].size = ftab[i].size;
        ftmp[i].time = ftab[i].time;
    }
    *treeList = (TypeTree **) malloc((fos->size+1)*sizeof(TypeTree*));
    (*treeList)[0] = newTree(tree->size);
    (*treeList)[0]->root = 0;
    (*treeList)[0]->node[0].sibling = NOSUCH;
    (*treeList)[0]->minTime = tree->minTime;
    (*treeList)[0]->minTimeInt.inf = tree->minTimeInt.inf;
    (*treeList)[0]->minTimeInt.sup = tree->minTimeInt.sup;
    (*treeList)[0]->maxTime = tree->maxTime;
    (*treeList)[0]->info = (void*) malloc(tree->size*sizeof(TypeNodeStatus));
    if(tree->time[tree->root] == NO_TIME)
		((TypeNodeStatus*)(*treeList)[0]->info)[0] = noneNodeStatus;
	else
 		((TypeNodeStatus*)(*treeList)[0]->info)[0] = divergenceNodeStatus;
//    ((TypeNodeStatus*)(*treeList)[0]->info)[0]=noneNodeStatus;
   *size = 1;
    fillSplitTreeFossil((*treeList)[0], tree->root, tree, ftmp, *treeList, size);
    for(i=0; i<tree->size; i++)
        if(ftab[i].time != NULL)
            free((void*)ftab[i].time);
    free((void*)ftab);
    free((void*)ftmp);
    for(i=0; i<*size; i++)
        reallocTree((*treeList)[i]->size, (*treeList)[i]);
}

void freeUncertaintyLikeData(TypeUncertaintyLikeData* data) {
	if((data->type == UncertaintyTypeB || data->type == UncertaintyTypeC || data->type == UncertaintyTypeD || data->type == UncertaintyTypeBbase || data->type == UncertaintyTypeCbase || data->type == UncertaintyTypeDbase) &&  (data->dataStd.BCD.like != NULL))
		free((void*)data->dataStd.BCD.like);
}

