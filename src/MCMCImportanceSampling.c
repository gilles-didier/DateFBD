#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <pthread.h>
#include <signal.h>

#include "Utils.h"
#include "Uncertainty.h"
#include "MCMCImportanceSampling.h"


typedef struct MINIMIZATION_PARAM_DATA {
    TypeTree **tree;
    TypeFossilFeature **fos;
    int nTree;
} TypeMinimizationParamData;




typedef struct THREAD_PARAMETER_DENS {
	int *number;
	pthread_mutex_t *mutex_number, *mutex_sum;
	pthread_cond_t *cond_number;
	TypeTree *tree;
	TypeFossilFeature *ffe;
	int node;
	double *sum;
	TypeModelParam *param;
} TypeThreadParameterDens;

typedef struct THREAD_PARAMETER_DIST {
	int *number;
	pthread_mutex_t *mutex_number;
	pthread_cond_t *cond_number;
	TypeTree *tree;
	TypeFossilFeature *ffe;
	int node;
	TypeDistribution *logD;
	double *logCond;
	TypeModelParam *param;
} TypeThreadParameterDist;

typedef struct THREAD_PARAMETER_SINGLE {
	int *number;
	pthread_mutex_t *mutex_number;
	pthread_cond_t *cond_number;
	TypeTree *tree;
	TypeFossilFeature *ffe;
	int taxa;
	TypeDistribution *logD;
	TypeExtendedModelParam *pext;
} TypeThreadParameterExtSingleTree;

typedef struct THREAD_PARAMETER_EXT_MULTIPLE {
	int *number;
	pthread_mutex_t *mutex_number;
	pthread_cond_t *cond_number;
	int *clade;
	TypeTree *tree;
	TypeFossilFeature *ffe;
	TypeDistribution *logD;
	TypeModelParam *param;
	TypeExtendedModelParam *pext;
} TypeThreadParameterExtMultipleTree;


typedef struct THREAD_PARAMETER_SINGLE_QUANT {
	int *number;
	pthread_mutex_t *mutex_number;
	pthread_cond_t *cond_number;
	TypeTree *tree;
	TypeFossilFeature *ffe;
	int taxa;
	double order, *quant;
	TypeExtendedModelParam *pext;
} TypeThreadParameterExtQuantSingleTree;


static pthread_mutex_t mutexN = PTHREAD_MUTEX_INITIALIZER, mutexS =  PTHREAD_MUTEX_INITIALIZER;
static pthread_cond_t condN = PTHREAD_COND_INITIALIZER;
static int maxT = 40;
//static int maxT = 1;

static void *threadComputeExtSingleTree(void *data);
static void *threadComputeDens(void *data);
static void *threadComputeDistribution(void *data);
static void *threadComputeExtMultipleTree(void *data);
static void *threadComputeExtQuantSingleTree(void *data);

static void fillFolNode(TypeFossilFeature *f, int size, int *fol);
static void fillBranchNode(TypeFossilFeature *f, int size, int *branch);
static void  getProposalFossil(TypeFossilFeature *fp, TypeFossilIntFeature *fi, int *fol, double al, int *move, double *newTime);
static void changeFossil(TypeTree **tree, int nTree, TypeFossilFeature *fp, int *fol, int *branch, int move, double newTime);
static double update(TypeTree **tree, int nTree, TypeFossilFeature *fp, TypeFossilIntFeature *fi, int *fol, int *branch, double prob, double al, double propParam, TypeModelParam *param, TypeModelParam *windSize, double probSpe, double probExt);
static TypeFossilFeature *sampleFossil(TypeFossilIntFeature* feat, int size);
static int compareLocal(const void* a, const void* b);
static double getLogDensitySum(TypeTree **tree, int nTree, TypeFossilFeature *fp, TypeModelParam *param);
static double getLogDensitySumNoThread(TypeTree **tree, int nTree, TypeFossilFeature *fp, TypeModelParam *param);
static double updateFossilFixed(TypeTree **tree, int nTree, TypeFossilFeature *fp, double prob, TypeModelParam *param, TypeModelParam *windSize, double probSpe, double probExt);
double getLogProbExtinctCondSum(double t, double startTime, double maxTime, TypeExtendedModelParam *param, int size);
double getQuantileSum(double q, double tol, double startTime, double maxTime, TypeExtendedModelParam *param, int size);
double getLogProbExtinctCondCladeSum(double val, int *clade, double maxTime, TypeExtendedModelParam *param, double **time, int size);
double getQuantileCladeSum(int *clade, double q, double tol, double maxTime, TypeExtendedModelParam *param, double **time, int size);
static double getLogProbComp(int *cladeA, int *cladeB, double maxTime, TypeExtendedModelParam *param, double *time);
static double getProbComp(int *cladeA, int *cladeB, double *time, double tmax, TypeExtendedModelParam *pext);

static TypeFossilList *fossilList;

/*return log(a+b) in an accurate way*/
double logSum(double a,double b) {
	double max, min;
	if(a>b) {
		max = a;
		min = b;
	} else {
		min = a;
		max = b;
	}
	return log(max) + log1p(min/max);
}

/*return log(exp(a)+exp(b)) in an accurate way*/
double logSumLog(double a, double b) {
	double max, min;
	if(a == NEG_INFTY)
		return b;
	if(b == NEG_INFTY)
		return a;
	if(a>b) {
		max = a;
		min = b;
	} else {
		min = a;
		max = b;
	}
	return max + log1p(exp(min-max));
}


int compareLocal(const void* a, const void* b) {
    if(fossilList[*((int*)a)].time>fossilList[*((int*)b)].time)
        return 1;
    if(fossilList[*((int*)a)].time<fossilList[*((int*)b)].time)
        return -1;
    return 0;
}

TypeFossilFeature *sampleFossil(TypeFossilIntFeature* feat, int size) {
    TypeFossilFeature *sample;
    int i, n, f, *tmp;

    if(feat == NULL)
        return NULL;
    sample = (TypeFossilFeature*) malloc(sizeof(TypeFossilFeature));
    sample->size = feat->sizeFossil;
    sample->sizeBuf = feat->sizeFossil;
    sample->fossilList = (TypeFossilList*) malloc(sample->size*sizeof(TypeFossilList));
    sample->fossil = (int*) malloc(size*sizeof(int));
    sample->status = (TypeNodeStatus*) malloc(size*sizeof(int));
    for(i=0; i<size; i++)
        sample->fossil[i] = feat->fossilInt[i];
    for(i=0; i<sample->size; i++) {
        sample->fossilList[i].prec = feat->fossilIntList[i].prec;
        sample->fossilList[i].time = feat->fossilIntList[i].fossilInt.inf+UNIF_RAND*(feat->fossilIntList[i].fossilInt.sup-feat->fossilIntList[i].fossilInt.inf);
    }
    fossilList = sample->fossilList;
    tmp = (int*) malloc((sample->size+1)*sizeof(int));
    for(n=0; n<size; n++) {
        int ind = 0;
        for(f=sample->fossil[n]; f!=NOSUCH; f=sample->fossilList[f].prec)
            tmp[ind++] = f;
        if(ind>0) {
            qsort(tmp, ind, sizeof(int), compareLocal);
            sample->fossilList[tmp[0]].prec = NOSUCH;
            for(f=1; f<ind; f++)
                sample->fossilList[tmp[f]].prec = tmp[f-1];
            sample->fossil[n] = tmp[ind-1];
        }
    }
    free((void*) tmp);
    return sample;
}


void fillFolNode(TypeFossilFeature *f, int size, int *fol) {
	int n;
	for(n=0; n<size; n++) {
		int k, tmp;
		tmp = NOSUCH;
		for(k=f->fossil[n]; k!=NOSUCH; k=f->fossilList[k].prec) {
			fol[k] = tmp;
			tmp = k;
		}
	}
}

void fillBranchNode(TypeFossilFeature *f, int size, int *branch) {
	int n;
	for(n=0; n<size; n++) {
		int k;
		for(k=f->fossil[n]; k!=NOSUCH; k=f->fossilList[k].prec)
			branch[k] = n;
	}
}

TypeModelParam getProposalParam(TypeModelParam *param, TypeModelParam *windSize, double probSpe, double probExt) {
	TypeModelParam res;
	double p = UNIF_RAND;
	res = *param;
	if(p<probSpe) {
		res.birth = UNIF_RAND*windSize->birth+param->birth-0.5*windSize->birth;
			if(res.birth < 0.)
				res.birth = -res.birth;
	} else {
		if(p<probSpe+probExt) {
			res.death = UNIF_RAND*windSize->death+param->death-0.5*windSize->death;
			if(res.death < 0.)
				res.death = -res.death;
		} else {
			res.fossil = UNIF_RAND*windSize->fossil+param->fossil-0.5*windSize->fossil;
			if(res.fossil < 0.)
				res.fossil = -res.fossil;
		}
	}
	return res;
}

void getProposalFossil(TypeFossilFeature *fp, TypeFossilIntFeature *fi, int *fol, double al, int *move, double *newTime) {
	double length;
    *move = RANGE_RAND(fp->size);
    length = (fi->fossilIntList[*move].fossilInt.sup-fi->fossilIntList[*move].fossilInt.inf)*al/2.;
    *newTime = UNIF_RAND*(2.*length)+fp->fossilList[*move].time-length;
    if(*newTime<fi->fossilIntList[*move].fossilInt.inf)
		*newTime = 2.*fi->fossilIntList[*move].fossilInt.inf - *newTime;
    if(*newTime>fi->fossilIntList[*move].fossilInt.sup)
		*newTime = 2.*fi->fossilIntList[*move].fossilInt.sup - *newTime;
}

void changeFossil(TypeTree **tree, int nTree, TypeFossilFeature *fp, int *fol, int *branch, int move, double newTime) {
	int xA, xB, change, i;
	change = (tree[0]->time[branch[move]] == fp->fossilList[fp->fossil[branch[move]]].time);
	fp->fossilList[move].time = newTime;
	xA = NOSUCH;
	for(xB=fol[move]; xB!=NOSUCH && fp->fossilList[move].time>fp->fossilList[xB].time; xB=fol[xB])
		xA = xB;
	if(xB != fol[move]) {;
		if(fol[move] != NOSUCH)
			fp->fossilList[fol[move]].prec = fp->fossilList[move].prec;
		else
			fp->fossil[branch[move]] = fp->fossilList[move].prec;
		if(fp->fossilList[move].prec != NOSUCH)
			fol[fp->fossilList[move].prec] = fol[move];
		if(xB != NOSUCH)
			fp->fossilList[xB].prec = move;
		else
			fp->fossil[branch[move]] = move;
		if(xA != NOSUCH)
			fol[xA] = move;
		fol[move] = xB;
		fp->fossilList[move].prec = xA;
	}
	xA = NOSUCH;
	for(xB=fp->fossilList[move].prec; xB!=NOSUCH && fp->fossilList[move].time<fp->fossilList[xB].time; xB=fp->fossilList[xB].prec)
		xA = xB;
	if(xB != fp->fossilList[move].prec) {
		if(fol[move] != NOSUCH)
			fp->fossilList[fol[move]].prec = fp->fossilList[move].prec;
		else
			fp->fossil[branch[move]] = fp->fossilList[move].prec;
		if(fp->fossilList[move].prec != NOSUCH)
			fol[fp->fossilList[move].prec] = fol[move];
		if(xB != NOSUCH)
			fol[xB] = move;
		if(xA != NOSUCH)
			fp->fossilList[xA].prec = move;
		fol[move] = xA;
		fp->fossilList[move].prec = xB;
	}
	double min = getMinFossilTime(fp);
	if(tree[0]->minTimeInt.sup > min)
		for(i=0; i<nTree; i++)
			tree[i]->minTimeInt.sup = min;
	if(change)
		for(i=0; i<nTree; i++)
			tree[i]->time[branch[move]] = fp->fossilList[fp->fossil[branch[move]]].time;
}

double update(TypeTree **tree, int nTree, TypeFossilFeature *fp, TypeFossilIntFeature *fi, int *fol, int *branch, double prob, double al, double propParam, TypeModelParam *param, TypeModelParam *windSize, double probSpe, double probExt) {
	int move;
	double newTime, oldTime, newProb;
	if(UNIF_RAND < propParam) {
		TypeModelParam newParam = getProposalParam(param, windSize, probSpe, probExt);
		newProb = getLogDensitySum(tree, nTree, fp, &newParam);
		if(isnan(newProb))
			error("Param %.2le %.2le %.2le\n", newParam.birth, newParam.death, newParam.fossil);
		if(UNIF_RAND < exp(newProb-prob)) {
			*param = newParam;
			return newProb;
		} else {
			return prob;
		}
	} else { 	
		getProposalFossil(fp, fi, fol, al, &move, &newTime);
		oldTime = fp->fossilList[move].time;
		changeFossil(tree, nTree, fp, fol, branch, move, newTime);
		newProb = getLogDensitySum(tree, nTree, fp, param);
		if(isnan(newProb))
			error("fossil %d %.5lf -> %.5lf (%d %s) time %.5lf\n", move, oldTime, newTime, branch[move], tree[0]->name[branch[move]], tree[0]->time[branch[move]]);
		if(UNIF_RAND < exp(newProb-prob)) {
			return newProb;
		} else {
			changeFossil(tree, nTree, fp, fol, branch, move, oldTime);
			return prob;
		}
	}
	return prob;
}

double updateFossilFixed(TypeTree **tree, int nTree, TypeFossilFeature *fp, double prob, TypeModelParam *param, TypeModelParam *windSize, double probSpe, double probExt) {
	double newProb;
	TypeModelParam newParam = getProposalParam(param, windSize, probSpe, probExt);
	newProb = getLogDensitySumNoThread(tree, nTree, fp, &newParam);
	if(isnan(newProb))
		error("Param %.2le %.2le %.2le\n", newParam.birth, newParam.death, newParam.fossil);
	if(UNIF_RAND < exp(newProb-prob)) {
		*param = newParam;
		return newProb;
	} else {
		return prob;
	}
	return prob;
}
double getLogDensitySumNoThread(TypeTree **tree, int nTree, TypeFossilFeature *fp, TypeModelParam *param) {
	double sum;
	int i;
	sum = NEG_INFTY;
	for(i=0; i<nTree; i++)
		sum = logSumLog(sum, getLogLikelihoodTreeFossil(tree[i], fp, param));
	return sum;
}

double getLogDensitySum(TypeTree **tree, int nTree, TypeFossilFeature *fp, TypeModelParam *param) {
	double sum;
	int cont=1, nT=0, i = 0;
	cont = 1;
	sum = NEG_INFTY;
	i = 0;
	while(cont) {
		pthread_mutex_lock(&mutexN);
		while(i < nTree && nT < maxT) {
			pthread_t thread;	
			int ret = 0;
			TypeThreadParameterDens *tp;
			tp = (TypeThreadParameterDens*) malloc(sizeof(TypeThreadParameterDens));
			tp->tree = tree[i];
			tp->ffe = fp;
			tp->param = param;
			tp->sum = &sum;
			tp->mutex_sum = &mutexS;
			tp->mutex_number = &mutexN;
			tp->cond_number = &condN;
			tp->number = &nT;
			if((ret = pthread_create(&thread, NULL, threadComputeDens, (void*) tp)) == 0) {
				int err;
				if((err = pthread_detach(thread)) == 0) {
					nT++; i++;
				} else {
					warning("Error %d while detaching thread: %s\n", err, (char*) strerror(err));
				}
			} else
				warning("Error %d while creating thread: %s\n", ret, (char*) strerror(ret));
		}
		cont = (nT > 0);
		if(cont)
			pthread_cond_wait(&condN, &mutexN);
		pthread_mutex_unlock(&mutexN);
	}
	return sum;
}

void *threadComputeDens(void *data) {
	double tmp;
	tmp = getLogLikelihoodTreeFossil(((TypeThreadParameterDens*)data)->tree, ((TypeThreadParameterDens*)data)->ffe, ((TypeThreadParameterDens*)data)->param);
	pthread_mutex_lock(((TypeThreadParameterDens*)data)->mutex_sum);
		*(((TypeThreadParameterDens*)data)->sum) = logSumLog(*(((TypeThreadParameterDens*)data)->sum), tmp);
	pthread_mutex_unlock(((TypeThreadParameterDens*)data)->mutex_sum);
	pthread_mutex_lock(((TypeThreadParameterDens*)data)->mutex_number);
		(*((TypeThreadParameterDens*)data)->number)--;
		pthread_cond_signal(((TypeThreadParameterDens*)data)->cond_number);
	pthread_mutex_unlock(((TypeThreadParameterDens*)data)->mutex_number);
	free(data);
	return NULL;
}

TypeDistribution MCMCSamplingDist(FILE *fout, FILE *find, int *node, TypeTree **tree, TypeFossilIntFeature *fi, double step, double al, int burn, int gap, int iter, double prop, TypeModelParam *windSize, TypeModelParam *init, double probSpe, double probExt) {
	int *fol, *branch, i, j, s, n, nTree;
	double prob, min, max, *logCond, logSumCond;
	TypeFossilFeature *fp;
	TypeDistribution *logDTmp, logSTmp, dist;
	TypeModelParam param, *saveParam;
	
	param.birth = UNIF_RAND*init->birth;
	param.death = UNIF_RAND*init->death;
	param.fossil = UNIF_RAND*init->fossil;
	param.sampling = 1.;

	fp = sampleFossil(fi, tree[0]->size);
	for(nTree=0; tree[nTree] != NULL; nTree++) {
		for(n=0; n<tree[nTree]->size; n++) {
			if(tree[nTree]->node[n].child == NOSUCH) {
				switch(fi->status[n]) {
					case contempNodeStatus:
						tree[nTree]->time[n] = tree[nTree]->maxTime;
					break;
					case unknownNodeStatus:

					break;
					case extinctNodeStatus:
						tree[nTree]->time[n] = fp->fossilList[fp->fossil[n]].time;
					break;
					default:
	//						error("Node %d has no status\n", n);
	;
				}
			} else
				tree[nTree]->time[n] = NO_TIME;
		}
	}
	
	fol = (int*) malloc(fp->size*sizeof(int));
	fillFolNode(fp, tree[0]->size, fol);
	branch = (int*) malloc(fp->size*sizeof(int));
	fillBranchNode(fp, tree[0]->size, branch);
	min = getMaxFossilIntTimeToNode(node[0], tree[0], fi);
	max = getMinFossilIntTimeFromNode(node[0], tree[0], fi);
	min = floor(min/step)*step;
	max = (ceil(max/step))*step;
	double **save;
	int f;
	saveParam = (TypeModelParam*) malloc(iter*sizeof(TypeModelParam));
	save = (double**) malloc(fp->size*sizeof(double*));
	for(f=0; f<fp->size; f++)
		save[f] = (double*) malloc(iter*sizeof(double));
	prob = getLogDensitySum(tree, nTree, fp, &param);

	for(i=0; i<burn; i++) {
		prob = update(tree, nTree, fp, fi, fol, branch, prob, al, prop, &param, windSize, probSpe, probExt);
if(i%10 == 0) {
	fprintf(stderr, "it %-7d\r", i); fflush(stderr);
}
	}
	logSumCond = NEG_INFTY;
	logCond = (double*) malloc(nTree*sizeof(double));
	logSTmp.size = (int) ceil((max-min)/step)+1;
	logSTmp.item = (TypeDistributionItem*) malloc(logSTmp.size*sizeof(TypeDistributionItem));
	logDTmp = (TypeDistribution*) malloc(nTree*sizeof(TypeDistribution));
	for(i=0; i<nTree; i++) {
		logDTmp[i].size = logSTmp.size;
		logDTmp[i].item = (TypeDistributionItem*) malloc(logSTmp.size*sizeof(TypeDistributionItem));
	}
	dist.size = logSTmp.size;
	dist.item = (TypeDistributionItem*) malloc(dist.size*sizeof(TypeDistributionItem));
	for(j=0; j<logSTmp.size; j++) {
		for(i=0; i<nTree; i++)
			logDTmp[i].item[j].val = min + ((double)j)*step;
		logSTmp.item[j].val = min + ((double)j)*step;
		dist.item[j].val = min + ((double)j)*step;
		dist.item[j].dens = NEG_INFTY;
	}
	for(s=0; s<iter; s++) {
		int j;
		int cont=1, nT=0, i = 0;
		while(cont) {
			pthread_mutex_lock(&mutexN);
			while(i < nTree && nT < maxT) {
				pthread_t thread;	
				int ret = 0;
				TypeThreadParameterDist *tp;
				tp = (TypeThreadParameterDist*) malloc(sizeof(TypeThreadParameterDist));
				tp->tree = tree[i];
				tp->ffe = fp;
				tp->param = &param;
				tp->node = node[i];
				tp->logD = &(logDTmp[i]);
				tp->logCond = &(logCond[i]);
				tp->mutex_number = &mutexN;
				tp->cond_number = &condN;
				tp->number = &nT;
				if((ret = pthread_create(&thread, NULL, threadComputeDistribution, (void*) tp)) == 0) {
					int err;
					if((err = pthread_detach(thread)) == 0) {
						nT++; i++;
					} else {
						warning("Error %d while detaching thread: %s\n", err, (char*) strerror(err));
					}
				} else
					warning("Error %d while creating thread: %s\n", ret, (char*) strerror(ret));
			}
			cont = (nT > 0);
			if(cont)
				pthread_cond_wait(&condN, &mutexN);
			pthread_mutex_unlock(&mutexN);
		}
		logSumCond = NEG_INFTY;
		for(j=0; j<logSTmp.size; j++)
			logSTmp.item[j].dens = NEG_INFTY;
		for(i=0; i<nTree; i++) {
			logSumCond = logSumLog(logSumCond, logCond[i]);
			for(j=0; j<dist.size; j++)
				logSTmp.item[j].dens = logSumLog(logSTmp.item[j].dens, logDTmp[i].item[j].dens);
		}
		for(j=0; j<dist.size; j++)
			dist.item[j].dens = logSumLog(dist.item[j].dens, logSTmp.item[j].dens-logSumCond);
		saveParam[s] = param;
		for(f=0; f<fp->size; f++)
			save[f][s] = fp->fossilList[f].time;
		for(j=0; j<gap; j++) {
			prob = update(tree, nTree, fp, fi, fol, branch, prob, al, prop, &param, windSize, probSpe, probExt);
		}
if(s%10 == 0) {
	fprintf(stderr, "it %-7d\r", s); fflush(stderr);
}
	}
	for(j=0; j<dist.size; j++)
		dist.item[j].dens = exp(dist.item[j].dens-log((double)iter));
	for(s=0; s<iter; s++)
		fprintf(fout, "%d\t%lf\n",  s+burn+1, saveParam[s].birth);
	for(s=0; s<iter; s++)
		fprintf(fout, "%d\t%lf\n",  s+burn+1, saveParam[s].death);
	for(s=0; s<iter; s++)
		fprintf(fout, "%d\t%lf\n",  s+burn+1, saveParam[s].fossil);
	free((void*)saveParam);
	for(f=0; f<fp->size; f++)
		for(s=0; s<iter; s++)
			fprintf(fout, "%d\t%lf\n",  s+burn+1, save[f][s]);
	for(f=0; f<fp->size; f++)
		free((void*)save[f]);
	free((void*)save);
	int off = 0;
	fprintf(find, "birth\t%d\t%d\n", off*iter+1, (off+1)*iter); off++;
	fprintf(find, "death\t%d\t%d\n", off*iter+1, (off+1)*iter); off++;
	fprintf(find, "fossil\t%d\t%d\n", off*iter+1, (off+1)*iter); off++;
	for(f=0; f<fp->size; f++) {
		fprintf(find, "%d", f+1);
		if(branch[f] != NOSUCH && tree[0]->name[branch[f]] != NULL)
			fprintf(find, "_%s_%d_%d", tree[0]->name[branch[f]], (int) round(fi->fossilIntList[f].fossilInt.inf), (int) round(fi->fossilIntList[f].fossilInt.sup));
		fprintf(find, "\t%d\t%d\n", off*iter+1, (off+1)*iter); off++;
	}
	free((void*)branch);
	free((void*)fol);
	free((void*)logSTmp.item);
	for(i=0; i<nTree; i++)
		free((void*)logDTmp[i].item);
	free((void*)logDTmp);
	free((void*)logCond);
	return dist;
}


void *threadComputeDistribution(void *data) {
	fillLogDistribution(((TypeThreadParameterDist*)data)->logD, ((TypeThreadParameterDist*)data)->logCond, ((TypeThreadParameterDist*)data)->node, ((TypeThreadParameterDist*)data)->tree, ((TypeThreadParameterDist*)data)->ffe, ((TypeThreadParameterDist*)data)->param);
	pthread_mutex_lock(((TypeThreadParameterDist*)data)->mutex_number);
		(*((TypeThreadParameterDist*)data)->number)--;
		pthread_cond_signal(((TypeThreadParameterDist*)data)->cond_number);
	pthread_mutex_unlock(((TypeThreadParameterDist*)data)->mutex_number);
	free(data);
	return NULL;
}


void MCMCSamplingDiag(FILE *fout, FILE *find, TypeTree **tree, TypeFossilIntFeature *fi, double al, int burn, int gap, int iter, double prop, TypeModelParam *windSize, TypeModelParam *init, double probSpe, double probExt) {
	int *fol, *branch, f, i , s, n, *todo, nTodo = 0, nTree;
	double prob, **save, **saveTodo;
	TypeFossilFeature *fp;
	TypeModelParam param, *saveParam;
	
	param.birth = UNIF_RAND*init->birth;
	param.death = UNIF_RAND*init->death;
	param.fossil = UNIF_RAND*init->fossil;
	param.sampling = 1.;
	
	todo = (int*) malloc(sizeof(int)*(tree[0]->size/2+2));
	for(n=0; n<tree[0]->size; n++)
		if(tree[0]->node[n].child == NOSUCH && fi->status[n] == extinctNodeStatus)
			todo[nTodo++] = n;
	fp = sampleFossil(fi, tree[0]->size);
	for(nTree=0; tree[nTree] != NULL; nTree++) {
		for(n=0; n<tree[nTree]->size; n++) {
			if(tree[nTree]->node[n].child == NOSUCH) {
				switch(fi->status[n]) {
					case contempNodeStatus:
						tree[nTree]->time[n] = tree[nTree]->maxTime;
					break;
					case unknownNodeStatus:

					break;
					case extinctNodeStatus:
						tree[nTree]->time[n] = fp->fossilList[fp->fossil[n]].time;
					break;
					default:
	//						error("Node %d has no status\n", n);
	;
				}
			} else
				tree[nTree]->time[n] = NO_TIME;
		}
	}
	fol = (int*) malloc(fp->size*sizeof(int));
	fillFolNode(fp, tree[0]->size, fol);
	branch = (int*) malloc(fp->size*sizeof(int));
	fillBranchNode(fp, tree[0]->size, branch);
	saveParam = (TypeModelParam*) malloc(iter*sizeof(TypeModelParam));
	save = (double**) malloc(fp->size*sizeof(double*));
	for(f=0; f<fp->size; f++)
		save[f] = (double*) malloc(iter*sizeof(double));
	saveTodo = (double**) malloc(nTodo*sizeof(double*));
	for(n=0; n<nTodo; n++)
		saveTodo[n] = (double*) malloc(iter*sizeof(double));
	prob = getLogDensitySum(tree, nTree, fp, &param);
	for(i=0; i<burn; i++) {
		prob = update(tree, nTree, fp, fi, fol, branch, prob, al, prop, &param, windSize, probSpe, probExt);
if(i%10 == 0) {
	fprintf(stderr, "it %d\r", i); fflush(stderr);
}
	}
	for(s=0; s<iter; s++) {
		int j, f;
		saveParam[s] = param;
		for(f=0; f<fp->size; f++)
			save[f][s] = fp->fossilList[f].time;
		for(n=0; n<nTodo; n++)
			saveTodo[n][s] = tree[0]->time[todo[n]];
		for(j=0; j<gap; j++) {
			prob = update(tree, nTree, fp, fi, fol, branch, prob, al, prop, &param, windSize, probSpe, probExt);
		}
if(s%10 == 0) {
	fprintf(stderr, "it %d\r", s); fflush(stderr);
}
	}
	for(s=0; s<iter; s++)
		fprintf(fout, "%d\t%lf\n",  s+burn+1, saveParam[s].birth);
	for(s=0; s<iter; s++)
		fprintf(fout, "%d\t%lf\n",  s+burn+1, saveParam[s].death);
	for(s=0; s<iter; s++)
		fprintf(fout, "%d\t%lf\n",  s+burn+1, saveParam[s].fossil);
	free((void*)saveParam);
	for(f=0; f<fp->size; f++)
		for(s=0; s<iter; s++)
			fprintf(fout, "%d\t%lf\n",  s+burn+1, save[f][s]);
	for(f=0; f<fp->size; f++)
		free((void*)save[f]);
	free((void*)save);
	for(n=0; n<nTodo; n++)
		for(s=0; s<iter; s++)
			fprintf(fout, "%d\t%lf\n",  s+burn+1, saveTodo[n][s]);
	for(n=0; n<nTodo; n++)
		free((void*)saveTodo[n]);
	free((void*)saveTodo);
	int off = 0;
	fprintf(find, "birth\t%d\t%d\n", off*iter+1, (off+1)*iter); off++;
	fprintf(find, "death\t%d\t%d\n", off*iter+1, (off+1)*iter); off++;
	fprintf(find, "fossil\t%d\t%d\n", off*iter+1, (off+1)*iter); off++;
	for(f=0; f<fp->size; f++) {
		fprintf(find, "%d", f+1);
		if(branch[f] != NOSUCH && tree[0]->name[branch[f]] != NULL)
			fprintf(find, "_%s_%d_%d", tree[0]->name[branch[f]], (int) round(fi->fossilIntList[f].fossilInt.inf), (int) round(fi->fossilIntList[f].fossilInt.sup));
		fprintf(find, "\t%d\t%d\n", off*iter+1, (off+1)*iter); off++;
	}
	for(n=0; n<nTodo; n++) {
		if(tree[0]->name[todo[n]] != NULL)
			fprintf(find, "%s", tree[0]->name[todo[n]]);
		else
			fprintf(find, "node_%d", todo[n]);
		fprintf(find, "\t%d\t%d\n", off*iter+1, (off+1)*iter); off++;
	}
	free((void*)todo);
	free((void*)branch);
	free((void*)fol);
}

TypeDistribution *MCMCSamplingMultipleDistSingleTree(int *node, int sizeNode, TypeTree *tree, TypeFossilIntFeature *fi, double step, double al, int burn, int gap, int iter, double prop, TypeModelParam *windSize, TypeModelParam *init, double probSpe, double probExt) {
	int *fol, *branch, i, j, s, n;
	double prob, *logCond;
	TypeFossilFeature *fp;
	TypeDistribution *d, *logSTmp, *logDTmp;
	TypeModelParam param;
	
	param.birth = UNIF_RAND*init->birth;
	param.death = UNIF_RAND*init->death;
	param.fossil = UNIF_RAND*init->fossil;
	param.sampling = 1.;

	d = (TypeDistribution*) malloc(tree->size*sizeof(TypeDistribution));
	for(n=0; n<tree->size; n++) {
		d[n].size = 0;
		d[n].item = NULL;
	}
	fp = sampleFossil(fi, tree->size);
	for(n=0; n<tree->size; n++) {
		if(tree->node[n].child == NOSUCH) {
			switch(fi->status[n]) {
				case contempNodeStatus:
					tree->time[n] = tree->maxTime;
				break;
				case unknownNodeStatus:

				break;
				case extinctNodeStatus:
					tree->time[n] = fp->fossilList[fp->fossil[n]].time;
				break;
				default:
//						error("Node %d has no status\n", n);
;
			}
		} else
			tree->time[n] = NO_TIME;
	}
	fol = (int*) malloc(fp->size*sizeof(int));
	fillFolNode(fp, tree->size, fol);
	branch = (int*) malloc(fp->size*sizeof(int));
	fillBranchNode(fp, tree->size, branch);

	prob = getLogDensitySum(&tree, 1, fp, &param);
	for(i=0; i<burn; i++) {
		prob = update(&tree, 1, fp, fi, fol, branch, prob, al, prop, &param, windSize, probSpe, probExt);
if(i%10 == 0) {
	fprintf(stderr, "it %-7d\r", i); fflush(stderr);
}
	}
	logCond = (double*) malloc(sizeNode*sizeof(double));
	logSTmp = (TypeDistribution*) malloc(sizeNode*sizeof(TypeDistribution));
	logDTmp = (TypeDistribution*) malloc(sizeNode*sizeof(TypeDistribution));
	for(n=0; n<sizeNode; n++) {
		double min, max;
		min = getMaxFossilIntTimeToNode(node[n], tree, fi);
//		max = utils_MIN(getMinFossilIntTimeFromNode(node[n], tree, fi), getMinTimeFromNode(node[n], tree));
		max = getMinFossilIntTimeFromNode(node[n], tree, fi);
		min = floor(min/step)*step;
		max = (ceil(max/step))*step;
		logSTmp[n].size = (int) ceil((max-min)/step)+1;
		logSTmp[n].item = (TypeDistributionItem*) malloc(logSTmp[n].size*sizeof(TypeDistributionItem));
		logDTmp[n].size = logSTmp[n].size;
		logDTmp[n].item = (TypeDistributionItem*) malloc(logSTmp[n].size*sizeof(TypeDistributionItem));
		for(j=0; j<logSTmp[n].size; j++) {
			logDTmp[n].item[j].val = min + ((double)j)*step;
			logSTmp[n].item[j].val = min + ((double)j)*step;
			logSTmp[n].item[j].dens = NEG_INFTY;
		}
	}
	for(s=0; s<iter; s++) {
		int j;
		int cont=1, nT=0, i = 0;
		while(cont) {
			pthread_mutex_lock(&mutexN);
			while(i < sizeNode && nT < maxT) {
				pthread_t thread;	
				int ret = 0;
				TypeThreadParameterDist *tp;
				tp = (TypeThreadParameterDist*) malloc(sizeof(TypeThreadParameterDist));
				tp->tree = tree;
				tp->ffe = fp;
				tp->param = &param;
				tp->node = node[i];
				tp->logD = &(logDTmp[i]);
				tp->logCond = &(logCond[i]);
				tp->mutex_number = &mutexN;
				tp->cond_number = &condN;
				tp->number = &nT;
				if((ret = pthread_create(&thread, NULL, threadComputeDistribution, (void*) tp)) == 0) {
					int err;
					if((err = pthread_detach(thread)) == 0) {
						nT++; i++;
					} else {
						warning("Error %d while detaching thread: %s\n", err, (char*) strerror(err));
					}
				} else
					warning("Error %d while creating thread: %s\n", ret, (char*) strerror(ret));
			}
			cont = (nT > 0);
			if(cont)
				pthread_cond_wait(&condN, &mutexN);
			pthread_mutex_unlock(&mutexN);
		}
		for(n=0; n<sizeNode; n++) {
			for(j=0; j<logSTmp[n].size; j++)
				logSTmp[n].item[j].dens = logSumLog(logSTmp[n].item[j].dens, logDTmp[n].item[j].dens-logCond[n]);
		}
		for(j=0; j<gap; j++) {
			prob = update(&tree, 1, fp, fi, fol, branch, prob, al, prop, &param, windSize, probSpe, probExt);
		}
if(s%10 == 0) {
	fprintf(stderr, "it %-7d\r", s); fflush(stderr);
}
	}
	for(n=0; n<sizeNode; n++) {
		d[node[n]].size = logSTmp[n].size;
		d[node[n]].item = (TypeDistributionItem*) malloc(logSTmp[n].size*sizeof(TypeDistributionItem));
		for(j=0; j<logSTmp[n].size; j++) {
			d[node[n]].item[j].val = logSTmp[n].item[j].val;
			d[node[n]].item[j].dens = exp(logSTmp[n].item[j].dens-(log((double)iter)));
		}
		free((void*) logSTmp[n].item);
		free((void*) logDTmp[n].item);
	}
	free((void*) logSTmp);
	free((void*) logDTmp);
	free((void*) logCond);
	free((void*)branch);
	free((void*)fol);
	return d;
}


void *threadComputeExtSingleTree(void *data) {
	int j;
	for(j=0; j<((TypeThreadParameterExtSingleTree*)data)->logD->size; j++) {
		if(((TypeThreadParameterExtSingleTree*)data)->logD->item[j].val>=((TypeThreadParameterExtSingleTree*)data)->tree->time[((TypeThreadParameterExtSingleTree*)data)->taxa])
			((TypeThreadParameterExtSingleTree*)data)->logD->item[j].dens = getLogProbExtinctCond(((TypeThreadParameterExtSingleTree*)data)->logD->item[j].val, ((TypeThreadParameterExtSingleTree*)data)->tree->time[((TypeThreadParameterExtSingleTree*)data)->taxa], ((TypeThreadParameterExtSingleTree*)data)->tree->maxTime, ((TypeThreadParameterExtSingleTree*)data)->pext);
		else
			((TypeThreadParameterExtSingleTree*)data)->logD->item[j].dens = NEG_INFTY;
	}
	pthread_mutex_lock(((TypeThreadParameterExtSingleTree*)data)->mutex_number);
		(*((TypeThreadParameterExtSingleTree*)data)->number)--;
		pthread_cond_signal(((TypeThreadParameterExtSingleTree*)data)->cond_number);
	pthread_mutex_unlock(((TypeThreadParameterExtSingleTree*)data)->mutex_number);
	free(data);
	return NULL;
}

TypeDistribution *MCMCSamplingExtinctionMultipleDistSingleTree(int *node, int sizeNode, double maxDisplayTime, TypeTree *tree, TypeFossilIntFeature *fi, double step, double al, int burn, int gap, int iter, double prop, TypeModelParam *windSize, TypeModelParam *init, double probSpe, double probExt) {
	int *fol, *branch, i, j, s, n;
	double prob;
	TypeFossilFeature *fp;
	TypeDistribution *d, *logSTmp, *logDTmp;
	TypeModelParam param;
	TypeExtendedModelParam pext;
	
	param.birth = UNIF_RAND*init->birth;
	param.death = UNIF_RAND*init->death;
	param.fossil = UNIF_RAND*init->fossil;
	param.sampling = 1.;

	d = (TypeDistribution*) malloc(tree->size*sizeof(TypeDistribution));
	for(n=0; n<tree->size; n++) {
		d[n].size = 0;
		d[n].item = NULL;
	}
	fp = sampleFossil(fi, tree->size);
	for(n=0; n<tree->size; n++) {
		if(tree->node[n].child == NOSUCH) {
			switch(fi->status[n]) {
				case contempNodeStatus:
					tree->time[n] = tree->maxTime;
				break;
				case unknownNodeStatus:

				break;
				case extinctNodeStatus:
					tree->time[n] = fp->fossilList[fp->fossil[n]].time;
				break;
				default:
//						error("Node %d has no status\n", n);
;
			}
		} else
			tree->time[n] = NO_TIME;
	}
	fol = (int*) malloc(fp->size*sizeof(int));
	fillFolNode(fp, tree->size, fol);
	branch = (int*) malloc(fp->size*sizeof(int));
	fillBranchNode(fp, tree->size, branch);

	prob = getLogDensitySum(&tree, 1, fp, &param);
	for(i=0; i<burn; i++) {
		prob = update(&tree, 1, fp, fi, fol, branch, prob, al, prop, &param, windSize, probSpe, probExt);
if(i%10 == 0) {
	fprintf(stderr, "it %-7d\r", i); fflush(stderr);
}
	}
	logSTmp = (TypeDistribution*) malloc(sizeNode*sizeof(TypeDistribution));
	logDTmp = (TypeDistribution*) malloc(sizeNode*sizeof(TypeDistribution));
	for(n=0; n<sizeNode; n++) {
		double min, max;
		min = getMaxFossilIntTimeToNode(node[n], tree, fi);
//		max = tree->maxTime;
		max = maxDisplayTime;
		min = floor(min/step)*step;
		max = ceil(max/step)*step;
		logSTmp[n].size = (int) ceil((max-min)/step)+1;
		logSTmp[n].item = (TypeDistributionItem*) malloc(logSTmp[n].size*sizeof(TypeDistributionItem));
		logDTmp[n].size = logSTmp[n].size;
		logDTmp[n].item = (TypeDistributionItem*) malloc(logDTmp[n].size*sizeof(TypeDistributionItem));
		for(j=0; j<logSTmp[n].size; j++) {
			logDTmp[n].item[j].val = min + ((double)j)*step;
			logSTmp[n].item[j].val = logDTmp[n].item[j].val;
			logSTmp[n].item[j].dens = NEG_INFTY;
		}
	}
	for(s=0; s<iter; s++) {
		int j;
		int cont=1, nT=0, i = 0;
		pext =  getExtendedModelParam(&param);
		while(cont) {
			pthread_mutex_lock(&mutexN);
			while(i < sizeNode && nT < maxT) {
				pthread_t thread;	
				int ret = 0;
				TypeThreadParameterExtSingleTree *tp;
				tp = (TypeThreadParameterExtSingleTree*) malloc(sizeof(TypeThreadParameterExtSingleTree));
				tp->tree = tree;
				tp->ffe = fp;
				tp->pext = &pext;
				tp->taxa = node[i];
				tp->logD = &(logDTmp[i]);
				tp->mutex_number = &mutexN;
				tp->cond_number = &condN;
				tp->number = &nT;
				if((ret = pthread_create(&thread, NULL, threadComputeExtSingleTree, (void*) tp)) == 0) {
					int err;
					if((err = pthread_detach(thread)) == 0) {
						nT++; i++;
					} else {
						warning("Error %d while detaching thread: %s\n", err, (char*) strerror(err));
					}
				} else
					warning("Error %d while creating thread: %s\n", ret, (char*) strerror(ret));
			}
			cont = (nT > 0);
			if(cont)
				pthread_cond_wait(&condN, &mutexN);
			pthread_mutex_unlock(&mutexN);
		}
		for(n=0; n<sizeNode; n++) {
			for(j=0; j<logSTmp[n].size; j++)
				if(logDTmp[n].item[j].dens > NEG_INFTY)
					logSTmp[n].item[j].dens = logSumLog(logSTmp[n].item[j].dens, logDTmp[n].item[j].dens);
		}
		for(j=0; j<gap; j++) {
			prob = update(&tree, 1, fp, fi, fol, branch, prob, al, prop, &param, windSize, probSpe, probExt);
		}
if(s%10 == 0) {
	fprintf(stderr, "it %-7d\r", s); fflush(stderr);
}
	}

		for(n=0; n<sizeNode; n++) {
			FILE *fo;
			char outputDistribution[100];
			if(tree->name[node[n]] != NULL)
				sprintf(outputDistribution, "DisTmp_%s.csv", tree->name[node[n]]);
			else
				sprintf(outputDistribution, "DisTmp_%d.csv", node[n]);
			for(j=0; j<logSTmp[n].size; j++)
				logDTmp[n].item[j].dens = exp(logDTmp[n].item[j].dens);
			if((fo = fopen(outputDistribution, "w"))) {
				fprintDistribution(fo, logDTmp[n]);
				fclose(fo);
			}
		}
		
	for(n=0; n<sizeNode; n++) {
		d[node[n]].size = logSTmp[n].size;
		d[node[n]].item = (TypeDistributionItem*) malloc(logSTmp[n].size*sizeof(TypeDistributionItem));
		for(j=0; j<logSTmp[n].size; j++) {
			d[node[n]].item[j].val = logSTmp[n].item[j].val;
			if(logSTmp[n].item[j].dens > NEG_INFTY)
				d[node[n]].item[j].dens = exp(logSTmp[n].item[j].dens-log((double)iter));
			else
				d[node[n]].item[j].dens = 0.;
		}
		free((void*) logSTmp[n].item);
		free((void*) logDTmp[n].item);
	}
	free((void*) logSTmp);
	free((void*) logDTmp);
	free((void*)branch);
	free((void*)fol);
	return d;
}

TypeDistribution *MCMCSamplingExtinctionDist(FILE *fout, FILE *find, int **clade, TypeTree **tree, TypeFossilIntFeature *fi, double maxDisplayTime, double step, double al, int burn, int gap, int iter, double prop, TypeModelParam *windSize, TypeModelParam *init, double probSpe, double probExt) {
	int *fol, *branch, i, j, f, s, n, nTree = 0, nClade = 0, nExt = 0, *ext;
	double prob, **save, **saveExt, min;
	TypeFossilFeature *fp;
	TypeDistribution *d, *logSTmp;
	TypeModelParam param, *saveParam;
	TypeExtendedModelParam pext;
	
	param.birth = UNIF_RAND*init->birth;
	param.death = UNIF_RAND*init->death;
	param.fossil = UNIF_RAND*init->fossil;
	param.sampling = 1.;
	fp = sampleFossil(fi, tree[0]->size);
	for(nTree=0; tree[nTree] != NULL; nTree++) {
		for(n=0; n<tree[nTree]->size; n++) {
			if(tree[nTree]->node[n].child == NOSUCH) {
				switch(fi->status[n]) {
					case contempNodeStatus:
						tree[nTree]->time[n] = tree[nTree]->maxTime;
					break;
					case unknownNodeStatus:

					break;
					case extinctNodeStatus:
						tree[nTree]->time[n] = fp->fossilList[fp->fossil[n]].time;
					break;
					default:
	//						error("Node %d has no status\n", n);
	;
				}
			} else
				tree[nTree]->time[n] = NO_TIME;
		}
	}
	ext = (int*) malloc(tree[0]->size*sizeof(int));
	for(n=0; n<tree[0]->size; n++)
		if(tree[0]->node[n].child == NOSUCH && fi->status[n] == extinctNodeStatus)
			ext[nExt++] = n;
	fol = (int*) malloc(fp->size*sizeof(int));
	fillFolNode(fp, tree[0]->size, fol);
	branch = (int*) malloc(fp->size*sizeof(int));
	fillBranchNode(fp, tree[0]->size, branch);
	saveParam = (TypeModelParam*) malloc(iter*sizeof(TypeModelParam));
	save = (double**) malloc(fp->size*sizeof(double*));
	for(f=0; f<fp->size; f++)
		save[f] = (double*) malloc(iter*sizeof(double));
	saveExt = (double**) malloc(nExt*sizeof(double*));
	for(n=0; n<nExt; n++)
		saveExt[n] = (double*) malloc(iter*sizeof(double));
	prob = getLogDensitySum(tree, nTree, fp, &param);
	for(i=0; i<burn; i++) {
		prob = update(tree, nTree, fp, fi, fol, branch, prob, al, prop, &param, windSize, probSpe, probExt);
if(i%10 == 0) {
	fprintf(stderr, "it %-7d\r", i); fflush(stderr);
}
	}
	for(nClade=0; clade[nClade]!=NULL; nClade++)
		;
	min = getMaxFossilIntTimeToNode(ext[0], tree[0], fi);
	for(n=1; n<nExt; n++) {
		if(getMaxFossilIntTimeToNode(ext[n], tree[0], fi)<min)
			min = getMaxFossilIntTimeToNode(ext[n], tree[0], fi);
	}
	min = floor(min/step)*step;
	logSTmp = (TypeDistribution*) malloc(nClade*sizeof(TypeDistribution));
	logSTmp[0].size = (int) ceil((ceil(maxDisplayTime/step)*step-min)/step)+1;
	for(n=0; n<nClade; n++) {
		logSTmp[n].size = logSTmp[0].size;
		logSTmp[n].item = (TypeDistributionItem*) malloc(logSTmp[n].size*sizeof(TypeDistributionItem));
		for(j=0; j<logSTmp[n].size; j++) {
			logSTmp[n].item[j].val = min + ((double)j)*step;
			logSTmp[n].item[j].dens = NEG_INFTY;
		}
	}
	for(s=0; s<iter; s++) {
		int j;
		pext =  getExtendedModelParam(&param);
		int cont=1, nT=0, i = 0;
		while(cont) {
			pthread_mutex_lock(&mutexN);
			while(i < nClade && nT < maxT) {
				pthread_t thread;	
				int ret = 0;
				TypeThreadParameterExtMultipleTree *tp;
				tp = (TypeThreadParameterExtMultipleTree*) malloc(sizeof(TypeThreadParameterExtMultipleTree));
				tp->tree = tree[0];
				tp->ffe = fp;
				tp->pext = &pext;
				tp->param = &param;
				tp->clade = clade[i];
				tp->logD = &(logSTmp[i]);
				tp->mutex_number = &mutexN;
				tp->cond_number = &condN;
				tp->number = &nT;
				if((ret = pthread_create(&thread, NULL, threadComputeExtMultipleTree, (void*) tp)) == 0) {
					int err;
					if((err = pthread_detach(thread)) == 0) {
						nT++; i++;
					} else {
						warning("Error %d while detaching thread: %s\n", err, (char*) strerror(err));
					}
				} else
					warning("Error %d while creating thread: %s\n", ret, (char*) strerror(ret));
			}
			cont = (nT > 0);
			if(cont)
				pthread_cond_wait(&condN, &mutexN);
			pthread_mutex_unlock(&mutexN);
		}
		saveParam[s] = param;
		for(f=0; f<fp->size; f++)
			save[f][s] = fp->fossilList[f].time;
		for(n=0; n<nExt; n++)
			saveExt[n][s] = tree[0]->time[ext[n]];
		for(j=0; j<gap; j++) {
			prob = update(tree, nTree, fp, fi, fol, branch, prob, al, prop, &param, windSize, probSpe, probExt);
		}
if(s%10 == 0) {
	fprintf(stderr, "it %-7d\r", s); fflush(stderr);
}
	}
	d = (TypeDistribution*) malloc(nExt*sizeof(TypeDistribution));
	for(n=0; n<nClade; n++) {
		d[n].size = logSTmp[n].size;
		d[n].item = (TypeDistributionItem*) malloc(logSTmp[n].size*sizeof(TypeDistributionItem));
		for(j=0; j<logSTmp[n].size; j++) {
			d[n].item[j].val = logSTmp[n].item[j].val;
			if(logSTmp[n].item[j].dens > NEG_INFTY)
				d[n].item[j].dens = exp(logSTmp[n].item[j].dens-log((double)iter));
			else
				d[n].item[j].dens = 0.;
		}
		free((void*) logSTmp[n].item);
	}
	free((void*) logSTmp);
	for(s=0; s<iter; s++)
		fprintf(fout, "%d\t%lf\n",  s+burn+1, saveParam[s].birth);
	for(s=0; s<iter; s++)
		fprintf(fout, "%d\t%lf\n",  s+burn+1, saveParam[s].death);
	for(s=0; s<iter; s++)
		fprintf(fout, "%d\t%lf\n",  s+burn+1, saveParam[s].fossil);
	free((void*)saveParam);
	for(f=0; f<fp->size; f++)
		for(s=0; s<iter; s++)
			fprintf(fout, "%d\t%lf\n",  s+burn+1, save[f][s]);
	for(f=0; f<fp->size; f++)
		free((void*)save[f]);
	free((void*)save);
	for(n=0; n<nExt; n++)
		for(s=0; s<iter; s++)
			fprintf(fout, "%d\t%lf\n",  s+burn+1, saveExt[n][s]);
	for(n=0; n<nExt; n++)
		free((void*)saveExt[n]);
	free((void*)saveExt);
	int off = 0;
	fprintf(find, "birth\t%d\t%d\n", off*iter+1, (off+1)*iter); off++;
	fprintf(find, "death\t%d\t%d\n", off*iter+1, (off+1)*iter); off++;
	fprintf(find, "fossil\t%d\t%d\n", off*iter+1, (off+1)*iter); off++;
	for(f=0; f<fp->size; f++) {
		fprintf(find, "%d", f+1);
		if(branch[f] != NOSUCH && tree[0]->name[branch[f]] != NULL)
			fprintf(find, "_%s_%d_%d", tree[0]->name[branch[f]], (int) round(fi->fossilIntList[f].fossilInt.inf), (int) round(fi->fossilIntList[f].fossilInt.sup));
		fprintf(find, "\t%d\t%d\n", off*iter+1, (off+1)*iter); off++;
	}
	for(n=0; n<nExt; n++) {
		if(tree[0]->name[ext[n]] != NULL)
			fprintf(find, "%s", tree[0]->name[ext[n]]);
		else
			fprintf(find, "node_%d", ext[n]);
		fprintf(find, "\t%d\t%d\n", off*iter+1, (off+1)*iter); off++;
	}
	free((void*)branch);
	free((void*)fol);
	return d;
}


void *threadComputeExtMultipleTree(void *data) {
	int j;
	for(j=0; j<((TypeThreadParameterExtMultipleTree*)data)->logD->size; j++) {
		int i;
		double tmp = 0.;
		for(i=0; ((TypeThreadParameterExtMultipleTree*)data)->clade[i]!=END_LIST_INT && ((TypeThreadParameterExtMultipleTree*)data)->logD->item[j].val>=((TypeThreadParameterExtMultipleTree*)data)->tree->time[((TypeThreadParameterExtMultipleTree*)data)->clade[i]]; i++)
			tmp += getLogProbExtinctCond(((TypeThreadParameterExtMultipleTree*)data)->logD->item[j].val, ((TypeThreadParameterExtMultipleTree*)data)->tree->time[((TypeThreadParameterExtMultipleTree*)data)->clade[i]], ((TypeThreadParameterExtMultipleTree*)data)->tree->maxTime, ((TypeThreadParameterExtMultipleTree*)data)->pext);
		if(((TypeThreadParameterExtMultipleTree*)data)->clade[i] != END_LIST_INT)
			tmp = NEG_INFTY;
		if(tmp != NEG_INFTY)
			((TypeThreadParameterExtMultipleTree*)data)->logD->item[j].dens = logSumLog(((TypeThreadParameterExtMultipleTree*)data)->logD->item[j].dens, tmp);

	}
	pthread_mutex_lock(((TypeThreadParameterExtMultipleTree*)data)->mutex_number);
		(*((TypeThreadParameterExtMultipleTree*)data)->number)--;
		pthread_cond_signal(((TypeThreadParameterExtMultipleTree*)data)->cond_number);
	pthread_mutex_unlock(((TypeThreadParameterExtMultipleTree*)data)->mutex_number);
	free(data);
	return NULL;
}


TypeDistribution *MCMCSamplingExtinctionMultipleDistSingleTreeFossilFixed(int *node, int sizeNode, double maxDisplayTime, TypeTree *tree, TypeFossilFeature *fp, double step, int burn, int gap, int iter, TypeModelParam *windSize, TypeModelParam *init, double probSpe, double probExt) {
	int i, j, s, n;
	double prob;
	TypeDistribution *d, *logSTmp, *logDTmp;
	TypeModelParam param;
	TypeExtendedModelParam pext;
	
	param.birth = UNIF_RAND*init->birth;
	param.death = UNIF_RAND*init->death;
	param.fossil = UNIF_RAND*init->fossil;
	param.sampling = 1.;

	d = (TypeDistribution*) malloc(sizeNode*sizeof(TypeDistribution));
	for(n=0; n<sizeNode; n++) {
		d[n].size = 0;
		d[n].item = NULL;
	}
	for(n=0; n<tree->size; n++)
		if(tree->node[n].child != NOSUCH)
			tree->time[n] = NO_TIME;
	prob = getLogDensitySum(&tree, 1, fp, &param);
	for(i=0; i<burn; i++) {
		prob = updateFossilFixed(&tree, 1, fp, prob, &param, windSize, probSpe, probExt);
if(i%10 == 0) {
	fprintf(stderr, "it %-7d\r", i); fflush(stderr);
}
	}
	logSTmp = (TypeDistribution*) malloc(sizeNode*sizeof(TypeDistribution));
	logDTmp = (TypeDistribution*) malloc(sizeNode*sizeof(TypeDistribution));
	for(n=0; n<sizeNode; n++) {
		double min, max;
		min = getMaxFossilTimeToNode(node[n], tree, fp);
//		max = tree->maxTime;
		max = maxDisplayTime;
		min = floor(min/step)*step;
		max = ceil(max/step)*step;
		logSTmp[n].size = (int) ceil((max-min)/step)+1;
		logSTmp[n].item = (TypeDistributionItem*) malloc(logSTmp[n].size*sizeof(TypeDistributionItem));
		logDTmp[n].size = logSTmp[n].size;
		logDTmp[n].item = (TypeDistributionItem*) malloc(logDTmp[n].size*sizeof(TypeDistributionItem));
		for(j=0; j<logSTmp[n].size; j++) {
			logDTmp[n].item[j].val = min + ((double)j)*step;
			logSTmp[n].item[j].val = logDTmp[n].item[j].val;
			logSTmp[n].item[j].dens = NEG_INFTY;
		}
	}
	for(s=0; s<iter; s++) {
		int j;
		int cont=1, nT=0, i = 0;
		pext =  getExtendedModelParam(&param);
		while(cont) {
			pthread_mutex_lock(&mutexN);
			while(i < sizeNode && nT < maxT) {
				pthread_t thread;	
				int ret = 0;
				TypeThreadParameterExtSingleTree *tp;
				tp = (TypeThreadParameterExtSingleTree*) malloc(sizeof(TypeThreadParameterExtSingleTree));
				tp->tree = tree;
				tp->ffe = fp;
				tp->pext = &pext;
				tp->taxa = node[i];
				tp->logD = &(logDTmp[i]);
				tp->mutex_number = &mutexN;
				tp->cond_number = &condN;
				tp->number = &nT;
				if((ret = pthread_create(&thread, NULL, threadComputeExtSingleTree, (void*) tp)) == 0) {
					int err;
					if((err = pthread_detach(thread)) == 0) {
						nT++; i++;
					} else {
						warning("Error %d while detaching thread: %s\n", err, (char*) strerror(err));
					}
				} else
					warning("Error %d while creating thread: %s\n", ret, (char*) strerror(ret));
			}
			cont = (nT > 0);
			if(cont)
				pthread_cond_wait(&condN, &mutexN);
			pthread_mutex_unlock(&mutexN);
		}
		for(n=0; n<sizeNode; n++) {
			for(j=0; j<logSTmp[n].size; j++)
				if(logDTmp[n].item[j].dens > NEG_INFTY)
					logSTmp[n].item[j].dens = logSumLog(logSTmp[n].item[j].dens, logDTmp[n].item[j].dens);
		}
		for(j=0; j<gap; j++) {
			prob = updateFossilFixed(&tree, 1, fp, prob, &param, windSize, probSpe, probExt);
		}
if(s%10 == 0) {
	fprintf(stderr, "it %-7d\r", s); fflush(stderr);
}
	}
	for(n=0; n<sizeNode; n++) {
		d[n].size = logSTmp[n].size;
		d[n].item = (TypeDistributionItem*) malloc(logSTmp[n].size*sizeof(TypeDistributionItem));
		for(j=0; j<logSTmp[n].size; j++) {
			d[n].item[j].val = logSTmp[n].item[j].val;
			if(logSTmp[n].item[j].dens > NEG_INFTY)
				d[n].item[j].dens = exp(logSTmp[n].item[j].dens-log((double)iter));
			else
				d[n].item[j].dens = 0.;
		}
		free((void*) logSTmp[n].item);
		free((void*) logDTmp[n].item);
	}
	free((void*) logSTmp);
	free((void*) logDTmp);
	return d;
}

TypeDistribution *MCMCSamplingExtinctionMultipleDistSingleTreeFossilFixedX(FILE *fout, FILE *find, int *node, int sizeNode, double maxDisplayTime, TypeTree *tree, TypeFossilFeature *fp, double step, int burn, int gap, int iter, TypeModelParam *windSize, TypeModelParam *init, double probSpe, double probExt) {
	int i, j, s, n;
	double prob;
	TypeDistribution *d, *logSTmp, *logDTmp;
	TypeModelParam param, *saveParam;
	TypeExtendedModelParam pext;
	
	param.birth = UNIF_RAND*init->birth;
	param.death = UNIF_RAND*init->death;
	param.fossil = UNIF_RAND*init->fossil;
	param.sampling = 1.;

	saveParam = (TypeModelParam*) malloc(iter*sizeof(TypeModelParam));
	d = (TypeDistribution*) malloc(sizeNode*sizeof(TypeDistribution));
	for(n=0; n<sizeNode; n++) {
		d[n].size = 0;
		d[n].item = NULL;
	}
	for(n=0; n<tree->size; n++)
		if(tree->node[n].child != NOSUCH)
			tree->time[n] = NO_TIME;
	prob = getLogDensitySum(&tree, 1, fp, &param);
	for(i=0; i<burn; i++) {
		prob = updateFossilFixed(&tree, 1, fp, prob, &param, windSize, probSpe, probExt);
if(i%10 == 0) {
	fprintf(stderr, "it %-7d\r", i); fflush(stderr);
}
	}
	logSTmp = (TypeDistribution*) malloc(sizeNode*sizeof(TypeDistribution));
	logDTmp = (TypeDistribution*) malloc(sizeNode*sizeof(TypeDistribution));
	for(n=0; n<sizeNode; n++) {
		double min, max;
		min = getMaxFossilTimeToNode(node[n], tree, fp);
//		max = tree->maxTime;
		max = maxDisplayTime;
		min = floor(min/step)*step;
		max = ceil(max/step)*step;
		logSTmp[n].size = (int) ceil((max-min)/step)+1;
		logSTmp[n].item = (TypeDistributionItem*) malloc(logSTmp[n].size*sizeof(TypeDistributionItem));
		logDTmp[n].size = logSTmp[n].size;
		logDTmp[n].item = (TypeDistributionItem*) malloc(logDTmp[n].size*sizeof(TypeDistributionItem));
		for(j=0; j<logSTmp[n].size; j++) {
			logDTmp[n].item[j].val = min + ((double)j)*step;
			logSTmp[n].item[j].val = logDTmp[n].item[j].val;
			logSTmp[n].item[j].dens = NEG_INFTY;
		}
	}
	for(s=0; s<iter; s++) {
		int j;
		int cont=1, nT=0, i = 0;
		pext =  getExtendedModelParam(&param);
		while(cont) {
			pthread_mutex_lock(&mutexN);
			while(i < sizeNode && nT < maxT) {
				pthread_t thread;	
				int ret = 0;
				TypeThreadParameterExtSingleTree *tp;
				tp = (TypeThreadParameterExtSingleTree*) malloc(sizeof(TypeThreadParameterExtSingleTree));
				tp->tree = tree;
				tp->ffe = fp;
				tp->pext = &pext;
				tp->taxa = node[i];
				tp->logD = &(logDTmp[i]);
				tp->mutex_number = &mutexN;
				tp->cond_number = &condN;
				tp->number = &nT;
				if((ret = pthread_create(&thread, NULL, threadComputeExtSingleTree, (void*) tp)) == 0) {
					int err;
					if((err = pthread_detach(thread)) == 0) {
						nT++; i++;
					} else {
						warning("Error %d while detaching thread: %s\n", err, (char*) strerror(err));
					}
				} else
					warning("Error %d while creating thread: %s\n", ret, (char*) strerror(ret));
			}
			cont = (nT > 0);
			if(cont)
				pthread_cond_wait(&condN, &mutexN);
			pthread_mutex_unlock(&mutexN);
		}
		for(n=0; n<sizeNode; n++) {
			for(j=0; j<logSTmp[n].size; j++)
				if(logDTmp[n].item[j].dens > NEG_INFTY)
					logSTmp[n].item[j].dens = logSumLog(logSTmp[n].item[j].dens, logDTmp[n].item[j].dens);
		}
		saveParam[s] = param;
		for(j=0; j<gap; j++) {
			prob = updateFossilFixed(&tree, 1, fp, prob, &param, windSize, probSpe, probExt);
		}
if(s%10 == 0) {
	fprintf(stderr, "it %-7d\r", s); fflush(stderr);
}
	}
	for(s=0; s<iter; s++)
		fprintf(fout, "%d\t%lf\n",  s+burn+1, saveParam[s].birth);
	for(s=0; s<iter; s++)
		fprintf(fout, "%d\t%lf\n",  s+burn+1, saveParam[s].death);
	for(s=0; s<iter; s++)
		fprintf(fout, "%d\t%lf\n",  s+burn+1, saveParam[s].fossil);
	free((void*)saveParam);
	int off = 0;
	fprintf(find, "birth\t%d\t%d\n", off*iter+1, (off+1)*iter); off++;
	fprintf(find, "death\t%d\t%d\n", off*iter+1, (off+1)*iter); off++;
	fprintf(find, "fossil\t%d\t%d\n", off*iter+1, (off+1)*iter); off++;
	for(n=0; n<sizeNode; n++) {
		d[n].size = logSTmp[n].size;
		d[n].item = (TypeDistributionItem*) malloc(logSTmp[n].size*sizeof(TypeDistributionItem));
		for(j=0; j<logSTmp[n].size; j++) {
			d[n].item[j].val = logSTmp[n].item[j].val;
			if(logSTmp[n].item[j].dens > NEG_INFTY)
				d[n].item[j].dens = exp(logSTmp[n].item[j].dens-log((double)iter));
			else
				d[n].item[j].dens = 0.;
		}
		free((void*) logSTmp[n].item);
		free((void*) logDTmp[n].item);
	}
	free((void*) logSTmp);
	free((void*) logDTmp);
	return d;
}

void *threadComputeExtQuantSingleTree(void *data) {
	*(((TypeThreadParameterExtQuantSingleTree*)data)->quant) = getQuantile(((TypeThreadParameterExtQuantSingleTree*)data)->order, 0.01, ((TypeThreadParameterExtQuantSingleTree*)data)->tree->time[((TypeThreadParameterExtQuantSingleTree*)data)->taxa], ((TypeThreadParameterExtQuantSingleTree*)data)->tree->maxTime, ((TypeThreadParameterExtQuantSingleTree*)data)->pext);
	pthread_mutex_lock(((TypeThreadParameterExtQuantSingleTree*)data)->mutex_number);
		(*((TypeThreadParameterExtQuantSingleTree*)data)->number)--;
		pthread_cond_signal(((TypeThreadParameterExtQuantSingleTree*)data)->cond_number);
	pthread_mutex_unlock(((TypeThreadParameterExtQuantSingleTree*)data)->mutex_number);
	free(data);
	return NULL;
}

double *MCMCSamplingExtinctionMultipleDistSingleTreeFossilFixedQuantileMean(FILE *fout, FILE *find, double order, int *node, int sizeNode, double maxDisplayTime, TypeTree *tree, TypeFossilFeature *fp, double step, int burn, int gap, int iter, TypeModelParam *windSize, TypeModelParam *init, double probSpe, double probExt) {
	int i, s, n;
	double prob, *quant, **quantTab;
	TypeModelParam param, *saveParam;
	TypeExtendedModelParam pext;
	
	param.birth = UNIF_RAND*init->birth;
	param.death = UNIF_RAND*init->death;
	param.fossil = UNIF_RAND*init->fossil;
	param.sampling = 1.;

	saveParam = (TypeModelParam*) malloc(iter*sizeof(TypeModelParam));
	quantTab = (double**) malloc(iter*sizeof(double*));
	for(i=0; i<iter; i++)
		quantTab[i] = (double*) malloc(sizeNode*sizeof(double));
	quant = (double*) malloc(sizeNode*sizeof(double));
	for(n=0; n<sizeNode; n++)
		quant[n] = 0.;
	for(n=0; n<tree->size; n++)
		if(tree->node[n].child != NOSUCH)
			tree->time[n] = NO_TIME;
	prob = getLogDensitySum(&tree, 1, fp, &param);
	for(i=0; i<burn; i++) {
		prob = updateFossilFixed(&tree, 1, fp, prob, &param, windSize, probSpe, probExt);
if(i%10 == 0) {
	fprintf(stderr, "it %-7d\r", i); fflush(stderr);
}
	}
	for(s=0; s<iter; s++) {
		int j;
		int cont=1, nT=0, i = 0;
		pext =  getExtendedModelParam(&param);
		while(cont) {
			pthread_mutex_lock(&mutexN);
			while(i < sizeNode && nT < maxT) {
				pthread_t thread;	
				int ret = 0;
				TypeThreadParameterExtQuantSingleTree *tp;
				tp = (TypeThreadParameterExtQuantSingleTree*) malloc(sizeof(TypeThreadParameterExtQuantSingleTree));
				tp->tree = tree;
				tp->ffe = fp;
				tp->pext = &pext;
				tp->taxa = node[i];
				tp->order = order;
				tp->quant = &(quantTab[s][i]);
				tp->mutex_number = &mutexN;
				tp->cond_number = &condN;
				tp->number = &nT;
				if((ret = pthread_create(&thread, NULL, threadComputeExtQuantSingleTree, (void*) tp)) == 0) {
					int err;
					if((err = pthread_detach(thread)) == 0) {
						nT++; i++;
					} else {
						warning("Error %d while detaching thread: %s\n", err, (char*) strerror(err));
					}
				} else
					warning("Error %d while creating thread: %s\n", ret, (char*) strerror(ret));
			}
			cont = (nT > 0);
			if(cont)
				pthread_cond_wait(&condN, &mutexN);
			pthread_mutex_unlock(&mutexN);
		}
		for(n=0; n<sizeNode; n++) {
			if(isnan(quant[n]+quantTab[s][n])) {
				printf("NAN %lf %lf %lf\n", quant[n], quantTab[s][n], tree->time[node[n]]);
			}
			quant[n] += quantTab[s][n];
			if(isnan(quantTab[s][n])) {
				printf("NAN %lf\n", tree->time[node[n]]);
			}
			if(isnan(quant[n])) {
				printf("NAN %lf\n", tree->time[node[n]]);
				exit(0);
			}
		}
		saveParam[s] = param;
		for(j=0; j<gap; j++) {
			prob = updateFossilFixed(&tree, 1, fp, prob, &param, windSize, probSpe, probExt);
		}
if(s%10 == 0) {
	fprintf(stderr, "it %-7d\r", s); fflush(stderr);
}
	}
	for(n=0; n<sizeNode; n++)
		quant[n] /= (double) iter;
	for(s=0; s<iter; s++)
		fprintf(fout, "%d\t%lf\n",  s+burn+1, saveParam[s].birth);
	for(s=0; s<iter; s++)
		fprintf(fout, "%d\t%lf\n",  s+burn+1, saveParam[s].death);
	for(s=0; s<iter; s++)
		fprintf(fout, "%d\t%lf\n",  s+burn+1, saveParam[s].fossil);
	free((void*)saveParam);
	for(n=0; n<sizeNode; n++)
		for(s=0; s<iter; s++)
			fprintf(fout, "%d\t%lf\n",  s+burn+1, quantTab[s][n]);
	for(s=0; s<iter; s++)
		free((void*)quantTab[s]);
	free((void*)quantTab);
	int off = 0;
	fprintf(find, "birth\t%d\t%d\n", off*iter+1, (off+1)*iter); off++;
	fprintf(find, "death\t%d\t%d\n", off*iter+1, (off+1)*iter); off++;
	fprintf(find, "fossil\t%d\t%d\n", off*iter+1, (off+1)*iter); off++;
	for(n=0; n<sizeNode; n++) {
		if(tree->name != NULL && tree->name[node[n]] != NULL)
			fprintf(find, "Quantile_%s\t%d\t%d\n", tree->name[node[n]], off*iter+1, (off+1)*iter);
		else
			fprintf(find, "Quantile_%d\t%d\t%d\n", node[n], off*iter+1, (off+1)*iter);
		off++;
	}
	return quant;
}

double *MCMCSamplingExtinctionMultipleDistSingleTreeFossilFixedQuantileQuantile(FILE *fout, FILE *find, double order, int *node, int sizeNode, double maxDisplayTime, TypeTree *tree, TypeFossilFeature *fp, double step, int burn, int gap, int iter, TypeModelParam *windSize, TypeModelParam *init, double probSpe, double probExt) {
	int i, s, n;
	double prob, *quant, **quantTab, *logLike, sumLog = NEG_INFTY;
	TypeModelParam param, *saveParam;
	TypeExtendedModelParam pext;
	
	param.birth = UNIF_RAND*init->birth;
	param.death = UNIF_RAND*init->death;
	param.fossil = UNIF_RAND*init->fossil;
	param.sampling = 1.;

	saveParam = (TypeModelParam*) malloc(iter*sizeof(TypeModelParam));
	logLike = (double*) malloc(iter*sizeof(double));
	quantTab = (double**) malloc(sizeNode*sizeof(double*));
	for(n=0; n<sizeNode; n++)
		quantTab[n] = (double*) malloc(iter*sizeof(double));
	quant = (double*) malloc(sizeNode*sizeof(double));
	for(n=0; n<tree->size; n++)
		if(tree->node[n].child != NOSUCH)
			tree->time[n] = NO_TIME;
	prob = getLogDensitySum(&tree, 1, fp, &param);
	for(i=0; i<burn; i++) {
		prob = updateFossilFixed(&tree, 1, fp, prob, &param, windSize, probSpe, probExt);
if(i%10 == 0) {
	fprintf(stderr, "it %-7d\r", i); fflush(stderr);
}
	}
	for(s=0; s<iter; s++) {
		int j;
		int cont=1, nT=0, i = 0;
		pext =  getExtendedModelParam(&param);
		logLike[s] = prob;
		sumLog = logSumLog(sumLog, logLike[s]);
		while(cont) {
			pthread_mutex_lock(&mutexN);
			while(i < sizeNode && nT < maxT) {
				pthread_t thread;	
				int ret = 0;
				TypeThreadParameterExtQuantSingleTree *tp;
				tp = (TypeThreadParameterExtQuantSingleTree*) malloc(sizeof(TypeThreadParameterExtQuantSingleTree));
				tp->tree = tree;
				tp->ffe = fp;
				tp->pext = &pext;
				tp->taxa = node[i];
				tp->order = order;
				tp->quant = &(quantTab[i][s]);
				tp->mutex_number = &mutexN;
				tp->cond_number = &condN;
				tp->number = &nT;
				if((ret = pthread_create(&thread, NULL, threadComputeExtQuantSingleTree, (void*) tp)) == 0) {
					int err;
					if((err = pthread_detach(thread)) == 0) {
						nT++; i++;
					} else {
						warning("Error %d while detaching thread: %s\n", err, (char*) strerror(err));
					}
				} else
					warning("Error %d while creating thread: %s\n", ret, (char*) strerror(ret));
			}
			cont = (nT > 0);
			if(cont)
				pthread_cond_wait(&condN, &mutexN);
			pthread_mutex_unlock(&mutexN);
		}
		saveParam[s] = param;
		for(j=0; j<gap; j++) {
			prob = updateFossilFixed(&tree, 1, fp, prob, &param, windSize, probSpe, probExt);
		}
if(s%10 == 0) {
	fprintf(stderr, "it %-7d\r", s); fflush(stderr);
}
	}
	for(s=0; s<iter; s++)
			logLike[s] -= sumLog;
	for(n=0; n<sizeNode; n++) {
		size_t *index = qsortindex(quantTab[n], iter, sizeof(double), compareDouble);
		double logCumul = 0.;
		for(s=iter-1; exp(logCumul) < 1.-order && s >= 0; s--)
			logSumLog(logCumul, logLike[index[s]]);
		quant[n] = quantTab[n][s];
		free((void*)index);
	}
	for(s=0; s<iter; s++)
		fprintf(fout, "%d\t%lf\n",  s+burn+1, saveParam[s].birth);
	for(s=0; s<iter; s++)
		fprintf(fout, "%d\t%lf\n",  s+burn+1, saveParam[s].death);
	for(s=0; s<iter; s++)
		fprintf(fout, "%d\t%lf\n",  s+burn+1, saveParam[s].fossil);
	free((void*)saveParam);
	for(n=0; n<sizeNode; n++) {
		for(s=0; s<iter; s++)
			fprintf(fout, "%d\t%lf\n",  s+burn+1, quantTab[n][s]);
		free((void*)quantTab[n]);
	}
	free((void*)quantTab);
	int off = 0;
	fprintf(find, "birth\t%d\t%d\n", off*iter+1, (off+1)*iter); off++;
	fprintf(find, "death\t%d\t%d\n", off*iter+1, (off+1)*iter); off++;
	fprintf(find, "fossil\t%d\t%d\n", off*iter+1, (off+1)*iter); off++;
	for(n=0; n<sizeNode; n++) {
		if(tree->name != NULL && tree->name[node[n]] != NULL)
			fprintf(find, "Quantile_%s\t%d\t%d\n", tree->name[node[n]], off*iter+1, (off+1)*iter);
		else
			fprintf(find, "Quantile_%d\t%d\t%d\n", node[n], off*iter+1, (off+1)*iter);
		off++;
	}
	return quant;
}

/*quant first sizeNode entries are "means", last sizeNode entries are "quantiles"*/
double *MCMCSamplingExtinctionMultipleDistSingleTreeFossilFixedQuantileQuantileMean(int tmp, FILE *fout, FILE *find, double order, int *node, int sizeNode, TypeTree *tree, TypeFossilFeature *fp, int burn, int gap, int iter, TypeModelParam *windSize, TypeModelParam *init, double probSpe, double probExt) {
	int i, s, n;
	double prob, *quant, **quantTab, *logLike, sumLog = NEG_INFTY;
	TypeModelParam param, *saveParam;
	TypeExtendedModelParam pext, *savePext;
	
	param.birth = UNIF_RAND*init->birth;
	param.death = UNIF_RAND*init->death;
	param.fossil = UNIF_RAND*init->fossil;
	param.sampling = 1.;
	if(fout != NULL && find != NULL)
		saveParam = (TypeModelParam*) malloc(iter*sizeof(TypeModelParam));
	savePext = (TypeExtendedModelParam*) malloc(iter*sizeof(TypeExtendedModelParam));
	logLike = (double*) malloc(iter*sizeof(double));
	quantTab = (double**) malloc(sizeNode*sizeof(double*));
	for(n=0; n<sizeNode; n++)
		quantTab[n] = (double*) malloc(iter*sizeof(double));
	quant = (double*) malloc(3*sizeNode*sizeof(double));
	for(n=0; n<sizeNode; n++)
		quant[n] = 0.;
	for(n=0; n<tree->size; n++)
		if(tree->node[n].child != NOSUCH)
			tree->time[n] = NO_TIME;
	prob = getLogDensitySumNoThread(&tree, 1, fp, &param);
	for(i=0; i<burn; i++) {
		prob = updateFossilFixed(&tree, 1, fp, prob, &param, windSize, probSpe, probExt);
//if(i%10 == 0) {
	//fprintf(stderr, "it %-7d\r", i); fflush(stderr);
//}
	}
	for(s=0; s<iter; s++) {
		int j;
		pext =  getExtendedModelParam(&param);
		logLike[s] = prob;
		sumLog = logSumLog(sumLog, logLike[s]);
		for(n=0; n<sizeNode; n++) {
			quantTab[n][s] = getQuantile(order, 0.01, tree->time[node[n]], tree->maxTime, &pext);
			quant[n] += quantTab[n][s];
		}
		if(fout != NULL && find != NULL)
			saveParam[s] = param;
		savePext[s] = pext;
		for(j=0; j<gap; j++) {
			prob = updateFossilFixed(&tree, 1, fp, prob, &param, windSize, probSpe, probExt);
		}
//if(s%10 == 0) {
	//fprintf(stderr, "it %-7d\r", s); fflush(stderr);
//}
	}
	for(n=0; n<sizeNode; n++)
		quant[n] /= (double) iter;
	for(s=0; s<iter; s++)
			logLike[s] -= sumLog;
	for(n=0; n<sizeNode; n++) {
		size_t *index = qsortindex(quantTab[n], iter, sizeof(double), compareDouble);
		double logCumul = 0.;
		for(s=iter-1; exp(logCumul) < 1.-order && s >= 0; s--)
			logSumLog(logCumul, logLike[index[s]]);
		quant[n+sizeNode] = quantTab[n][s];
		free((void*)index);
	}
	for(n=0; n<sizeNode; n++)
		quant[n+2*sizeNode] = getQuantileSum(order, 0.01, tree->time[node[n]], tree->maxTime, savePext, iter);
	free((void*)savePext);
	if(fout != NULL && find != NULL) {
		for(s=0; s<iter; s++)
			fprintf(fout, "%d\t%lf\n",  s+burn+1, saveParam[s].birth);
		for(s=0; s<iter; s++)
			fprintf(fout, "%d\t%lf\n",  s+burn+1, saveParam[s].death);
		for(s=0; s<iter; s++)
			fprintf(fout, "%d\t%lf\n",  s+burn+1, saveParam[s].fossil);
		free((void*)saveParam);
		for(n=0; n<sizeNode; n++) {
			for(s=0; s<iter; s++)
				fprintf(fout, "%d\t%lf\n",  s+burn+1, quantTab[n][s]);
		}
		int off = 0;
		fprintf(find, "birth\t%d\t%d\n", off*iter+1, (off+1)*iter); off++;
		fprintf(find, "death\t%d\t%d\n", off*iter+1, (off+1)*iter); off++;
		fprintf(find, "fossil\t%d\t%d\n", off*iter+1, (off+1)*iter); off++;
		for(n=0; n<sizeNode; n++) {
			if(tree->name != NULL && tree->name[node[n]] != NULL)
				fprintf(find, "Quantile_%s\t%d\t%d\n", tree->name[node[n]], off*iter+1, (off+1)*iter);
			else
				fprintf(find, "Quantile_%d\t%d\t%d\n", node[n], off*iter+1, (off+1)*iter);
			off++;
		}
	}
	for(n=0; n<sizeNode; n++)
		free((void*)quantTab[n]);
	free((void*)quantTab);
	return quant;
}

double getLogProbExtinctCondSum(double t, double startTime, double maxTime, TypeExtendedModelParam *param, int size) {
	double sumLog = NEG_INFTY;
	int i;
	for(i=0; i<size; i++)
		sumLog = logSumLog(sumLog, getLogProbExtinctCond(t, startTime, maxTime, &(param[i])));
	return sumLog-log((double)size);
}

double getQuantileSum(double q, double tol, double startTime, double maxTime, TypeExtendedModelParam *param, int size) {
	double a, b;
	a = startTime;
	b = maxTime;
	while(b-a>tol) {
		double m = (a+b)/2.;
		if(exp(getLogProbExtinctCondSum(m, startTime, maxTime, param, size))>q)
			b = m;
		else
			a = m;
	}
	return (a+b)/2.;
}


double *MCMCSamplingExtinctionDistQuantMean(FILE *fout, FILE *find, double order, int **clade, TypeTree **tree, TypeFossilIntFeature *fi, double maxDisplayTime, double al, int burn, int gap, int iter, double prop, TypeModelParam *windSize, TypeModelParam *init, double probSpe, double probExt) {
	int *fol, *branch, i, j, f, s, n, nTree = 0, nClade = 0, nExt = 0, maxExt, *ext;
	double prob, **save, **saveExt, *quant;
	TypeFossilFeature *fp;
	TypeModelParam param, *saveParam;
	TypeExtendedModelParam pext, *savePext;
	
	param.birth = UNIF_RAND*init->birth;
	param.death = UNIF_RAND*init->death;
	param.fossil = UNIF_RAND*init->fossil;
	param.sampling = 1.;
	fp = sampleFossil(fi, tree[0]->size);
	for(nTree=0; tree[nTree] != NULL; nTree++) {
		for(n=0; n<tree[nTree]->size; n++) {
			if(tree[nTree]->node[n].child == NOSUCH) {
				switch(fi->status[n]) {
					case contempNodeStatus:
						tree[nTree]->time[n] = tree[nTree]->maxTime;
					break;
					case unknownNodeStatus:

					break;
					case extinctNodeStatus:
						tree[nTree]->time[n] = fp->fossilList[fp->fossil[n]].time;
					break;
					default:
	//						error("Node %d has no status\n", n);
	;
				}
			} else
				tree[nTree]->time[n] = NO_TIME;
		}
	}
	savePext = (TypeExtendedModelParam*) malloc(iter*sizeof(TypeExtendedModelParam));
	ext = (int*) malloc(tree[0]->size*sizeof(int));
	maxExt = 0;
	for(n=0; n<tree[0]->size; n++)
		if(tree[0]->node[n].child == NOSUCH && fi->status[n] == extinctNodeStatus) {
			ext[nExt++] = n;
			if(n>maxExt)
				maxExt = n;
		}
	saveExt = (double**) malloc((maxExt+1)*sizeof(double*));
	for(n=0; n<=maxExt; n++)
		saveExt[n] = (double*) malloc(iter*sizeof(double));
	fol = (int*) malloc(fp->size*sizeof(int));
	fillFolNode(fp, tree[0]->size, fol);
	branch = (int*) malloc(fp->size*sizeof(int));
	fillBranchNode(fp, tree[0]->size, branch);
	if(fout != NULL && find != NULL) {
		saveParam = (TypeModelParam*) malloc(iter*sizeof(TypeModelParam));
		save = (double**) malloc(fp->size*sizeof(double*));
		for(f=0; f<fp->size; f++)
			save[f] = (double*) malloc(iter*sizeof(double));
	}
	prob = getLogDensitySum(tree, nTree, fp, &param);
	for(i=0; i<burn; i++) {
		prob = update(tree, nTree, fp, fi, fol, branch, prob, al, prop, &param, windSize, probSpe, probExt);
if(i%10 == 0) {
	fprintf(stderr, "it %-7d\r", i); fflush(stderr);
}
	}
	for(nClade=0; clade[nClade]!=NULL; nClade++)
		;
	quant = (double*) malloc(nClade*sizeof(double));
	for(s=0; s<iter; s++) {
		pext =  getExtendedModelParam(&param);
		savePext[s] = pext;
		if(fout != NULL && find != NULL) {
			saveParam[s] = param;
			for(f=0; f<fp->size; f++)
				save[f][s] = fp->fossilList[f].time;
		}
		for(n=0; n<nExt; n++) 
			saveExt[ext[n]][s] = tree[0]->time[ext[n]];
		for(j=0; j<gap; j++) {
			prob = update(tree, nTree, fp, fi, fol, branch, prob, al, prop, &param, windSize, probSpe, probExt);
		}
if(s%10 == 0) {
	fprintf(stderr, "it %-7d\r", s); fflush(stderr);
}
	}
	for(n=0; n<nClade; n++) {
		int j;
		for(j=0; clade[n][j] != END_LIST_INT; j++)
			printf("clade[%d] = %d %s\n", j, clade[n][j], tree[0]->name[clade[n][j]]);
		quant[n] = getQuantileCladeSum(clade[n], order, 0.01, tree[0]->maxTime, savePext, saveExt, iter);
	}
	free((void*)savePext);
	if(fout != NULL && find != NULL) {
		for(s=0; s<iter; s++)
			fprintf(fout, "%d\t%lf\n",  s+burn+1, saveParam[s].birth);
		for(s=0; s<iter; s++)
			fprintf(fout, "%d\t%lf\n",  s+burn+1, saveParam[s].death);
		for(s=0; s<iter; s++)
			fprintf(fout, "%d\t%lf\n",  s+burn+1, saveParam[s].fossil);
		free((void*)saveParam);
		for(f=0; f<fp->size; f++)
			for(s=0; s<iter; s++)
				fprintf(fout, "%d\t%lf\n",  s+burn+1, save[f][s]);
		for(f=0; f<fp->size; f++)
			free((void*)save[f]);
		free((void*)save);
		for(n=0; n<nExt; n++)
			for(s=0; s<iter; s++)
				fprintf(fout, "%d\t%lf\n",  s+burn+1, saveExt[ext[n]][s]);
		int off = 0;
		fprintf(find, "birth\t%d\t%d\n", off*iter+1, (off+1)*iter); off++;
		fprintf(find, "death\t%d\t%d\n", off*iter+1, (off+1)*iter); off++;
		fprintf(find, "fossil\t%d\t%d\n", off*iter+1, (off+1)*iter); off++;
		for(f=0; f<fp->size; f++) {
			fprintf(find, "%d", f+1);
			if(branch[f] != NOSUCH && tree[0]->name[branch[f]] != NULL)
				fprintf(find, "_%s_%d_%d", tree[0]->name[branch[f]], (int) round(fi->fossilIntList[f].fossilInt.inf), (int) round(fi->fossilIntList[f].fossilInt.sup));
			fprintf(find, "\t%d\t%d\n", off*iter+1, (off+1)*iter); off++;
		}
		for(n=0; n<nExt; n++) {
			if(tree[0]->name[ext[n]] != NULL)
				fprintf(find, "%s", tree[0]->name[ext[n]]);
			else
				fprintf(find, "node_%d", ext[n]);
			fprintf(find, "\t%d\t%d\n", off*iter+1, (off+1)*iter); off++;
		}
	}
	for(n=0; n<=maxExt; n++)
		free((void*)saveExt[n]);
	free((void*)saveExt);
	free((void*)branch);
	free((void*)fol);
	return quant;
}


double getLogProbExtinctCondCladeSum(double val, int *clade, double maxTime, TypeExtendedModelParam *param, double **time, int size) {
	int i;
	double sumLog = NEG_INFTY;
	for(i=0; i<size; i++) {
		int j;
		double logProb = 0.;
		for(j=0; clade[j] != END_LIST_INT; j++)
			logProb += getLogProbExtinctCond(val, time[clade[j]][i], maxTime, &(param[i]));
		sumLog = logSumLog(sumLog, logProb);
	}
	return sumLog-log((double)size);
}

double getQuantileCladeSum(int *clade, double q, double tol, double maxTime, TypeExtendedModelParam *param, double **time, int size) {
	double a, b;
	int i, j;
	a = time[clade[0]][0];
	for(i=0; i<size; i++)
		for(j=0; clade[j] != END_LIST_INT; j++)
			if(time[clade[j]][i]<a)
				a = time[clade[j]][i];
	b = maxTime;
	while(b-a>tol) {
		double m = (a+b)/2.;
		if(exp(getLogProbExtinctCondCladeSum(m, clade, maxTime, param, time, size))>q)
			b = m;
		else
			a = m;
	}
	return (a+b)/2.;
}



double MCMCSamplingExtinctionComp(FILE *fout, FILE *find, int *cladeA, int *cladeB, TypeTree **tree, TypeFossilIntFeature *fi, double al, int burn, int gap, int iter, double prop, TypeModelParam *windSize, TypeModelParam *init, double probSpe, double probExt) {
	int *fol, *branch, i, j, f, s, n, nTree = 0, nExt = 0, maxExt, *ext;
	double prob, **save, **saveText, sumLog;
	TypeFossilFeature *fp;
	TypeModelParam param, *saveParam;
	TypeExtendedModelParam pext, *savePext;
	
	param.birth = UNIF_RAND*init->birth;
	param.death = UNIF_RAND*init->death;
	param.fossil = UNIF_RAND*init->fossil;
	param.sampling = 1.;
	fp = sampleFossil(fi, tree[0]->size);
	for(nTree=0; tree[nTree] != NULL; nTree++) {
		for(n=0; n<tree[nTree]->size; n++) {
			if(tree[nTree]->node[n].child == NOSUCH) {
				switch(fi->status[n]) {
					case contempNodeStatus:
						tree[nTree]->time[n] = tree[nTree]->maxTime;
					break;
					case unknownNodeStatus:

					break;
					case extinctNodeStatus:
						tree[nTree]->time[n] = fp->fossilList[fp->fossil[n]].time;
					break;
					default:
						;
				}
			} else
				tree[nTree]->time[n] = NO_TIME;
		}
	}
	savePext = (TypeExtendedModelParam*) malloc(iter*sizeof(TypeExtendedModelParam));
	ext = (int*) malloc(tree[0]->size*sizeof(int));
	maxExt = 0;
	for(n=0; n<tree[0]->size; n++)
		if(tree[0]->node[n].child == NOSUCH && fi->status[n] == extinctNodeStatus) {
			ext[nExt++] = n;
			if(n>maxExt)
				maxExt = n;
		}
	saveText = (double**) malloc(iter*sizeof(double*));
	for(s=0; s<iter; s++)
		saveText[s] = (double*) malloc((maxExt+1)*sizeof(double));
	fol = (int*) malloc(fp->size*sizeof(int));
	fillFolNode(fp, tree[0]->size, fol);
	branch = (int*) malloc(fp->size*sizeof(int));
	fillBranchNode(fp, tree[0]->size, branch);
	if(fout != NULL && find != NULL) {
		saveParam = (TypeModelParam*) malloc(iter*sizeof(TypeModelParam));
		save = (double**) malloc(fp->size*sizeof(double*));
		for(f=0; f<fp->size; f++)
			save[f] = (double*) malloc(iter*sizeof(double));
	}
	prob = getLogDensitySum(tree, nTree, fp, &param);
	for(i=0; i<burn; i++) {
		prob = update(tree, nTree, fp, fi, fol, branch, prob, al, prop, &param, windSize, probSpe, probExt);
if(i%10 == 0) {
	fprintf(stderr, "it %-7d\r", i); fflush(stderr);
}
	}
	sumLog = NEG_INFTY;
	for(s=0; s<iter; s++) {
		pext =  getExtendedModelParam(&param);
		savePext[s] = pext;
		if(fout != NULL && find != NULL) {
			saveParam[s] = param;
			for(f=0; f<fp->size; f++)
				save[f][s] = fp->fossilList[f].time;
		}
		for(n=0; n<nExt; n++) 
			saveText[s][ext[n]] = tree[0]->time[ext[n]];
		sumLog = logSumLog(sumLog, getLogProbComp(cladeA, cladeB, tree[0]->maxTime, &(pext), tree[0]->time));
//		sumLog = logSumLog(sumLog, log(getProbComp(cladeA, cladeB, tree[0]->time, tree[0]->maxTime, &(pext))));
		for(j=0; j<gap; j++) {
			prob = update(tree, nTree, fp, fi, fol, branch, prob, al, prop, &param, windSize, probSpe, probExt);
		}
if(s%10 == 0) {
	fprintf(stderr, "it %-7d\r", s); fflush(stderr);
}
	}
	free((void*)savePext);
	if(fout != NULL && find != NULL) {
		for(s=0; s<iter; s++)
			fprintf(fout, "%d\t%lf\n",  s+burn+1, saveParam[s].birth);
		for(s=0; s<iter; s++)
			fprintf(fout, "%d\t%lf\n",  s+burn+1, saveParam[s].death);
		for(s=0; s<iter; s++)
			fprintf(fout, "%d\t%lf\n",  s+burn+1, saveParam[s].fossil);
		free((void*)saveParam);
		for(f=0; f<fp->size; f++)
			for(s=0; s<iter; s++)
				fprintf(fout, "%d\t%lf\n",  s+burn+1, save[f][s]);
		for(f=0; f<fp->size; f++)
			free((void*)save[f]);
		free((void*)save);
		for(n=0; n<nExt; n++)
			for(s=0; s<iter; s++)
				fprintf(fout, "%d\t%lf\n",  s+burn+1, saveText[s][ext[n]]);
		int off = 0;
		fprintf(find, "birth\t%d\t%d\n", off*iter+1, (off+1)*iter); off++;
		fprintf(find, "death\t%d\t%d\n", off*iter+1, (off+1)*iter); off++;
		fprintf(find, "fossil\t%d\t%d\n", off*iter+1, (off+1)*iter); off++;
		for(f=0; f<fp->size; f++) {
			fprintf(find, "%d", f+1);
			if(branch[f] != NOSUCH && tree[0]->name[branch[f]] != NULL)
				fprintf(find, "_%s_%d_%d", tree[0]->name[branch[f]], (int) round(fi->fossilIntList[f].fossilInt.inf), (int) round(fi->fossilIntList[f].fossilInt.sup));
			fprintf(find, "\t%d\t%d\n", off*iter+1, (off+1)*iter); off++;
		}
		for(n=0; n<nExt; n++) {
			if(tree[0]->name[ext[n]] != NULL)
				fprintf(find, "%s", tree[0]->name[ext[n]]);
			else
				fprintf(find, "node_%d", ext[n]);
			fprintf(find, "\t%d\t%d\n", off*iter+1, (off+1)*iter); off++;
		}
	}
	free((void*)ext);
	freeFossilFeature(fp);
	for(s=0; s<iter; s++)
		free((void*)saveText[s]);
	free((void*)saveText);
	free((void*)branch);
	free((void*)fol);
	return sumLog-log((double)iter);
}

double sumExpSimpsonDistributionX(TypeDistribution d) {
	int i, n;
	double sum;
	if(d.size == 0)
		return 0;
	if(d.size == 1)
		return d.item[0].dens;
	n = d.size/2;
	sum = (d.item[0].dens>NEG_INFTY)?exp(d.item[0].dens):0.;
	for(i=1; i<n-1; i++)
		sum += 2.*((d.item[2*i].dens>NEG_INFTY)?exp(d.item[2*i].dens):0.);
	for(i=1; i<n; i++)
		sum += 4.*((d.item[2*i-1].dens>NEG_INFTY)?exp(d.item[2*i-1].dens):0.);
	sum += (d.item[2*n].dens>NEG_INFTY)?exp(d.item[2*n].dens):0.;
	sum *= (d.item[1].val-d.item[0].val)/3.;
	return sum;
}

double getLogProbComp(int *cladeA, int *cladeB, double maxTime, TypeExtendedModelParam *param, double *time) {
	double a, b, step = 0.1, result;
	int j;
	TypeDistribution d;

	a = time[cladeA[0]];
	for(j=1; cladeA[j] != END_LIST_INT; j++)
		if(time[cladeA[j]]>a)
			a = time[cladeA[j]];
	a += 0.00000001;
	b = time[cladeB[0]];
	for(j=1; cladeB[j] != END_LIST_INT; j++)
		if(time[cladeB[j]]>b)
			b = time[cladeB[j]];
	b += 0.00000001;
	d.size = (int) ceil((maxTime-a)/step);
	d.item = (TypeDistributionItem*) malloc(d.size*sizeof(TypeDistributionItem));
	for(j=0; j<d.size; j++) {
			d.item[j].val = a + ((double)j)*step;
			if(d.item[j].val>b)
				d.item[j].dens = exp(getLogDensClade(d.item[j].val, cladeA, time, maxTime, param)+getLogProbClade(d.item[j].val, cladeB, time, maxTime, param));
			else
				d.item[j].dens = 0.;
	}
	//FILE *fo;
	//if((fo = fopen("test_dist.csv", "w"))) {
		//fprintDistribution(fo, d);
		//fclose(fo);
	//}
	result = sumSimpsonDistribution(d);
	free((void*)d.item);
	return log(result);
}

double getProbComp(int *cladeA, int *cladeB, double *time, double tmax, TypeExtendedModelParam *pext) {
	int i, j, A, B, *uni;
	double maxA, maxB, max, res = 0., **tabD;
	
	for(A=0; cladeA[A] != END_LIST_INT; A++)
		;
	for(B=0; cladeB[B] != END_LIST_INT; B++)
		;
	maxA = time[cladeA[0]];
	for(i=1; cladeA[i] != END_LIST_INT; i++)
		if(time[cladeA[i]]>maxA)
			maxA = time[cladeA[i]];
	maxB = time[cladeB[0]];
	for(j=1; cladeB[j] != END_LIST_INT; j++)
		if(time[cladeB[j]]>maxB)
			maxB = time[cladeB[j]];
	if(maxA>maxB)
		max = maxA;
	else
		max = maxB;
	uni = (int*) malloc((A+B)*sizeof(int));
	for(i=0; i<A; i++)
		uni[i] = cladeA[i];
	for(j=0; j<B; j++)
		uni[j+A] = cladeB[j];
	tabD = (double**) malloc((A+B)*sizeof(double*));	
	for(i=0; i<A+B; i++) {
		int k;
		double tmp = 1.;
		tabD[i] = (double*) malloc(B*sizeof(double));
		for(k=0; k<i; k++)
			tmp *= (pext->beta-pext->alpha*exp(pext->omega*(time[uni[k]]-time[uni[i]])))/(1.-exp(pext->omega*(time[uni[k]]-time[uni[i]])));
		for(k=i+1; k<A+B; k++)
			tmp *= (pext->beta-pext->alpha*exp(pext->omega*(time[uni[k]]-time[uni[i]])))/(1.-exp(pext->omega*(time[uni[k]]-time[uni[i]])));
		for(j=0; j<B; j++)
			if(cladeB[j] != uni[i])
				tabD[i][j] = (tmp*(1.-exp(pext->omega*(time[cladeB[j]]-time[uni[i]]))))/(pext->beta-pext->alpha*exp(pext->omega*(time[cladeB[j]]-time[uni[i]])));
	}
	for(j=0; j<B; j++) {
		double termX, termY, tmp = 0.;
		for(i=0; i<A+B; i++)
			if(uni[i] != cladeB[j])
				tmp += (tabD[i][j])/(1.-exp(pext->omega*(time[cladeB[j]]-time[uni[i]])));
		termX = (pext->beta*(pext->bma*tmp-pow(pext->beta, A+B-1.)))*((1/(pext->beta-pext->alpha*exp(pext->omega*(tmax-time[cladeB[j]]))))-(1/(pext->beta-pext->alpha*exp(pext->omega*(max-time[cladeB[j]])))));
		tmp = 0.;
		for(i=0; i<A+B; i++)
			if(uni[i] != cladeB[j])
				tmp += (tabD[i][j]*exp(pext->omega*(time[cladeB[j]]-time[uni[i]]))*(log1p(pext->alpha*(1-exp(pext->omega*(time[cladeB[j]]-time[uni[i]])))/(pext->beta*exp(pext->omega*(time[cladeB[j]]-max))-pext->alpha))-log1p(pext->alpha*(1-exp(pext->omega*(time[cladeB[j]]-time[uni[i]])))/(pext->beta*exp(pext->omega*(time[cladeB[j]]-tmax))-pext->alpha))))/pow(1.-exp(pext->omega*(time[cladeB[j]]-time[uni[i]])), 2.);
		termY = pext->bma*tmp;
		res += termX+termY;
printf("termX %lf %lf\n", termX, termY);
	}
	for(i=0; i<A+B; i++)
		free((void*) tabD[i]);
	free((void*) tabD);
	free((void*) uni);
	res *= pext->bma;
printf("res %.2lf\n", res/exp(getLogProbExtinctClade(cladeA, time, tmax, pext)+getLogProbExtinctClade(cladeB, time, tmax, pext)));
	return res/exp(getLogProbExtinctClade(cladeA, time, tmax, pext)+getLogProbExtinctClade(cladeB, time, tmax, pext));
}

double MCMCSamplingSingleTreeFossilDensity(TypeTree *tree, TypeFossilIntFeature *fi, double al, int burn, int gap, int iter, double prop, TypeModelParam *windSize, TypeModelParam *init, double probSpe, double probExt) {
	int *fol, *branch, i, s, n;
	double prob, sumLog = NEG_INFTY;
	TypeFossilFeature *fp;
	TypeModelParam param;
	
	param.birth = UNIF_RAND*init->birth;
	param.death = UNIF_RAND*init->death;
	param.fossil = UNIF_RAND*init->fossil;
	param.sampling = 1.;


	fp = sampleFossil(fi, tree->size);
	for(n=0; n<tree->size; n++) {
		if(tree->node[n].child == NOSUCH) {
			switch(fi->status[n]) {
				case contempNodeStatus:
					tree->time[n] = tree->maxTime;
				break;
				case unknownNodeStatus:

				break;
				case extinctNodeStatus:
					tree->time[n] = fp->fossilList[fp->fossil[n]].time;
				break;
				default:
//						error("Node %d has no status\n", n);
;
			}
		} else
			tree->time[n] = NO_TIME;
	}
	fol = (int*) malloc(fp->size*sizeof(int));
	fillFolNode(fp, tree->size, fol);
	branch = (int*) malloc(fp->size*sizeof(int));
	fillBranchNode(fp, tree->size, branch);
	prob = getLogDensitySumNoThread(&tree, 1, fp, &param);

	for(i=0; i<burn; i++) {
		prob = update(&tree, 1, fp, fi, fol, branch, prob, al, prop, &param, windSize, probSpe, probExt);
if(i%10 == 0) {
	fprintf(stderr, "it %-7d\r", i); fflush(stderr);
}
	}
	for(s=0; s<iter; s++) {
		int j;
		sumLog = logSumLog(sumLog, prob);
		for(j=0; j<gap; j++) {
			prob = update(&tree, 1, fp, fi, fol, branch, prob, al, prop, &param, windSize, probSpe, probExt);
		}
if(s%10 == 0) {
	fprintf(stderr, "it %-7d\r", s); fflush(stderr);
}
	}
	freeFossilFeature(fp);
	free((void*)branch);
	free((void*)fol);
	return sumLog-log((double) iter);
}



