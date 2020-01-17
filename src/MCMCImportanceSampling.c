#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <pthread.h>
#include <signal.h>
#include <nlopt.h>

#include "Utils.h"
#include "MyError.h"
#include "Uncertainty.h"
#include "MCMCImportanceSampling.h"


typedef struct MINIMIZATION_PARAM_DATA {
    TypeTree **tree;
    TypeFossilFeature **fos;
    int nTree;
} TypeMinimizationParamData;


typedef struct THREAD_PARAMETER_DENS {
	int *number;
	pthread_mutex_t *mutex_number;
	pthread_cond_t *cond_number;
	TypeTree *tree;
	TypeFossilFeature *ffe;
	int node;
	double time, *densTot, *densTime;
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


static pthread_mutex_t mutexN = PTHREAD_MUTEX_INITIALIZER;
static pthread_cond_t condN = PTHREAD_COND_INITIALIZER;
static int maxT = 40;

static void *threadComputeDens(void *data);
static void *threadComputeDistribution(void *data);

static void fillFolNode(TypeFossilFeature *f, int size, int *fol);
static void fillBranchNode(TypeFossilFeature *f, int size, int *branch);
static void  getProposalFossil(TypeFossilFeature *fp, TypeFossilIntFeature *fi, int *fol, double al, int *move, double *newTime);
static void changeFossil(TypeTree **tree, int nTree, TypeFossilFeature **fp, int *fol, int **branch, int move, double newTime);
static double update(TypeTree **tree, int nTree, TypeFossilFeature **fp, TypeFossilIntFeature **fi, int *fol, int **branch, double prob, double al, double propParam, TypeModelParam *param, TypeModelParam *windSize, double probSpe, double probExt);
static TypeFossilFeature *sampleFossil(TypeFossilIntFeature* feat, int size);
static int compareLocal(const void* a, const void* b);
static double getLogDensitySum(TypeTree **tree, int nTree, TypeFossilFeature **fp, TypeModelParam *param);

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
//printf("fossil %d t %.2lf [%.2lf, %.2lf]\n", i, sample->fossilList[i].time, feat->fossilIntList[i].fossilInt.inf, feat->fossilIntList[i].fossilInt.sup);
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

void  getProposalFossil(TypeFossilFeature *fp, TypeFossilIntFeature *fi, int *fol, double al, int *move, double *newTime) {
	double length;
    *move = RANGE_RAND(fp->size);
    length = (fi->fossilIntList[*move].fossilInt.sup-fi->fossilIntList[*move].fossilInt.inf)*al/2.;
    *newTime = UNIF_RAND*(2.*length)+fp->fossilList[*move].time-length;
    if(*newTime<fi->fossilIntList[*move].fossilInt.inf)
		*newTime = 2.*fi->fossilIntList[*move].fossilInt.inf - *newTime;
    if(*newTime>fi->fossilIntList[*move].fossilInt.sup)
		*newTime = 2.*fi->fossilIntList[*move].fossilInt.sup - *newTime;
}

void changeFossil(TypeTree **tree, int nTree, TypeFossilFeature **fp, int *fol, int **branch, int move, double newTime) {
	int xA, xB, change, i;
	change = (tree[0]->time[branch[0][move]] == fp[0]->fossilList[fp[0]->fossil[branch[0][move]]].time);
	fp[0]->fossilList[move].time = newTime;
	xA = NOSUCH;
	for(xB=fol[move]; xB!=NOSUCH && fp[0]->fossilList[move].time>fp[0]->fossilList[xB].time; xB=fol[xB])
		xA = xB;
	if(xB != fol[move]) {;
		if(fol[move] != NOSUCH)
			fp[0]->fossilList[fol[move]].prec = fp[0]->fossilList[move].prec;
		else
			fp[0]->fossil[branch[0][move]] = fp[0]->fossilList[move].prec;
		if(fp[0]->fossilList[move].prec != NOSUCH)
			fol[fp[0]->fossilList[move].prec] = fol[move];
		if(xB != NOSUCH)
			fp[0]->fossilList[xB].prec = move;
		else
			for(i=0; i<nTree; i++)
				fp[i]->fossil[branch[i][move]] = move;
		if(xA != NOSUCH)
			fol[xA] = move;
		fol[move] = xB;
		fp[0]->fossilList[move].prec = xA;
	}
	xA = NOSUCH;
	for(xB=fp[0]->fossilList[move].prec; xB!=NOSUCH && fp[0]->fossilList[move].time<fp[0]->fossilList[xB].time; xB=fp[0]->fossilList[xB].prec)
		xA = xB;
	if(xB != fp[0]->fossilList[move].prec) {
		if(fol[move] != NOSUCH)
			fp[0]->fossilList[fol[move]].prec = fp[0]->fossilList[move].prec;
		else
			for(i=0; i<nTree; i++)
				fp[i]->fossil[branch[i][move]] = fp[0]->fossilList[move].prec;
		if(fp[0]->fossilList[move].prec != NOSUCH)
			fol[fp[0]->fossilList[move].prec] = fol[move];
		if(xB != NOSUCH)
			fol[xB] = move;
		if(xA != NOSUCH)
			fp[0]->fossilList[xA].prec = move;
		fol[move] = xA;
		fp[0]->fossilList[move].prec = xB;
	}
	double min = getMinFossilTime(fp[0]);
	if(tree[0]->minTimeInt.sup > min)
		for(i=0; i<nTree; i++)
			tree[i]->minTimeInt.sup = min;
	if(change)
		for(i=0; i<nTree; i++)
			tree[i]->time[branch[i][move]] = fp[i]->fossilList[fp[i]->fossil[branch[i][move]]].time;
}

double update(TypeTree **tree, int nTree, TypeFossilFeature **fp, TypeFossilIntFeature **fi, int *fol, int **branch, double prob, double al, double propParam, TypeModelParam *param, TypeModelParam *windSize, double probSpe, double probExt) {
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
		getProposalFossil(fp[0], fi[0], fol, al, &move, &newTime);
		oldTime = fp[0]->fossilList[move].time;
		changeFossil(tree, nTree, fp, fol, branch, move, newTime);
		newProb = getLogDensitySum(tree, nTree, fp, param);
		if(isnan(newProb))
			error("fossil %d %.5lf -> %.5lf (%d %s) time %.5lf\n", move, oldTime, newTime, branch[0][move], tree[0]->name[branch[0][move]], tree[0]->time[branch[0][move]]);
		if(UNIF_RAND < exp(newProb-prob)) {
			return newProb;
		} else {
			changeFossil(tree, nTree, fp, fol, branch, move, oldTime);
			return prob;
		}
	}
	return prob;
}

double getLogDensitySum(TypeTree **tree, int nTree, TypeFossilFeature **fp, TypeModelParam *param) {
	double sum, *dens;
	int cont=1, nT=0, i = 0;
	cont = 1;
	dens = (double*) malloc(nTree*sizeof(double));
	for(i=0; i<nTree; i++)
		dens[i] = 0;
	i = 0;
	while(cont) {
		pthread_mutex_lock(&mutexN);
		while(i < nTree && nT < maxT) {
			pthread_t thread;	
			int ret = 0;
			TypeThreadParameterDens *tp;
			tp = (TypeThreadParameterDens*) malloc(sizeof(TypeThreadParameterDens));
			tp->tree = tree[i];
			tp->ffe = fp[i];
			tp->param = param;
			tp->densTime = &(dens[i]);
			tp->mutex_number = &mutexN;
			tp->cond_number = &condN;
			tp->number = &nT;
			if((ret = pthread_create(&thread, NULL, threadComputeDens, (void*) tp)) == 0) {
				int err;
				if((err = pthread_detach(thread)) == 0) {
					nT++; i++;
				} else {
					fprintf (stderr, "Error %d while detaching thread: %s\n", err, (char*) strerror(err));
				}
			} else
				fprintf (stderr, "Error %d while creating thread: %s\n", ret, (char*) strerror(ret));
		}
		cont = (nT > 0);
		if(cont)
			pthread_cond_wait(&condN, &mutexN);
		pthread_mutex_unlock(&mutexN);
	}
	sum = NEG_INFTY;
	for(i=0; i<nTree; i++)
		sum = logSumLog(sum, dens[i]);
	free((void*)dens);
	return sum;
}

void *threadComputeDens(void *data) {
	*(((TypeThreadParameterDens*)data)->densTime) = getLogLikelihoodTreeFossil(((TypeThreadParameterDens*)data)->tree, ((TypeThreadParameterDens*)data)->ffe, ((TypeThreadParameterDens*)data)->param);
	pthread_mutex_lock(((TypeThreadParameterDens*)data)->mutex_number);
		(*((TypeThreadParameterDens*)data)->number)--;
		pthread_cond_signal(((TypeThreadParameterDens*)data)->cond_number);
	pthread_mutex_unlock(((TypeThreadParameterDens*)data)->mutex_number);
	free(data);
	return NULL;
}

TypeDistribution MCMCSamplingDist(FILE *fout, FILE *find, TypeTree **tree, int nTree, int *node, TypeFossilIntFeature **fi, double step, double al, int burn, int gap, int iter, double prop, TypeModelParam *windSize, TypeModelParam *init, double probSpe, double probExt) {
	int *fol, **branch, i, j, s, n,  **indexNode;
	double prob, min, max, *logCond, logSumCond;
	TypeFossilFeature **fp;
	TypeDistribution *logDTmp, logSTmp, dist;
	TypeLexiTree *dict;
	TypeModelParam param, *saveParam;
	
	param.birth = UNIF_RAND*init->birth;
	param.death = UNIF_RAND*init->death;
	param.fossil = UNIF_RAND*init->fossil;
	param.sampling = 1.;

	//param.birth = 0.5;
	//param.death = 0.45;
	//param.fossil = 1.25;
	//param.sampling = 1.;


	dict = newLexiTree();
	for(n=0; n<tree[0]->size; n++)
		if(tree[0]->name && tree[0]->name[n]) {
			if(addWordLexi(tree[0]->name[n], n, dict)>=0)
				warning("Warning! duplicate identifier '%s'\n", tree[0]->name[n]);
		}
	indexNode = (int**) malloc(nTree*sizeof(int*));
	for(i=0; i<nTree; i++) {
		indexNode[i] = (int*) malloc(tree[i]->size*sizeof(int));
		for(n=0; n<tree[i]->size; n++)
			if(tree[i]->name && tree[i]->name[n])
				indexNode[i][n] = findWordLexi(tree[i]->name[n], dict);
			else
				indexNode[i][n] = NOSUCH;
	}
	freeLexiTree(dict);
	fp = (TypeFossilFeature**) malloc(nTree*sizeof(TypeFossilFeature*));
	fp[0] = sampleFossil(fi[0], tree[0]->size);
	for(i=1; i<nTree; i++) {
		fp[i] = (TypeFossilFeature*) malloc(sizeof(TypeFossilFeature));
		fp[i]->size = fp[0]->size;
		fp[i]->sizeBuf = fp[0]->sizeBuf;
		fp[i]->fossilList = fp[0]->fossilList;
		fp[i]->fossil = (int*) malloc(tree[i]->size*sizeof(int));
		fp[i]->status = NULL;
		for(n=0; n<tree[i]->size; n++) {
			if(indexNode[i][n] != NOSUCH)
				fp[i]->fossil[n] = fp[0]->fossil[indexNode[i][n]];
			else
				fp[i]->fossil[n] = NOSUCH;
		}
	}
	for(i=0; i<nTree; i++) {
		for(n=0; n<tree[i]->size; n++) {
			if(tree[i]->node[n].child == NOSUCH) {
				switch(fi[i]->status[n]) {
					case contempNodeStatus:
						tree[i]->time[n] = tree[i]->maxTime;
					break;
					case unknownNodeStatus:

					break;
					case extinctNodeStatus:
						tree[i]->time[n] = fp[i]->fossilList[fp[i]->fossil[n]].time;
					break;
					default:
	//						error("Node %d has no status\n", n);
	;
				}
			} else
				tree[i]->time[n] = NO_TIME;
		}
	}
	fol = (int*) malloc(fp[0]->size*sizeof(int));
	fillFolNode(fp[0], tree[0]->size, fol);
	branch = (int**) malloc(nTree*sizeof(int*));
	for(i=0; i<nTree; i++) {
		branch[i] = (int*) malloc(fp[i]->size*sizeof(int));
		fillBranchNode(fp[i], tree[i]->size, branch[i]);
	}
	min = getMaxFossilIntTimeToNode(node[0], tree[0], fi[0]);
	max = getMinFossilIntTimeFromNode(node[0], tree[0], fi[0]);
	min = floor(min/step)*step;
	max = (ceil(max/step))*step;
	double **save;
	int f;
	saveParam = (TypeModelParam*) malloc(iter*sizeof(TypeModelParam));
	save = (double**) malloc(fp[0]->size*sizeof(double*));
	for(f=0; f<fp[0]->size; f++)
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
				tp->ffe = fp[i];
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
						fprintf (stderr, "Error %d while detaching thread: %s\n", err, (char*) strerror(err));
					}
				} else
					fprintf (stderr, "Error %d while creating thread: %s\n", ret, (char*) strerror(ret));
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
		for(f=0; f<fp[0]->size; f++)
			save[f][s] = fp[0]->fossilList[f].time;
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
	for(f=0; f<fp[0]->size; f++)
		for(s=0; s<iter; s++)
			fprintf(fout, "%d\t%lf\n",  s+burn+1, save[f][s]);
	for(f=0; f<fp[0]->size; f++)
		free((void*)save[f]);
	free((void*)save);
	int off = 0;
	fprintf(find, "birth\t%d\t%d\n", off*iter+1, (off+1)*iter); off++;
	fprintf(find, "death\t%d\t%d\n", off*iter+1, (off+1)*iter); off++;
	fprintf(find, "fossil\t%d\t%d\n", off*iter+1, (off+1)*iter); off++;
	for(f=0; f<fp[0]->size; f++) {
		fprintf(find, "%d", f+1);
		if(branch[0][f] != NOSUCH && tree[0]->name[branch[0][f]] != NULL)
			fprintf(find, "_%s_%d_%d", tree[0]->name[branch[0][f]], (int) round(fi[0]->fossilIntList[f].fossilInt.inf), (int) round(fi[0]->fossilIntList[f].fossilInt.sup));
		fprintf(find, "\t%d\t%d\n", off*iter+1, (off+1)*iter); off++;
	}
	for(i=0; i<nTree; i++) {
		free((void*)indexNode[i]);
		free((void*)branch[i]);
	}
	free((void*)indexNode);
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


void MCMCSamplingDiag(FILE *fout, FILE *find, TypeTree **tree, int nTree, TypeFossilIntFeature **fi, double al, int burn, int gap, int iter, double prop, TypeModelParam *windSize, TypeModelParam *init, double probSpe, double probExt) {
	int *fol, **branch, f, i , s, n,  **indexNode;
	double prob, **save;
	TypeFossilFeature **fp;
	TypeLexiTree *dict;
	TypeModelParam param, *saveParam;
	
	param.birth = UNIF_RAND*init->birth;
	param.death = UNIF_RAND*init->death;
	param.fossil = UNIF_RAND*init->fossil;
	param.sampling = 1.;
	
	dict = newLexiTree();
	for(n=0; n<tree[0]->size; n++)
		if(tree[0]->name && tree[0]->name[n]) {
			if(addWordLexi(tree[0]->name[n], n, dict)>=0)
				warning("Warning! duplicate identifier '%s'\n", tree[0]->name[n]);
		}
	indexNode = (int**) malloc(nTree*sizeof(int*));
	for(i=0; i<nTree; i++) {
		indexNode[i] = (int*) malloc(tree[i]->size*sizeof(int));
		for(n=0; n<tree[i]->size; n++)
			if(tree[i]->name && tree[i]->name[n])
				indexNode[i][n] = findWordLexi(tree[i]->name[n], dict);
			else
				indexNode[i][n] = NOSUCH;
	}
	freeLexiTree(dict);
	fp = (TypeFossilFeature**) malloc(nTree*sizeof(TypeFossilFeature*));
	fp[0] = sampleFossil(fi[0], tree[0]->size);
	for(i=1; i<nTree; i++) {
		fp[i] = (TypeFossilFeature*) malloc(sizeof(TypeFossilFeature));
		fp[i]->size = fp[0]->size;
		fp[i]->sizeBuf = fp[0]->sizeBuf;
		fp[i]->fossilList = fp[0]->fossilList;
		fp[i]->fossil = (int*) malloc(tree[i]->size*sizeof(int));
//		fp[i]->status = (TypeNodeStatus*) malloc(tree[i]->size*sizeof(int));
		fp[i]->status = NULL;
		for(n=0; n<tree[i]->size; n++) {
			if(indexNode[i][n] != NOSUCH)
				fp[i]->fossil[n] = fp[0]->fossil[indexNode[i][n]];
			else
				fp[i]->fossil[n] = NOSUCH;
		}
	}
	for(i=0; i<nTree; i++) {
		for(n=0; n<tree[i]->size; n++) {
			if(tree[i]->node[n].child == NOSUCH) {
				switch(fi[i]->status[n]) {
					case contempNodeStatus:
						tree[i]->time[n] = tree[i]->maxTime;
					break;
					case unknownNodeStatus:

					break;
					case extinctNodeStatus:
						tree[i]->time[n] = fp[i]->fossilList[fp[i]->fossil[n]].time;
					break;
					default:
	//						error("Node %d has no status\n", n);
	;
				}
			} else
				tree[i]->time[n] = NO_TIME;
		}
	}
	fol = (int*) malloc(fp[0]->size*sizeof(int));
	fillFolNode(fp[0], tree[0]->size, fol);
	branch = (int**) malloc(nTree*sizeof(int*));
	for(i=0; i<nTree; i++) {
		branch[i] = (int*) malloc(fp[i]->size*sizeof(int));
		fillBranchNode(fp[i], tree[i]->size, branch[i]);
	}
	saveParam = (TypeModelParam*) malloc(iter*sizeof(TypeModelParam));
	save = (double**) malloc(fp[0]->size*sizeof(double*));
	for(f=0; f<fp[0]->size; f++)
		save[f] = (double*) malloc(iter*sizeof(double));
	prob = getLogDensitySum(tree, nTree, fp, &param);
	for(i=0; i<burn; i++) {
		prob = update(tree, nTree, fp, fi, fol, branch, prob, al, prop, &param, windSize, probSpe, probExt);
if(i%10 == 0) {
	fprintf(stderr, "it %d\r", i); fflush(stderr);
}
	}
//fprintf(stderr, "acceptance ratio %lf\n", 100.*((double)countChange)/((double)burn));
printf("end burn\n");
	for(s=0; s<iter; s++) {
		int j, f;
		saveParam[s] = param;
		for(f=0; f<fp[0]->size; f++)
			save[f][s] = fp[0]->fossilList[f].time;
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
	for(f=0; f<fp[0]->size; f++)
		for(s=0; s<iter; s++)
			fprintf(fout, "%d\t%lf\n",  s+burn+1, save[f][s]);
	for(f=0; f<fp[0]->size; f++)
		free((void*)save[f]);
	free((void*)save);
	int off = 0;
	fprintf(find, "birth\t%d\t%d\n", off*iter+1, (off+1)*iter); off++;
	fprintf(find, "death\t%d\t%d\n", off*iter+1, (off+1)*iter); off++;
	fprintf(find, "fossil\t%d\t%d\n", off*iter+1, (off+1)*iter); off++;
	for(f=0; f<fp[0]->size; f++) {
		fprintf(find, "%d", f+1);
		if(branch[0][f] != NOSUCH && tree[0]->name[branch[0][f]] != NULL)
			fprintf(find, "_%s_%d_%d", tree[0]->name[branch[0][f]], (int) round(fi[0]->fossilIntList[f].fossilInt.inf), (int) round(fi[0]->fossilIntList[f].fossilInt.sup));
		fprintf(find, "\t%d\t%d\n", off*iter+1, (off+1)*iter); off++;
	}
	for(i=0; i<nTree; i++) {
		free((void*)indexNode[i]);
		free((void*)branch[i]);
	}
	free((void*)indexNode);
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

	prob = getLogDensitySum(&tree, 1, &fp, &param);
	for(i=0; i<burn; i++) {
		prob = update(&tree, 1, &fp, &fi, fol, &branch, prob, al, prop, &param, windSize, probSpe, probExt);
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
						fprintf (stderr, "Error %d while detaching thread: %s\n", err, (char*) strerror(err));
					}
				} else
					fprintf (stderr, "Error %d while creating thread: %s\n", ret, (char*) strerror(ret));
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
			prob = update(&tree, 1, &fp, &fi, fol, &branch, prob, al, prop, &param, windSize, probSpe, probExt);
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
