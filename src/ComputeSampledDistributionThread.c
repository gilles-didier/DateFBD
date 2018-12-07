#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <pthread.h>
#include <signal.h>

#include <gsl/gsl_rng.h>

#include "Utils.h"
#include "Tree.h"
#include "SimulTree.h"
#include "Fossil.h"
#include "FossilInt.h"
#include "Model.h"
#include "Uncertainty.h"
#include "SimulFossil.h"
#include "MinimizeNLOpt.h"
#include "Distribution.h"

#ifdef DO_PS
#endif

#define STRING_SIZE 300
#define INC_CUT 5
#define SEQ_SIZE 30
#define EXT_OUTPUT "_added.phy"
#define MAX_PRED 7.5E+07

#define SIZE_BUFFER_CHAR 300
#define INFTY 1E99
#define RINFTY 1E99
#define DEF 10
#define MIN_VAL 0.000001
#define DELTA 0.000001

#define MINVAL 0.01
#define TRIAL 10
#define FIX_VAL(x) (((x)<=0)?MINVAL:(x))

#define MAX_ITER 1000

#define M_MAX 6
#define M_MAX_F 4
#define MIN_TREE 20
#define PREFIX "table"
#define MAX_TRIALS 1000

#define HELPMESSAGE "--------------------------\n\nNAME\n	dist - Computation of the divergence time distibution associated to a given clade from a set of possible trees, the fossil ages, and the diversification rates\n	\nSYNOPSIS\n	dist [OPTIONS] <input Tree File> <input Fossil File> <input List Clade> [output File]\n\nDESCRIPTION\n	Compute the distribution of the divergence time associated to the clade corresponding to the list of tips given in <input List Clade> by sampling into the trees contianed in <input Tree File> (in Newick format) into the fossil ranges provided in <input Fossil File> (in csv format) and output the distribution as a .csv table <output File>.csv\n\n	Options are\n	-z <input Tree File>\n		output the tree in text debug format in the console and exit \n	-o <origin bound inf> [origin bound sup]\n		set the origin range; if origin bound sup is not given it is set to the most ancient fossil age\n	-e <end bound inf> <end bound sup>\n		set the end time range\n	-p <speciation rate> <extinction rate> <fossilization rate>\n		set the speciation, extinction and fossilization rates\n	-s <number>\n		set the number of samples\n	-d\n		return the distribution (otherwise the density is returned by default)\n	-u <value>\n		set the step discretizing the time distribution to <value>\n	-s <number>\n		set the number of thread running in parallell\n	-h\n		display help\n\n--------------------------\n\n"



typedef struct THREAD_PARAMETER {
	int *number;
	pthread_mutex_t *mutex_number;
	pthread_cond_t *cond_number;
	TypeTree *tree;
	TypeFossilFeature *ffe;
	int def, node;
	double min, max, step;
	TypeDistribution *logD;
	double *logCond;
	TypeModelParam *param;
} TypeThreadParameter;


static char **readList(FILE *f);
static void *threadComputeDistribution(void *data);

int main(int argc, char **argv) {	
	char *inputFileNameTree, *inputFileNameList, *inputFileNameFossil = NULL, *outputPrefix = PREFIX, outputDistribution[STRING_SIZE], option[256];
	FILE *fi, *fl, *ff, *fo;
	int i, j, nSamp = 100, def, outDens = 1, maxT = 2;
	double minTimeIntInf = NO_TIME, minTimeIntSup = NO_TIME, maxTimeIntInf = 0., maxTimeIntSup = 0., step = 0.1;
	TypeModelParam param = {.birth=1.3, .death = 1., .fossil = 1., .sampling = 1.};
	pthread_mutex_t mutexN = PTHREAD_MUTEX_INITIALIZER;
	pthread_cond_t condN = PTHREAD_COND_INITIALIZER;

	
	for(i=0; i<256; i++)
		option[i] = 0;
	   
	for(i=1; i<argc && *(argv[i]) == '-'; i++) {
		for(j=1; argv[i][j] != '\0'; j++)
			option[argv[i][j]] = 1;
		if(option['z']) {
			option['z'] = 0;
			if((i+1)<argc && (fi = fopen(argv[++i], "r"))) {
				TypeTree *tree;
				tree = readTree(fi);
				toBinary(tree);
				printTreeDebug(stdout, tree->root, tree, tree->name);
				reorderTree(tree->name, tree);
				if(tree->minTime == NO_TIME || tree->minTime == 0.)
					tree->minTime = tree->time[tree->root]*0.9;
				if(tree->maxTime == NO_TIME) {
					int n;
					tree->maxTime = 0.;
					for(n=0; n<tree->size; n++)
						if(tree->time[n]>tree->maxTime)
							tree->maxTime = tree->time[n];
				}
				printTreeDebug(stdout, tree->root, tree, tree->name);
			} else {
				fprintf(stderr, "Error while reading %s.\n", argv[i]);
				exit(1);
			}
			exit(0);
		}
		if(option['o']) {
			option['o'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%le", &minTimeIntInf) == 1)
				i++;
			else
				exitProg(ErrorArgument, "at least one value is expected after -o");
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &minTimeIntSup) == 1)
				i++;
			else
				minTimeIntSup = NO_TIME;
		}
		if(option['e']) {
			option['e'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%le", &maxTimeIntInf) == 1)
				i++;
			else
				exitProg(ErrorArgument, "at least one value is expected after -o");
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &maxTimeIntSup) == 1)
				i++;
			else
				maxTimeIntSup = maxTimeIntInf;
		}
		if(option['p']) {
			option['p'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%le", &(param.birth)) == 1)
				i++;
			else
				exitProg(ErrorArgument, "3 values are expected after -p");
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &(param.death)) == 1)
				i++;
			else
				exitProg(ErrorArgument, "3 values are expected after -p");
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &(param.fossil)) == 1)
				i++;
			else
				exitProg(ErrorArgument, "3 values are expected after -p");

		}
		if(option['s']) {
			option['s'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%d", &nSamp) == 1)
				i++;
			else
				exitProg(ErrorArgument, "a number is expected after -s");
		}
		if(option['d']) {
			option['d'] = 0;
			outDens = 0;
		}
		if(option['u']) {
			option['u'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%le", &(step)) == 1)
				i++;
			else
				exitProg(ErrorArgument, "a number is expected after -u");
		}
		if(option['t']) {
			option['t'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%d", &maxT) == 1)
				i++;
			else
				exitProg(ErrorArgument, "a number is required after option -t");
		}
		if(option['h']) {
			option['h'] = 0;
			printf("%s\n", HELPMESSAGE);
			exit(0);
		}
	}
	if(i<argc) {
		inputFileNameTree = argv[i++];
	} else {
		fprintf(stderr, "Please provide the name of a file containing a phylogenetic tree in Newick format\n");
		exit(1);
	}
	if(i<argc) {
		inputFileNameFossil = argv[i++];
	} else {
		fprintf(stderr, "Please provide the name of a file containing a fossil list\n");
		exit(1);
	}
	if(i<argc) {
		inputFileNameList = argv[i++];
	} else {
		fprintf(stderr, "Please provide the name of a file containing a list of leaf names of a clade\n");
		exit(1);
	}
	if(i<argc)
		outputPrefix = argv[i++];
	if((fi = fopen(inputFileNameTree, "r")) && (fl = fopen(inputFileNameList, "r")) && (ff = fopen(inputFileNameFossil, "r"))) {
		TypeTree **tree;
		char **list;
		TypeFossilIntFeature **fos;
		TypeDistribution *logD, dist;
		double *logCond, min, max;
		int i, n, sizeTree, s, *node, nT, cont;
		gsl_rng *rg = gsl_rng_alloc(gsl_rng_random_glibc2);
		
		list = readList(fl);
		fclose(fl);
        tree = readTrees(fi);
        fclose(fi);
        if(tree[0] == NULL) {
			fprintf(stderr, "Error: no tree\n");
			return 1;
		}
		for(i=0; tree[i]!=NULL; i++)
			;
		fos = (TypeFossilIntFeature**) malloc(i*sizeof(TypeFossilIntFeature*));
		node = (int*) malloc(i*sizeof(int));
		sizeTree = 0;
		for(i=0; tree[i]!=NULL; i++) {
			int n;
			toBinary(tree[i]);
			if(tree[i]->name!=NULL)
				for(n=0; n<tree[i]->size; n++)
					if(tree[i]->name[n]!=NULL)
						fixSpace(tree[i]->name[n]);
			n = getClade(list, tree[i]);
			if(n != NOSUCH) {
				node[sizeTree] = n;
				tree[sizeTree] = tree[i];
				rewind(ff);
				fos[sizeTree] = getFossilIntFeature(ff, tree[sizeTree]->name, tree[sizeTree]->size);
				if(getMaxFossilIntTime(fos[sizeTree]) > 0.)
					negateFossilInt(fos[sizeTree]);
				fixStatus(tree[sizeTree], fos[sizeTree]);
				sizeTree++;
			}
		}
		if(sizeTree == 0) {
			fprintf(stderr, "Error: no tree (size)\n");
			exit(0);
		}
		logD = (TypeDistribution*) malloc(nSamp*sizeof(TypeDistribution));
		logCond = (double*) malloc(nSamp*sizeof(double));
		double minFossilTime = getMinFossilIntTime(fos[0]);
		if(minTimeIntInf == NO_TIME || minTimeIntInf > minFossilTime) {
			if(minFossilTime<0)
				minTimeIntInf = 1.2*minFossilTime;
			else
				minTimeIntInf = 0.8*minFossilTime;
		}
		min = minTimeIntInf;
		max = getMinFossilIntTimeFromNode(node[0], tree[0], fos[0]);
		min = floor(min/step)*step;
		max = ceil(max/step)*step;
		def = 1+round((max-min)/step);
		s = 0;
		nT = 0;
		cont = 1;
		while(cont) {
			pthread_mutex_lock(&mutexN);
			while(s < nSamp && nT < maxT) {
				pthread_t thread;	
				int ret = 0;
				TypeFossilFeature *ffe;
				TypeTree *treeTmp;
				TypeThreadParameter *tp;
				int it;
fprintf(stdout, "\rSampling %d/%d", s, nSamp); fflush(stdout);
				it = gsl_rng_uniform_int(rg, (unsigned long int) sizeTree);
				treeTmp = cpyTree(tree[it]);
				treeTmp->minTime = minTimeIntInf;
				treeTmp->minTimeInt.inf = minTimeIntInf;
				treeTmp->minTimeInt.sup = minTimeIntSup;
				treeTmp->maxTime = maxTimeIntInf+gsl_rng_uniform(rg)*(maxTimeIntSup-maxTimeIntInf);
				ffe = sampleFossilInt(fos[it], treeTmp->size);
				for(n=0; n<treeTmp->size; n++) {
					if(treeTmp->node[n].child == NOSUCH) {
						int f;
						double maxTmp;
						switch(fos[it]->status[n]) {
							case contempNodeStatus:
								treeTmp->time[n] = treeTmp->maxTime;
							break;
							case unknownNodeStatus:
								
							break;
							case extinctNodeStatus:
								maxTmp = NEG_INFTY;
								for(f=ffe->fossil[n]; f!=NOSUCH; f=ffe->fossilList[f].prec)
									if(ffe->fossilList[f].time>maxTmp)
										maxTmp = ffe->fossilList[f].time; 
								treeTmp->time[n] = maxTmp;
							break;
							default:
								fprintf(stderr, "Node %d has no status\n", n);
								return 1;
						}
					} else
						treeTmp->time[n] = NO_TIME;
				}
				logD[s].size = def;
				logD[s].item = (TypeDistributionItem*) malloc(logD[s].size*sizeof(TypeDistributionItem));
				tp = (TypeThreadParameter*) malloc(sizeof(TypeThreadParameter));
				tp->tree = treeTmp;
				tp->ffe = ffe;
				tp->param = &param;
				tp->min = min;
				tp->max = max;
				tp->step = step;
				tp->def = def;
				tp->node = node[it];
				tp->logD = &(logD[s]);
				tp->logCond = &(logCond[s]);
				tp->mutex_number = &mutexN;
				tp->cond_number = &condN;
				tp->number = &nT;
				if((ret = pthread_create(&thread, NULL, threadComputeDistribution, (void*) tp)) == 0) {
					int err;
					if((err = pthread_detach(thread)) == 0) {
						nT++; s++;
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
		double sumCond, logSumCond, offset;
		offset = logCond[0];
		sumCond = 1.;
		for(s=1; s<nSamp; s++)
			if(!(isinf(logCond[s]) || isnan(logCond[s]))) {
				if(logCond[s]>offset) { /*compute max in offset just to avoid numerical precision issues*/
					sumCond *= exp(offset-logCond[s]);
					offset = logCond[s];
					sumCond++;
				} else
					sumCond += exp(logCond[s]-offset);
			}
		logSumCond = log(sumCond)+offset;
		dist.size = def;
		dist.item = (TypeDistributionItem*) malloc(dist.size*sizeof(TypeDistributionItem));
		for(i=0; i<dist.size; i++) {
			double sumI;
			dist.item[i].val = min + i*step;
			offset = logD[0].item[i].dens;
			sumI = 1.;
			for(s=1; s<nSamp; s++) {
				if(!(isinf(logD[s].item[i].dens) || isnan(logD[s].item[i].dens))) {
					if(logD[s].item[i].dens>offset) { /*compute max in offset just to avoid numerical precision issues*/
						sumI *= exp(offset-logD[s].item[i].dens);
						offset = logD[s].item[i].dens;
						sumI++;
					} else
						sumI += exp(logD[s].item[i].dens-offset);
				}
			}
			dist.item[i].dens = exp(log(sumI)+offset-logSumCond);
		}
		sprintf(outputDistribution, "%s.csv", outputPrefix);
		if((fo = fopen(outputDistribution, "w"))) {
			if(outDens) {
				deriveDistribution(&dist);
				fprintDistribution(fo, dist);
			} else
				fprintDistribution(fo, dist);
			fclose(fo);
		}
		if(dist.item != NULL)
			free((void*)dist.item);
		for(s=0; s<nSamp; s++)
			if(logD[s].item != NULL)
				free((void*)logD[s].item);
		gsl_rng_free(rg);
		free((void*)node);
	} else {
		fprintf(stderr, "Cannot read %s or %s or %s\n", inputFileNameTree, inputFileNameList, inputFileNameFossil);
		exit(1);
	}
	return 0;
}

void *threadComputeDistribution(void *data) {
	int i;
	for(i=0; i<((TypeThreadParameter*)data)->logD->size; i++)
		((TypeThreadParameter*)data)->logD->item[i].val = ((TypeThreadParameter*)data)->min + i*((TypeThreadParameter*)data)->step;
	fillLogDistribution(((TypeThreadParameter*)data)->logD, ((TypeThreadParameter*)data)->logCond, ((TypeThreadParameter*)data)->node, ((TypeThreadParameter*)data)->tree, ((TypeThreadParameter*)data)->ffe, ((TypeThreadParameter*)data)->param);
	freeTree(((TypeThreadParameter*)data)->tree);
	freeFossilFeature(((TypeThreadParameter*)data)->ffe);
	pthread_mutex_lock(((TypeThreadParameter*)data)->mutex_number);
		(*((TypeThreadParameter*)data)->number)--;
		pthread_cond_signal(((TypeThreadParameter*)data)->cond_number);
	pthread_mutex_unlock(((TypeThreadParameter*)data)->mutex_number);
	free(data);
	return NULL;
}

#define MAX_SIZE_TMP 50
#define INC_BUFFER 50
#define IS_SEP(c) (c == ' ' || c == '\t' || c == '\r' || c == '\n' || c == ';')

char **readList(FILE *f) {
	char c, tmp[MAX_SIZE_TMP+1], **list;
	int size, sizeBuffer;

	sizeBuffer = INC_BUFFER;
	list= (char**) malloc(sizeBuffer*sizeof(char*));
	size = 0;
	do {
		c = getc(f);
	} while(c!=EOF && IS_SEP(c)); 
	while(c != EOF) {
		int i;
		i = 0;
		while(i<MAX_SIZE_TMP && c !=EOF && !IS_SEP(c)) {
			tmp[i] = c;
			c = getc(f);
			i++;
		}
		tmp[i++] = '\0';
		if(i == MAX_SIZE_TMP) {
			fprintf(stderr, "Ident too long (%s) ...", tmp);
			exit(1);
		}
		if(i>1) {
			if(size >= sizeBuffer) {
				sizeBuffer += INC_BUFFER;
				list = (char**) realloc((void *) list, sizeBuffer*sizeof(char*));
			}
			list[size] = (char*) malloc((strlen(tmp)+1)*sizeof(char));
			strcpy(list[size], tmp);
			size++;
		}
		while(c!=EOF && IS_SEP(c))
			c=getc(f);
	}
	if(size >= sizeBuffer) {
		sizeBuffer += INC_BUFFER;
		list = (char**) realloc((void *) list, sizeBuffer*sizeof(char*));
	}
	list[size++] = NULL;
	return list;
}
