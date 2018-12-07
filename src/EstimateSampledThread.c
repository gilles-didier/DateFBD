#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <pthread.h>
#include <signal.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multimin.h>

#include "Utils.h"
#include "Tree.h"
#include "SimulTree.h"
#include "Fossil.h"
#include "FossilInt.h"
#include "Model.h"
#include "Uncertainty.h"
#include "SimulFossil.h"
#include "MinimizeNLOpt.h"


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

#define MAX_nSamp 1000

#define M_MAX 6
#define M_MAX_F 4
#define MIN_TREE 20
#define MAX_TRIALS 1000

#define HELPMESSAGE "--------------------------\n\nNAME\n	esti - estimates the speciation, extinction and fossilization rates \n	\nSYNOPSIS\n	esti [OPTIONS] <input Tree File> <input Fossil File> <input List Clade> [output File]\n\nDESCRIPTION\n	Estimate the speciation, extinction and fossilization rates by sampling into the trees contained in <input Tree File> (in Newick format), into the fossil ranges provided in <input Fossil File> (in csv format) and output the result as a text file <output File>\n\n	Options are\n	-o <origin bound inf> [origin bound sup]\n		set the origin range; if origin bound sup is not given it is set to the most ancient fossil age\n	-e <end bound inf> <end bound sup>\n		set the end time range\n	-s <number>\n		set the number of samples\n	-t <number>\n		set the number of thread running in parallell\n	-n <option File>\n		set the options of the numerical optimization to that provided in <option File>\n	-h\n		display help\n\n--------------------------\n\n"

typedef struct THREAD_PARAMETER {
	int *number;
	pthread_mutex_t *mutex_number;
	pthread_cond_t *cond_number;
	TypeTree *tree;
	TypeFossilFeature *ffe;
	int def, node;
	double min, max, step;
	TypeNLOptOption *nloptOption;
	TypeEstimation *estim;
} TypeThreadParameter;

static void *threadEstimateParameters(void *data);
static void fprintReport(FILE *f, TypeEstimation estimationMean, TypeEstimation estimationStd, TypeNLOptOption nloptOption);

void fprintReport(FILE *f, TypeEstimation estimationMean, TypeEstimation estimationStd, TypeNLOptOption nloptOption) {
	time_t t = time(NULL);
	struct tm tm = *localtime(&t);
	char buffer[50];
	fprintf(f, "Diversification execution report\n%d/%d/%d at %d:%d\n\n", tm.tm_mon + 1, tm.tm_mday, tm.tm_year + 1900, tm.tm_hour, tm.tm_min);
	fprintf(f, "\nEstimations\tMean\tStandard Deviation\n");
	sprintf(buffer, "%.6lE", estimationMean.param.birth);
	fprintf(f, "Speciation rate:\t%s", buffer);
	sprintf(buffer, "%.6lE", estimationStd.param.birth);
	fprintf(f, "\t%s\n", buffer);
	sprintf(buffer, "%.6lE", estimationMean.param.death);
	fprintf(f, "Extinction rate:\t%s", buffer);
	sprintf(buffer, "%.6lE", estimationStd.param.death);
	fprintf(f, "\t%s\n", buffer);
	sprintf(buffer, "%.6lE", estimationMean.param.fossil);
	fprintf(f, "Fossil find rate:\t%s", buffer);
	sprintf(buffer, "%.6lE", estimationStd.param.fossil);
	fprintf(f, "\t%s\n", buffer);
	sprintf(buffer, "%.6lE", estimationMean.logLikelihood);
	fprintf(f, "Log-likelihood:\t%s", buffer);
	sprintf(buffer, "%.6lE", estimationStd.logLikelihood);
	fprintf(f, "\t%s\n", buffer);
	fprintf(f, "\nOptimization settings\n");
	fprintNLoptOption(f, &nloptOption);
}

int main(int argc, char **argv) {	
	char *inputFileNameTree, *inputFileNameFossil, *outputFileName, option[256];
	FILE *fi, *fo, *ff;
	int i, j, nSamp = 10, maxT = 2;
	double minTimeIntInf = NO_TIME, minTimeIntSup = NO_TIME, maxTimeIntInf = 0., maxTimeIntSup = 0.;
	gsl_rng *rg = gsl_rng_alloc(gsl_rng_random_glibc2);
	TypeNLOptOption nloptOption;
	pthread_mutex_t mutexN = PTHREAD_MUTEX_INITIALIZER;
	pthread_cond_t condN = PTHREAD_COND_INITIALIZER;

	nloptOption.infSpe = 0.;
	nloptOption.supSpe = 1.;
	nloptOption.infExt = 0.;
	nloptOption.supExt = 1.;
	nloptOption.infFos = 0.;
	nloptOption.supFos = 1.;
	nloptOption.trials = 10;
	nloptOption.tolOptim = 0.00001;

	for(i=0; i<256; i++)
	option[i] = 0;

	for(i=1; i<argc && *(argv[i]) == '-'; i++) {
		for(j=1; argv[i][j] != '\0'; j++)
			option[argv[i][j]] = 1;
		if(option['n']) {
			option['n'] = 0;
			if((i+1)<argc) {
				FILE *fopt;
				if((fopt = fopen(argv[++i], "r"))) {
					fscanNLoptOptionTag(fopt, &nloptOption);
					fclose(fopt);
				} else {
					fprintf(stderr, "Can't open file %s\n", argv[++i]);
					exit(EXIT_FAILURE);
				}
			} else {
				fprintf(stderr, "File name missing after option -o\n");
				exit(EXIT_FAILURE);
			}
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
		if(option['s']) {
			option['s'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%d", &nSamp) == 1)
				i++;
			else
				exitProg(ErrorArgument, "A number expected after -p");
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
	if(i<argc)
		outputFileName = argv[i++];
	else
		outputFileName = "out.txt";
	if((fi = fopen(inputFileNameTree, "r")) && (ff = fopen(inputFileNameFossil, "r"))) {
		TypeTree **tree;
		TypeFossilIntFeature **fos;
		TypeEstimation mean, std, *estimation;
		int i, n, s, sizeTree, nT, cont;
		tree = readTrees(fi);
		fclose(fi);
		if(tree[0] == NULL) {
			fprintf(stderr, "Error: no tree\n");
			exit(1);
		}
		for(sizeTree=0; tree[sizeTree]!=NULL; sizeTree++) {
			if(tree[sizeTree]->name!=NULL)
				for(i=0; i<tree[sizeTree]->size; i++)
					if(tree[sizeTree]->name[i]!=NULL)
						fixSpace(tree[sizeTree]->name[i]);
		}	
		fos = (TypeFossilIntFeature**) malloc(sizeTree*sizeof(TypeFossilIntFeature*));
		for(i=0; i<sizeTree; i++) {
			toBinary(tree[i]);
			rewind(ff);
			fos[i] = getFossilIntFeature(ff, tree[i]->name, tree[i]->size);
			if(getMaxFossilIntTime(fos[i]) > 0.)
			negateFossilInt(fos[i]);
			fixStatus(tree[i], fos[i]);
		}
		estimation = (TypeEstimation*) malloc(nSamp*sizeof(TypeEstimation));

		double minFossilTime = getMinFossilIntTime(fos[0]);
		if(minTimeIntInf == NO_TIME || minTimeIntInf > minFossilTime) {
			if(minFossilTime<0)
				minTimeIntInf = 1.2*minFossilTime;
			else
				minTimeIntInf = 0.8*minFossilTime;
		}		s = 0;
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
				tp = (TypeThreadParameter*) malloc(sizeof(TypeThreadParameter));
				tp->tree = treeTmp;
				tp->ffe = ffe;
				tp->nloptOption = &nloptOption;
				tp->estim = &(estimation[s]);
				tp->mutex_number = &mutexN;
				tp->cond_number = &condN;
				tp->number = &nT;
				if((ret = pthread_create(&thread, NULL, threadEstimateParameters, (void*) tp)) == 0) {
					int err;
					if((err = pthread_detach(thread)) == 0) {
						nT++; s++;
					} else {
						fprintf (stderr, "Error %d while detaching thread: %s\n", err, (char*) strerror(err));
	//					pthread_kill(thread, 0);
					}
				} else
					fprintf (stderr, "Error %d while creating thread: %s\n", ret, (char*) strerror(ret));
			}
			cont = (nT > 0);
			if(cont)
				pthread_cond_wait(&condN, &mutexN);
			pthread_mutex_unlock(&mutexN);
		}
		printf("\n");
		int size = 0;
		mean.param.birth = 0.;
		mean.param.death = 0.;
		mean.param.fossil = 0.;
		mean.logLikelihood = 0.;
		for(s=0; s<nSamp; s++) {
			if(!(isnan(estimation[s].param.birth) || isnan(estimation[s].param.death) || isnan(estimation[s].param.fossil) || isnan(estimation[s].logLikelihood)) && estimation[s].logLikelihood > -1.00e+99) {
				mean.param.birth += estimation[s].param.birth;
				mean.param.death += estimation[s].param.death;
				mean.param.fossil += estimation[s].param.fossil;
				mean.logLikelihood += estimation[s].logLikelihood;
				size++;
			}
		}
		mean.param.birth /= (double) size;
		mean.param.death /= (double) size;
		mean.param.fossil /= (double) size;
		mean.logLikelihood /= (double) size;
		std.param.birth = 0.;
		std.param.death = 0.;
		std.param.fossil = 0.;
		std.logLikelihood = 0.;
		if(size>1) {
			for(s=0; s<nSamp; s++) {
				if(!(isnan(estimation[s].param.birth) || isnan(estimation[s].param.death) || isnan(estimation[s].param.fossil) || isnan(estimation[s].logLikelihood)) && estimation[s].logLikelihood > -1.00e+99) {
					std.param.birth += pow(estimation[s].param.birth-mean.param.birth, 2.);
					std.param.death += pow(estimation[s].param.death-mean.param.death, 2.);
					std.param.fossil += pow(estimation[s].param.fossil-mean.param.fossil, 2.);
					std.logLikelihood += pow(estimation[s].logLikelihood-mean.logLikelihood, 2.);
				}
			}
			std.param.birth /= (double) size-1.;
			std.param.death /= (double) size-1.;
			std.param.fossil /= (double) size-1.;
			std.logLikelihood /= (double) size-1.;
			std.param.birth = sqrt(std.param.birth);
			std.param.death = sqrt(std.param.death);
			std.param.fossil = sqrt(std.param.fossil);
			std.logLikelihood = sqrt(std.logLikelihood);
		}
		if((fo = fopen(outputFileName, "w"))) {
			fprintReport(fo, mean, std, nloptOption);
			fclose(fo);
		} else {
			fprintf(stderr, "Cannot open %s\n", outputFileName);
			exit(1);
		}
	} else {
		fprintf(stderr, "Cannot read %s or %s\n", inputFileNameTree, inputFileNameFossil);
		exit(1);
	}
	return 0;
}

void *threadEstimateParameters(void *data) {
	if(!minimizeBirthDeathFossilFromTreeFossil(getLogLikelihoodTreeFossil, ((TypeThreadParameter*)data)->tree, ((TypeThreadParameter*)data)->ffe, ((TypeThreadParameter*)data)->nloptOption, ((TypeThreadParameter*)data)->estim)) {
		fprintf(stderr, "Estimation issue\n");
		((TypeThreadParameter*)data)->estim->param.birth = sqrt(-1.);
	}
	freeTree(((TypeThreadParameter*)data)->tree);
	freeFossilFeature(((TypeThreadParameter*)data)->ffe);
	pthread_mutex_lock(((TypeThreadParameter*)data)->mutex_number);
		(*((TypeThreadParameter*)data)->number)--;
		pthread_cond_signal(((TypeThreadParameter*)data)->cond_number);
	pthread_mutex_unlock(((TypeThreadParameter*)data)->mutex_number);
	free(data);
	return NULL;
}
