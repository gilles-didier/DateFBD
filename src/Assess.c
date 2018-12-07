#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <pthread.h>
#include <omp.h>
#include <gsl/gsl_rng.h>

#include "Utils.h"
#include "Tree.h"
#include "SimulTree.h"
#include "Fossil.h"
#include "Model.h"
#include "Uncertainty.h"
#include "SimulFossil.h"
#include "MinimizeNLOpt.h"
#include "GNUFile.h"

#ifdef DO_PS
#endif

#define STRING_SIZE 300
#define INC_CUT 5
#define SEQ_SIZE 30
#define EXT_OUTPUT "_added.phy"


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
#define MAX_NF 10

#define M_MAX 6
#define M_MAX_F 4
#define MIN_TREE 20
#define MAX_PRED 7.5E+07


#define N_METHODS 2

//./asso -b 1.5 -d 1. -f 1. -f 2. -f 3. -x 8 -l 5 7 0.5 -i 100 -m 10 -M 1000 -c 5.E+07 BB
//./asso -b 1.5 -d 1. -f 1. -f 2. -f 3. -x 40 -l 5 7 0.5 -i 500 -m 10 -M 1000 -c 5.E+07 AA

#define HELP_MESSAGE "\nusage: assess [options] [<output file>]\n\nassess simulates random trees and fossils finds, estimates speciation and extinction rates and returns the mean absolute error\n\n\nOptions:\n\t-h : display help\n\t-b <birth>\t: set birth rate\n\t-d <death>\t: set death rate\n\t-f <fossil>\t: set fossil find rate\n\t-m <min>\t: set minimum number of contemporary species of a simulation to be considered\n\t-M <size>\t: set maximum size of a simulation to be considered\n\t-i <niter>\t: set the number of simulations\n\t-l <start> <end> <step> : the range of the end time of the diversification (start is always 0)\n\t-o <options file name>\t: load the settings of the optimizer. <options file name> has to be in the format:\n\t\t:SPE [0;1] :EXT [0;1] :FOS [0:1] :TRI 10 :TOL 1.E-7 :ITE 500\n\t-x <number>\t: set the number of threads\n\t-c <complexity>\t: set the max complexity index of a simulation to be considered\n\n\n"

typedef double TypeErrorFunction(double, double);

typedef struct OBSERVATION {
	double length;
	int number;
} TypeObservation;


typedef struct NODE_FOSSIL {
	int node, fossil;
	double time;
} TypeNodeFossil;

typedef struct NODE_FOSSIL_TABLE {
	int size;
	TypeNodeFossil *table;
} TypeNodeFossilTable;


static double AbsoluteError(double x, double y);

typedef struct {
	double birth, death, *fossil, start, step, **resBirth, **resDeath, **resFossil;
	int **resStat;
	int nf, niter, size, *iter, srbd, srf;
	int minContemp, maxSizeTree;
	double maxComplexity;
	TypeNLOptOption nloptOption;
} TypeThreadArgument;


typedef struct {
	int number, maxThreads, *cur;
	gsl_rng *r;
	TypeThreadArgument arg;
	pthread_mutex_t mutex_cur;
	pthread_mutex_t mutex_rand;
	pthread_mutex_t mutex_number;
	pthread_mutex_t mutex_arg;
	pthread_cond_t cond_number;
} TypeThreadData;

 
TypeThreadData set =
{
	.number = 0,
	.maxThreads = 1,
	.mutex_cur = PTHREAD_MUTEX_INITIALIZER,
	.mutex_rand = PTHREAD_MUTEX_INITIALIZER,
	.mutex_number = PTHREAD_MUTEX_INITIALIZER,
	.mutex_arg = PTHREAD_MUTEX_INITIALIZER,
	.cond_number = PTHREAD_COND_INITIALIZER
};

static void *simulThread(void *data) {
	TypeTree *tree = (TypeTree*) data;
	double *bufferB = (double*) malloc(set.arg.srbd*sizeof(double));
	double *bufferD = (double*) malloc(set.arg.srbd*sizeof(double));
	double *bufferF = (double*) malloc(set.arg.srf*sizeof(double));
	int *bufferS = (int*) malloc(3*sizeof(int)), i = (int) round((tree->maxTime-set.arg.start)/set.arg.step);;
	double totTime;
	double pred = 0.;
	TypeTree *tree1, *tree2;
	TypeListEvent *event;
	TypeEstimation estimation;
	int ok = 1;
	int m = 0, u, mc;
	event = getEventSequenceBD(tree);
	int nBirth, nDeath, nFossil;
	getStatEvent(event, &nBirth,  &nDeath,  &nFossil);
	totTime = getTotalTimeEvent(event);
	freeListEvent(event);
	bufferS[0] = nBirth;
	bufferS[1] = nDeath;
	bufferS[2] = nBirth-nDeath;
	bufferB[m] = ((double)nBirth)/totTime;
	bufferD[m] = ((double)nDeath)/totTime;
	m++;
	mc = m; /*keep the comtemporary for later since it is almost always OK*/
	m++;
	for(u=0; ok && u<set.arg.nf; u++) {
		TypeFossilFeature *fos;
		int n, i;	
		pthread_mutex_lock (&set.mutex_rand);
		fos = addFossils(set.r, set.arg.fossil[u], tree);
		pthread_mutex_unlock (&set.mutex_rand);
		event = getEventSequenceFBD(tree, fos);
		getStatEvent(event, &nBirth,  &nDeath,  &nFossil);
		freeListEvent(event);
		bufferF[u] = ((double)nFossil)/totTime;
		tree1 = pruneFossil(tree, fos);
		freeFossilFeature(fos);
		tree2 = fixBinaryFossil(tree1,  (TypeFossilFeature*) tree1->info);
		freeFossilFeature((TypeFossilFeature*) tree1->info);
		freeTree(tree1);
		if(1) {
			double fossilRate;
			i = 0;
			event = getEventSequenceFBD(tree2, (TypeFossilFeature*) tree2->info);
			fossilRate = 0.;
			for(n=0; n<event->size; n++)
				if(event->list[n].type == 'f')
					fossilRate++;
			fossilRate /= getTotalTimeEvent(event);
			estimation.param.fossil = fossilRate;
			estimation.param.fossil = 0.;
			for(n=0; n<event->size; n++)
				if(event->list[n].type == 'f')
					estimation.param.fossil++;
			estimation.param.fossil /= getTotalTimeEvent(event);
			if(ok && minimizeBirthDeathFromEventList(logProbEventFBD, event, &(set.arg.nloptOption), &estimation)>=0) {
				bufferB[i*set.arg.nf+m+u] = estimation.param.birth;
				bufferD[i*set.arg.nf+m+u] = estimation.param.death;
				bufferF[(i+1)*set.arg.nf+u] = estimation.param.fossil;
			} else
				ok = 0;
			freeListEvent(event);
			i++;
			double *save;
			save = tree2->time;
			tree2->time = (double*) malloc(tree2->size*sizeof(double));
			for(n=0; n<tree2->size; n++)
				if(tree2->node[n].child == NOSUCH)
					tree2->time[n] = save[n];
				else
					tree2->time[n] = NO_TIME;
			if(ok && minimizeBirthDeathFossilFromTreeFossil(getLogLikelihoodTreeFossil, tree2, (TypeFossilFeature*) tree2->info, &(set.arg.nloptOption), &estimation)) {
				bufferB[i*set.arg.nf+m+u] = estimation.param.birth;
				bufferD[i*set.arg.nf+m+u] = estimation.param.death;
				bufferF[(i+1)*set.arg.nf+u] = estimation.param.fossil;
			} else
				ok = 0;
			free((void*)tree2->time);
			tree2->time = save;
		} else
			ok = 0;
		freeFossilFeature((TypeFossilFeature*) tree2->info);
		freeTree(tree2);
	}
	if(ok) {
		tree1 = pruneContemp(tree);
		tree2 = fixBinary(tree1);
		freeTree(tree1);
		event = getEventSequenceBD(tree2);
		freeTree(tree2);
		if(minimizeBirthDeathFromEventList(logProbEventBD, event, &(set.arg.nloptOption), &estimation)>=0) {
			bufferB[mc] = estimation.param.birth;
			bufferD[mc] = estimation.param.death;
		} else
			ok = 0;
		freeListEvent(event);
	}
	if(ok) {
		pthread_mutex_lock (&set.mutex_arg);
		if(set.arg.iter[i]<set.arg.niter) {
			int m;
			for(m=0; m<3; m++)
				set.arg.resStat[m][i] = bufferS[m];
			for(m=0; m<set.arg.srbd; m++) {
				set.arg.resBirth[m][i] += AbsoluteError(set.arg.birth, bufferB[m]);
				set.arg.resDeath[m][i] += AbsoluteError(set.arg.death, bufferD[m]);
			}
			for(m=0; m<=N_METHODS; m++)
				for(u=0; u<set.arg.nf; u++)
					set.arg.resFossil[m*set.arg.nf+u][i] += AbsoluteError(set.arg.fossil[u], bufferF[m*set.arg.nf+u]);
			set.arg.iter[i]++;
		}
		pthread_mutex_unlock (&set.mutex_arg);
fprintf(stderr, "ok thread %d\n", i);
	} else
fprintf(stderr, "ko thread %d (%.1lE)\n", i, pred);
	freeTree(tree);
	free((void*)bufferB);
	free((void*)bufferD);
	free((void*)bufferF);
	free((void*)bufferS);
	pthread_mutex_lock (&set.mutex_cur);
	set.cur[i]--;
fprintf(stderr, "fin thread %d -> %d\n", i, set.cur[i]);
	pthread_mutex_unlock (&set.mutex_cur);
	pthread_mutex_lock (&set.mutex_number);
	set.number--;
	pthread_cond_signal (&set.cond_number);
	pthread_mutex_unlock (&set.mutex_number);
	return NULL;
}


static void *mainThread (void *data) {
	int i, cont = 1;
	set.cur = (int*) malloc(set.arg.size*sizeof(int));
	for(i=0; i<set.arg.size; i++)
		set.cur[i] = 0;
	cont = 1;
	i = 0;
	while(cont) {
		pthread_mutex_lock(&set.mutex_number);
		pthread_mutex_lock(&set.mutex_cur);
		pthread_mutex_lock(&set.mutex_arg);
		while(i < set.arg.size && set.number < set.maxThreads) {
			pthread_t thread;
			TypeTree *tree = NULL;
			double maxTime = set.arg.start+((double)i)*set.arg.step;
			do {
				if(tree != NULL)
					freeTree(tree);
				pthread_mutex_lock (&set.mutex_rand);
				tree = simulTree(set.r, set.arg.birth, set.arg.death, maxTime);
				pthread_mutex_unlock (&set.mutex_rand);
			} while(!(tree != NULL && (countContemp(tree))>=set.arg.minContemp && tree->size<=set.arg.maxSizeTree));
fprintf (stderr, "Launching thread %d on simulated tree with %.1lf, %.1lf during %.1lf processing %d, completed %d (%d)\n", set.number, set.arg.birth, set.arg.death, maxTime, set.cur[i], set.arg.iter[i], set.arg.niter);
			int ret = 0;
			if((ret = pthread_create(&thread, NULL, simulThread, (void*) tree)) == 0) {
				int err;
				if((err = pthread_detach(thread)) == 0) {
					set.number++;
					set.cur[i]++;
				} else {
					fprintf (stderr, "Error %d while detaching thread: %s\n", err, (char*) strerror (err));
//					pthread_kill(thread, 0);
				}
			} else
				fprintf (stderr, "Error %d while creating thread: %s\n", ret, (char*) strerror (ret));
fprintf(stderr, "iter %d cur %d, comp %d", i, set.cur[i], set.arg.iter[i]);
			for(i=0; i<set.arg.size && (set.cur[i]+set.arg.iter[i])>=set.arg.niter; i++)
				;
fprintf(stderr, " -> %d\n", i);
		}
		cont = (set.number > 0) || (i < set.arg.size);
		pthread_mutex_unlock (& set.mutex_number);
		pthread_mutex_unlock(&set.mutex_arg);
		pthread_mutex_unlock(&set.mutex_cur);
		if(cont) {
			pthread_mutex_lock(&set.mutex_number);
			pthread_cond_wait (& set.cond_number, & set.mutex_number);
			pthread_mutex_unlock (& set.mutex_number);
			pthread_mutex_lock(&set.mutex_cur);
			pthread_mutex_lock(&set.mutex_arg);
			for(i=0; i<set.arg.size && (set.cur[i]+set.arg.iter[i])>=set.arg.niter; i++)
				;
			pthread_mutex_unlock(&set.mutex_arg);
			pthread_mutex_unlock(&set.mutex_cur);
			pthread_mutex_lock(&set.mutex_number);
			cont = (set.number > 0) || (i < set.arg.size);
			pthread_mutex_unlock (& set.mutex_number);
		}
	}
	pthread_mutex_lock(&set.mutex_cur);
	pthread_mutex_lock(&set.mutex_arg);
		for(i=0; i<set.arg.size && set.arg.iter[i]>=set.arg.niter; i++)
			printf("%d\titer %d, cur %d\n", i, set.arg.iter[i], set.cur[i]);
	pthread_mutex_unlock(&set.mutex_cur);
	pthread_mutex_unlock(&set.mutex_arg);
	free((void*)set.cur);
fprintf (stderr, "End %d\n", set.number);
	return NULL;
}

int main(int argc, char **argv) {
	int i;
	char  outputFileName[SIZE_BUFFER_CHAR], bufferOutput[SIZE_BUFFER_CHAR], *outtrunc, option[256];
	double birth = 2., death = 1., start=2, end=5, step=0.5;
	int niter = 20, m;
	double **smBirth, **smDeath, **smFossil, fossilTab[MAX_NF];
	int **smStat, nf = 0;
	time_t t1, t2;
	gsl_rng *rg = gsl_rng_alloc(gsl_rng_random_glibc2);
	TypeNLOptOption nloptOption;
	int minContemp = 10, maxSizeTree = 10000;
	double maxComplexity = 7.5E+07;
		
    nloptOption.infSpe = 0.;
    nloptOption.supSpe = 1.;
    nloptOption.infExt = 0.;
    nloptOption.supExt = 1.;
    nloptOption.infFos = 0.;
    nloptOption.supFos = 1.;
    nloptOption.trials = 10;
    nloptOption.tolOptim = 0.00001;
    nloptOption.maxIter = 500;
	
	for(i=0; i<256; i++)
		option[i] = 0;
	for(i=1; i<argc && *(argv[i]) == '-'; i++) {
		int j;
		for(j=1; argv[i][j] != '\0'; j++)
			option[argv[i][j]] = 1;
		if(option['b']) {
			option['b'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &birth) == 1)
				i++;
		}
		if(option['d']) {
			option['d'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &death) == 1)
				i++;
		}
		if(option['f']) {
			option['f'] = 0;
			if(nf < MAX_NF) {
				if((i+1)<argc && sscanf(argv[i+1], "%lf", &(fossilTab[nf])) == 1) {
					i++;
					nf++;
				} else
					printf("Can't read \"%s\" as fossil rate\n", argv[i+1]);
			} else
				printf("Can't handle more than %d fossil rates\n", MAX_NF);
		}
		if(option['i']) {
			option['i'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%d", &niter) == 1)
				i++;
		}
		if(option['m']) {
			option['m'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%d", &minContemp) == 1)
				i++;
		}
		if(option['M']) {
			option['M'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%d", &maxSizeTree) == 1)
				i++;
		}
		if(option['c']) {
			option['c'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%lE", &maxComplexity) == 1)
				i++;
		}
		if(option['l']) {
			option['l'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &start) == 1)
				i++;
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &end) == 1)
				i++;
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &step) == 1)
				i++;
		}
		if(option['x']) {
			option['x'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%d", &(set.maxThreads)) == 1)
				i++;
		}
		if(option['o']) {
			FILE *fopt;
			if((i+1)<argc) {
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
		if(option['h']) {
			printf("%s\n", HELP_MESSAGE);
			exit(EXIT_SUCCESS);
		}
	}

	if (i>=argc || sscanf(argv[i++], "%s", outputFileName) != 1)
		strcpy(outputFileName, "out");
	if((outtrunc = strrchr(outputFileName, '/')) != NULL)
		outtrunc++;
	else
		outtrunc = outputFileName;
	fprintf(stdout, "Run with optimizer options:\n");
	fprintNLoptOption(stdout, &nloptOption);
	fprintf(stdout, "\n\n");
t1 = time(&t1);
	{
		FILE *fmb, *fmd, *fmf, *fgnu;
		int size, srbd, srf, i;
		char BUFFER[N_METHODS*MAX_NF][SIZE_BUFFER_CHAR], bufferO[SIZE_BUFFER_CHAR], bufferD[SIZE_BUFFER_CHAR], **TITLE;
		TypeGNUInfo gnu;
		srbd = 2+N_METHODS*nf;
		srf = (N_METHODS+1)*nf;

		TITLE = (char**) malloc(utils_MAX(srbd, srf)*sizeof(char*));
/*Gnuplot general info*/
		gnu.type = EPS;
		gnu.xLabel = "time";
		gnu.yLabel = "absolute error";
		gnu.title = TITLE;
		gnu.outputFileName = bufferO;
		gnu.dataFileName = bufferD;
		TITLE[0] = "complete information";
		TITLE[1] = "contemporary lineages + divergence times";
		i = 0;
		for(m=0; m<nf; m++) {
			sprintf(BUFFER[2+i*nf+m], "contemporary lineages + divergence times + fossil finds with rate %.2lf", fossilTab[m]);
			TITLE[2+i*nf+m] = BUFFER[2+i*nf+m];
		}
		i++;
		for(m=0; m<nf; m++) {
			sprintf(BUFFER[2+i*nf+m], "contemporary lineages + fossil finds with rate %.2lf", fossilTab[m]);
			TITLE[2+i*nf+m] = BUFFER[2+i*nf+m];
		}
		gnu.nColumn = 2+N_METHODS*nf;
/*Birth files*/
		sprintf(bufferO, "%s_Abs_Mean_birth_r%.1lf_%d.eps", outtrunc, birth/death, niter);
		sprintf(bufferD, "%s_Abs_Mean_birth_r%.1lf_%d.csv", outputFileName, birth/death, niter);
		sprintf(bufferOutput, "%s_Abs_Mean_birth_r%.1lf_%d.gnu", outputFileName, birth/death, niter);
		if((fgnu = fopen(bufferOutput, "w"))) {
			fprintGNUFile(fgnu, &gnu);
			fclose(fgnu);
		} else {
			fprintf(stderr, "Can't open file %s\n", bufferOutput);
			exit(EXIT_FAILURE);
		}
		if(!(fmb = fopen(bufferD, "w"))) {
			fprintf(stderr, "Can't open file %s\n", bufferD);
			exit(EXIT_FAILURE);
		}
/*Death files*/
		sprintf(bufferO, "%s_Abs_Mean_death_r%.1lf_%d.eps", outtrunc, birth/death, niter);
		sprintf(bufferD, "%s_Abs_Mean_death_r%.1lf_%d.csv", outputFileName, birth/death, niter);
		sprintf(bufferOutput, "%s_Abs_Mean_death_r%.1lf_%d.gnu", outputFileName, birth/death, niter);
		if((fgnu = fopen(bufferOutput, "w"))) {
			fprintGNUFile(fgnu, &gnu);
			fclose(fgnu);
		} else {
			fprintf(stderr, "Can't open file %s\n", bufferOutput);
			exit(EXIT_FAILURE);
		}
		if(!(fmd = fopen(bufferD, "w"))) {
			fprintf(stderr, "Can't open file %s\n", bufferD);
			exit(EXIT_FAILURE);
		}
/*Fossil files*/
		gnu.nColumn = (N_METHODS+1)*nf;
		i = 0;
		for(m=0; m<nf; m++) {
			sprintf(BUFFER[i*nf+m], "complete information + fossil finds with rate %.2lf", fossilTab[m]);
			TITLE[i*nf+m] = BUFFER[i*nf+m];
		}
		i++;
		for(m=0; m<nf; m++) {
			sprintf(BUFFER[i*nf+m], "contemporary lineages + divergence times + fossil finds with rate %.2lf", fossilTab[m]);
			TITLE[i*nf+m] = BUFFER[i*nf+m];
		}
		i++;
		for(m=0; m<nf; m++) {
			sprintf(BUFFER[i*nf+m], "contemporary lineages + fossil finds with rate %.2lf", fossilTab[m]);
			TITLE[i*nf+m] = BUFFER[i*nf+m];
		}
		sprintf(bufferO, "%s_Abs_Mean_fossil_r%.1lf_%d.eps", outtrunc, birth/death, niter);
		sprintf(bufferD, "%s_Abs_Mean_fossil_r%.1lf_%d.csv", outputFileName, birth/death, niter);
		sprintf(bufferOutput, "%s_Abs_Mean_fossil_r%.1lf_%d.gnu", outputFileName, birth/death, niter);
		if((fgnu = fopen(bufferOutput, "w"))) {
			fprintGNUFile(fgnu, &gnu);
			fclose(fgnu);
		} else {
			fprintf(stderr, "Can't open file %s\n", bufferOutput);
			exit(EXIT_FAILURE);
		}
		if(!(fmf = fopen(bufferD, "w"))) {
			fprintf(stderr, "Can't open file %s\n", bufferD);
			exit(EXIT_FAILURE);
		}
		free((void*)TITLE);

		size = (int) floor((end-start+step/10.)/step)+1;
		smBirth = (double**) malloc(srbd*sizeof(double));
		smDeath = (double**) malloc(srbd*sizeof(double));
		smFossil = (double**) malloc(srf*sizeof(double));
		for(m=0; m<srbd; m++) {
			smBirth[m] = (double*) malloc(size*sizeof(double));
			smDeath[m] = (double*) malloc(size*sizeof(double));
			for(i=0; i<size; i++) {
				smBirth[m][i] = 0.;
				smDeath[m][i] = 0.;
			}
		}
		for(m=0; m<srf; m++) {
			smFossil[m] = (double*) malloc(size*sizeof(double));
			for(i=0; i<size; i++)
				smFossil[m][i] = 0.;
		}
		smStat = (int**) malloc(3*sizeof(int*));
		for(m=0; m<3; m++) {
			smStat[m] = (int*) malloc(size*sizeof(int));
			for(i=0; i<size; i++)
				smStat[m][i] = 0.;
		}
		set.r = rg;
		set.arg.minContemp = minContemp;
		set.arg.maxSizeTree = maxSizeTree;
		set.arg.maxComplexity = maxComplexity;
		set.arg.nloptOption = nloptOption;	
		set.arg.birth = birth;
		set.arg.death = death;
		set.arg.fossil = fossilTab;
		set.arg.start = start;
		set.arg.step = step;
		set.arg.resBirth = smBirth;
		set.arg.resDeath = smDeath;
		set.arg.resFossil = smFossil;
		set.arg.resStat = smStat;
		set.arg.nf = nf;
		set.arg.niter = niter;
		set.arg.size = size;
		set.arg.iter = (int*) malloc(size*sizeof(int));
		for(i=0; i<size; i++)
			set.arg.iter[i] = 0;
		set.arg.srbd = srbd;
		set.arg.srf = srf;
		int ret;
		pthread_t thread;
		if((ret = pthread_create (&thread, NULL, mainThread, NULL)) == 0) {
			pthread_join (thread, NULL);
		} else {
			fprintf (stderr, "Erreur %d %s\n", ret, (char*) strerror (ret));
			exit(1);
		}
		for(i=0; i<size; i++) {
			double maxTime;
			int m;
			maxTime = start+((double)i)*step;
			fprintf(fmb, "%lf", maxTime);
			for(m=0; m<srbd; m++)
				fprintf(fmb, "\t%lE", set.arg.resBirth[m][i]/((double)set.arg.iter[i]));
			fprintf(fmb, "\n");
			fprintf(fmd, "%lf", maxTime);
			for(m=0; m<srbd; m++)
				fprintf(fmd, "\t%lE", set.arg.resDeath[m][i]/((double)set.arg.iter[i]));
			fprintf(fmd, "\n");
			fprintf(fmf, "%lf", maxTime);
			for(m=0; m<srf; m++)
				fprintf(fmf, "\t%lE", set.arg.resFossil[m][i]/((double)set.arg.iter[i]));
			fprintf(fmf, "\n");
		}
		fclose(fmb);
		fclose(fmd);
		fclose(fmf);
		for(m=0; m<srbd; m++) {
			free((void*)smBirth[m]);
			free((void*)smDeath[m]);
		}
		for(m=0; m<srf; m++)
			free((void*)smFossil[m]);
		for(m=0; m<3; m++)
			free((void*)smStat[m]);
		free((void*)smBirth);
		free((void*)smDeath);
		free((void*)smFossil);
		free((void*)smStat);
	}
t2 = time(&t2);
printf("\nTotal time %.1lf\n", difftime(t2, t1));
	gsl_rng_free(rg);
	exit(EXIT_SUCCESS);
}


double AbsoluteError(double x, double y) {
	return fabs(x-y);
}


