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
#include "DrawTreeCairo.h"
#include "DrawTreePSTricks.h"
#include "DrawTreeTikz.h"
#include "DrawTreeGeneric.h"
#include "DrawFossilInt.h"
#include "DrawDensity.h"

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
#define MAX_NODE 100

#define HELPMESSAGE "--------------------------\n\nNAME\n	tree - Computation of the divergence time distibutions associated to a given set of nodes from a single tree, the fossil ages, and the diversification rates and draw a tree figure\nSYNOPSIS\n	tree [OPTIONS] <input Tree File> <input Fossil File> [output File]\n\nDESCRIPTION\n	Compute the divergence time distibutions associated to a given set of nodes from the tree of <input Tree File>, the fossil ages of <input Fossil File>, and the diversification rates and draw a tree figure\n\n	Options are\n	-z <input Tree File>\n		output the tree in text debug format in the console and exit \n	-o <origin bound inf> [origin bound sup]\n		set the origin range; if origin bound sup is not given it is set to the most ancient fossil age\n	-e <end bound inf> <end bound sup>\n		set the end time range\n	-p <speciation rate> <extinction rate> <fossilization rate>\n		set the speciation, extinction and fossilization rates\n	-s <number>\n		set the number of samples\n	-d\n		return the distribution (otherwise the density is returned by default)\n	-u <value>\n		set the step discretizing the time distribution to <value>\n	-t <number>\n		set the number of thread running in parallell\n	-n <number>\n		compute the divergence distribution associated to node <number>; option can be used several times; if it is not used all the divergence times are computed\n	-x <number>\n		set the graphic format of the output (option is required if one wants a graphic output)\n			-f 1 -> pdf\n			-f 2 -> postscript\n			-f 3 -> png\n			-f 4 -> svg\n			-f 5 -> LaTeX (psTricks)\n			-f 6 -> LaTeX (TikZ)\n	-h\n		display help\n\n--------------------------\n\n"

typedef struct DATA_DRAW_FOSSIL_DENSITY {
	TypeDataDrawFossilInt *dataFossil;
	TypeDataDrawDensity *dataDensity;
} TypeDataDrawFossilDensity;


void drawFossilDensity(int n , double x, double y, TypeInfoDrawTreeGeneric *info, void *data) {
	drawFossilInt(n, x, y, info, (void*) ((TypeDataDrawFossilDensity*)data)->dataFossil);
	drawDensity(n, x, y, info, (void*) ((TypeDataDrawFossilDensity*)data)->dataDensity);
}

typedef struct THREAD_PARAMETER {
	int *number;
	pthread_mutex_t *mutex_number;
	pthread_cond_t *cond_number;
	TypeTree *tree;
	TypeFossilFeature *ffe;
	int *todo, nTodo;
	double *min, *max, step;
	TypeDistribution *logD;
	double *logCond;
	TypeModelParam *param;
} TypeThreadParameter;


static void *threadComputeDistribution(void *data);


int main(int argc, char **argv) {	
	char *inputFileNameTree, *inputFileNameFossil = NULL, *outputPrefix = PREFIX, outputDistribution[STRING_SIZE], option[256], format = '1';
	FILE *fi, *fo;
	int i, j, nSamp = 100, node = NOSUCH, outDens = 1, maxT = 2, *todo, tab_todo[MAX_NODE], nTodo=0;
	double minTimeIntInf = NO_TIME, minTimeIntSup = NO_TIME, maxTimeIntInf = 0., maxTimeIntSup = 0., step = 0.1, stepFig = 10.;
	TypeModelParam param = {.birth=1.3, .death = 1., .fossil = 1., .sampling = 1.};
	pthread_mutex_t mutexN = PTHREAD_MUTEX_INITIALIZER;
	pthread_cond_t condN = PTHREAD_COND_INITIALIZER;

	
	for(i=0; i<256; i++)
		option[i] = 0;
	   
	for(i=1; i<argc && *(argv[i]) == '-'; i++) {
		for(j=1; argv[i][j] != '\0'; j++)
			option[argv[i][j]] = 1;
		if(option['x']) {
			option['x'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%c", &format) == 1)
				i++;
			else
				exitProg(ErrorArgument, "a character is required after option -f");
		}
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
				exitProg(ErrorArgument, "at least one value is expected after -e");
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
		if(option['d']) {
			option['d'] = 0;
			outDens = 0;
		}
		if(option['s']) {
			option['s'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%d", &nSamp) == 1)
				i++;
			else
				exitProg(ErrorArgument, "a number is expected after -s");
		}
		if(option['n']) {
			option['n'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%d", &(node)) == 1) {
				if(nTodo < MAX_NODE)
					tab_todo[nTodo++] = node;
				else {
					fprintf(stderr, "Warning: you get specify more than %d nodes to check.\n", MAX_NODE);
				}
				i++;
			} else
				exitProg(ErrorArgument, "a node index are expected after -n");
		}
		if(option['d']) {
			option['d'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &(step)) == 1)
				i++;
			else
				exitProg(ErrorArgument, "a number is expected after -d");
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
		if(option['q']) {
			option['q'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &stepFig) == 1)
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
		outputPrefix = argv[i++];
	if((fi = fopen(inputFileNameTree, "r"))) {
		TypeTree *tree;
		TypeFossilIntFeature *fos;
		TypeDistribution **logD, *d;
		double **logCond, *min, *max;
		int i, n, s, nT, cont;
		TypeInfoDrawTreeGeneric info;
		TypeAdditionalDrawTreeGeneric add;
		TypeDataDrawFossilInt data;
		TypeDataDrawDensity dataD;
		TypeDataDrawFossilDensity dataFD;
		gsl_rng *rg = gsl_rng_alloc(gsl_rng_random_glibc2);
		
        tree = readTree(fi);
        fclose(fi);
		toBinary(tree);
		reorderTreeSize(tree);
		if(inputFileNameFossil != NULL) {
			FILE *ff;
			if((ff = fopen(inputFileNameFossil, "r"))) {
				fos = getFossilIntFeature(ff, tree->name, tree->size);
			} else {
				fprintf(stderr, "Cannot read %s\n", inputFileNameFossil);
				exit(1);
			}
		} else
			fos = fosToFossilInt(tree);
		if(getMaxFossilIntTime(fos) > 0.)
			negateFossilInt(fos);
		fixStatus(tree, fos);
		if(nTodo>0)
			todo = tab_todo;
		else {
			int n;
			todo = (int*) malloc(sizeof(int)*tree->size/2);
			for(n=0; n<tree->size; n++)
				if(tree->node[n].child != NOSUCH)
					todo[nTodo++] = n;
		}
		for(n=0; n<tree->size; n++)
			if(tree->node[n].child != NOSUCH)
				tree->time[n] = NO_TIME;
		min = (double*) malloc(nTodo*sizeof(double));
		max = (double*) malloc(nTodo*sizeof(double));
		if(tree->minTimeInt.inf != NO_TIME)
			minTimeIntInf = tree->minTimeInt.inf;
		if(tree->minTimeInt.sup != NO_TIME)
			minTimeIntSup = tree->minTimeInt.sup;
		if(tree->maxTimeInt.inf != NO_TIME)
			maxTimeIntInf = tree->maxTimeInt.inf;
		if(tree->maxTimeInt.sup != NO_TIME)
			maxTimeIntSup = tree->maxTimeInt.sup;
		double minFossilTime = getMinFossilIntTime(fos);
		if(minTimeIntInf == NO_TIME && minTimeIntSup == NO_TIME) {
			double tmp = getMinFossilIntTime(fos);
			if(tmp<0) {
				minTimeIntInf = 1.2*minFossilTime;
				minTimeIntSup = 1.1*minFossilTime;
			} else {
				minTimeIntInf = 0.8*minFossilTime;
				minTimeIntSup = 0.9*minFossilTime;
			}
		}
		if(minTimeIntInf > minFossilTime)
			minTimeIntInf = minFossilTime;
		tree->maxTime = maxTimeIntSup;
		tree->minTime = minTimeIntInf;
		tree->minTimeInt.inf = minTimeIntInf;
		tree->minTimeInt.sup = minTimeIntSup;
		tree->maxTime = maxTimeIntSup;
		tree->maxTimeInt.inf = maxTimeIntInf;
		tree->maxTimeInt.sup = maxTimeIntSup;
		if(tree->parent==NULL)
			tree->parent = getParent(tree);
		for(n=0; n<nTodo; n++) {
			min[n] = getMaxFossilIntTimeToNode(todo[n], tree, fos);
			max[n] = utils_MIN(getMinFossilIntTimeFromNode(todo[n], tree, fos), getMinTimeFromNode(todo[n], tree));
			min[n] = floor(min[n]/step)*step;
			max[n] = (ceil(max[n]/step))*step;
		}			
		logD = (TypeDistribution**) malloc(nSamp*sizeof(TypeDistribution*));
		logCond = (double**) malloc(nSamp*sizeof(double*));
		for(s=0; s<nSamp; s++) {
			logD[s] = (TypeDistribution*) malloc(nTodo*sizeof(TypeDistribution));
			logCond[s] = (double*) malloc(nTodo*sizeof(double));
				for(n=0; n<nTodo; n++) {
					logD[s][n].size = ceil((max[n]-min[n])/step)+1;
					logD[s][n].item = (TypeDistributionItem*) malloc(logD[s][n].size*sizeof(TypeDistributionItem));
				}
		}
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
fprintf(stdout, "\rSampling %d/%d", s, nSamp); fflush(stdout);
				treeTmp = cpyTree(tree);
				treeTmp->minTime = minTimeIntInf;
				treeTmp->minTimeInt.inf = minTimeIntInf;
				treeTmp->minTimeInt.sup = minTimeIntSup;
				treeTmp->maxTime = maxTimeIntInf+gsl_rng_uniform(rg)*(maxTimeIntSup-maxTimeIntInf);
				ffe = sampleFossilInt(fos, treeTmp->size);
				for(n=0; n<treeTmp->size; n++) {
					if(treeTmp->node[n].child == NOSUCH) {
						int f;
						double maxTmp;
						switch(fos->status[n]) {
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
				tp->param = &param;
				tp->min = min;
				tp->max = max;
				tp->step = step;
				tp->todo = todo;
				tp->nTodo = nTodo;
				tp->logD = logD[s];
				tp->logCond = logCond[s];
				tp->mutex_number = &mutexN;
				tp->cond_number = &condN;
				tp->number = &nT;
				if((ret = pthread_create(&thread, NULL, threadComputeDistribution, (void*) tp)) == 0) {
					int err;
					if((err = pthread_detach(thread)) == 0) {
						nT++; s++;
					} else
						fprintf (stderr, "Error %d while detaching thread: %s\n", err, (char*) strerror(err));
				} else
					fprintf (stderr, "Error %d while creating thread: %s\n", ret, (char*) strerror(ret));
			}
			cont = (nT > 0);
			if(cont)
				pthread_cond_wait(&condN, &mutexN);
			pthread_mutex_unlock(&mutexN);
		}
		d = (TypeDistribution*) malloc(tree->size*sizeof(TypeDistribution));
		for(n=0; n<tree->size; n++) {
			d[n].size = 0;
			d[n].item = NULL;
		}
		for(n=0; n<nTodo; n++) {
			double sumCond, logSumCond, offset;
			offset = logCond[0][n];
			sumCond = 1.;
			for(s=1; s<nSamp; s++)
				if(!(isinf(logCond[s][n]) || isnan(logCond[s][n]))) {
					if(logCond[s][n]>offset) { /*compute max in offset just to avoid numerical precision issues*/
						sumCond *= exp(offset-logCond[s][n]);
						offset = logCond[s][n];
						sumCond++;
					} else
						sumCond += exp(logCond[s][n]-offset);
				}
			logSumCond = log(sumCond)+offset;
			d[todo[n]].size = logD[0][n].size;
			d[todo[n]].item = (TypeDistributionItem*) malloc(d[todo[n]].size*sizeof(TypeDistributionItem));
			for(i=0; i<d[todo[n]].size; i++) {
				double sumI;
				d[todo[n]].item[i].val = logD[0][n].item[i].val;
				offset = logD[0][n].item[i].dens;
				sumI = 1.;
				for(s=1; s<nSamp; s++) {
					if(!(isinf(logD[s][n].item[i].dens) || isnan(logD[s][n].item[i].dens))) {
						if(logD[s][n].item[i].dens>offset) { /*compute max in offset just to avoid numerical precision issues*/
							sumI *= exp(offset-logD[s][n].item[i].dens);
							offset = logD[s][n].item[i].dens;
							sumI++;
						} else
							sumI += exp(logD[s][n].item[i].dens-offset);
					}
				}
				d[todo[n]].item[i].dens = exp(log(sumI)+offset-logSumCond);
			}
			tree->time[todo[n]] = getMedian(d[todo[n]]);
			sprintf(outputDistribution, "%s_%d.csv", outputPrefix, todo[n]);
			if((fo = fopen(outputDistribution, "w"))) {
				if(outDens) {
					deriveDistribution(&(d[todo[n]]));
					fprintDistribution(fo, d[todo[n]]);
				} else
					fprintDistribution(fo, d[todo[n]]);
				fclose(fo);
			}
//			fprintf(stdout, "Median %.2lf\nMean %.2lf\nMode %.2lf\nQuantile Inf %.2lf\nQuantile sup %.2lf\n", getMedian(d[todo[n]]), getMean(d[todo[n]]), getMode(d[todo[n]]), getQuantileInf(d[todo[n]], 0.025), getQuantileSup(d[todo[n]], 0.025));
		}
		fprintf(stdout,"\nplot '%s_%d.csv' using 1:2 with lines title '%s %d' lw 3", outputPrefix, todo[0], "node", todo[0]);
		for(n=1; n<nTodo; n++)
			fprintf(stdout,", '%s_%d.csv' using 1:2 with lines title '%s %d' lw 3", outputPrefix, todo[n], "node", todo[n]);
		fprintf(stdout,"\n");
		if(format != '0') {
			char *tmp, outputFileNameG[SIZE_BUFFER_CHAR], *outputFileName;
			outputFileName = inputFileNameTree;
			for(n=0; n<tree->size; n++) {
				if(tree->node[n].child == NOSUCH) {
							int f;
							double max;
					switch(fos->status[n]) {
						case contempNodeStatus:
							tree->time[n] = tree->maxTime;
						break;
						case unknownNodeStatus:
							//tree->time[n] = NO_TIME;
						break;
						case extinctNodeStatus:
						//case unknownNodeStatus:
							max = NEG_INFTY;
							for(f=fos->fossilInt[n]; f!=NOSUCH; f=fos->fossilIntList[f].prec)
								if(fos->fossilIntList[f].fossilInt.sup>max)
									max = fos->fossilIntList[f].fossilInt.sup; 
							tree->time[n] = max;
						break;
						default:
							fprintf(stderr, "Node %d has no status\n", n);
							return 1;
					}
				}
			}
			if((tmp = strrchr(outputFileName, '.')) != NULL)
				tmp[0] = '\0';
			info.param.tmin = tree->minTime;
			info.param.tmax = tree->maxTime;
			info.param.ratio = 0.55;
			info.param.scaleStep = stepFig;
//			data.color = (TypeRGB) {.red = 1., .green = 1., .blue = 0.};
			data.color = (TypeRGB) {.red = 0.63, .green = 0.32, .blue = 0.18};
			data.radius = 5.;
			data.alpha = 0.5;
			data.fos = fos;
			dataD.color = (TypeRGB) {.red = 1., .green = 0., .blue = 0.};
			dataD.alpha = 0.5;
			dataD.dens = d;
			dataFD.dataFossil = &data;
			dataFD.dataDensity = &dataD;
			add.data = &dataFD;
			add.draw = drawFossilDensity;
			tree->maxTime = getMaximumTime(tree);
//printf("min %.2lf %s\n", tree->minTime, sprintRGBTikz(buffer, dataD.color));
			if(d[tree->root].size>0) {
				int i;
				for(i=0; i<d[tree->root].size && d[tree->root].item[i].dens<0.001; i++)
					;
				tree->minTime = d[tree->root].item[i].val;
			}
			fillUnknownTimesFossilInt(tree->minTime, tree->maxTime, tree, fos);
			switch(format) {
				case '1':
					sprintf(outputFileNameG, "%s_tree.pdf", outputFileName);
					setFunctPDF(&(info.funct));
					data.drawLineDot = drawLineDotCairo;
					dataD.fillPolygon = fillPolygonCairo;
					drawTreeFileGeneric(outputFileNameG, tree, &info, &add);
					break;
				case '2':
					sprintf(outputFileNameG, "%s_tree.ps", outputFileName);
					setFunctPS(&(info.funct));
					data.drawLineDot = drawLineDotCairo;
					dataD.fillPolygon = fillPolygonCairo;
					drawTreeFileGeneric(outputFileNameG, tree, &info, &add);
					break;
				case '3':
					sprintf(outputFileNameG, "%s_tree.png", outputFileName);
					setFunctPNG(&(info.funct));
					data.drawLineDot = drawLineDotCairo;
					dataD.fillPolygon = fillPolygonCairo;
					drawTreeFileGeneric(outputFileNameG, tree, &info, &add);
					break;
				case '4':
					sprintf(outputFileNameG, "%s_tree.svg", outputFileName);
					setFunctSVG(&(info.funct));
					data.drawLineDot = drawLineDotCairo;
					dataD.fillPolygon = fillPolygonCairo;
					drawTreeFileGeneric(outputFileNameG, tree, &info, &add);
					break;
				case '5':
					sprintf(outputFileNameG, "%s_tree_pst.tex", outputFileName);
					setFunctPSTricks(&(info.funct));
					data.drawLineDot = drawLineDotPSTricks;
					dataD.fillPolygon = fillPolygonPSTricks;
					drawTreeFileGeneric(outputFileNameG, tree, &info, &add);
					break;
				case '6':
					sprintf(outputFileNameG, "%s_tree_tikz.tex", outputFileName);
					setFunctTikz(&(info.funct));
					data.drawLineDot = drawLineDotTikz;
					dataD.fillPolygon = fillPolygonTikz;
					drawTreeFileGeneric(outputFileNameG, tree, &info, &add);
					break;
				default:
					;
			}
		}
		if(todo != NULL && todo != tab_todo)
			free((void*)todo);
		free((void*)min);
		free((void*)max);
		if(logD != NULL) {
			for(s=0; s<nSamp; s++)
				if(logD[s] != NULL) {
					for(n=0; n<nTodo; n++)
						if(logD[s][n].item != NULL)
							free((void*)logD[s][n].item);
					free((void*)logD[s]);
				}
			free((void*)logD);
		}
		if(logCond != NULL) {
			for(s=0; s<nSamp; s++)
				if(logCond[s] != NULL)
					free((void*)logCond[s]);
			free((void*)logCond);
		}
		if(d != NULL) {
			for(i=0; i<tree->size; i++)
				if(d[i].item != NULL)
					free((void*)d[i].item);
			free((void*)d);
		}
		gsl_rng_free(rg);
		freeTree(tree);
		freeFossilIntFeature(fos);
	} else {
		fprintf(stderr, "Cannot read %s\n", inputFileNameTree);
		exit(1);
	}
	return 0;
}

void *threadComputeDistribution(void *data) {
	int n, i;
	for(n=0; n<((TypeThreadParameter*)data)->nTodo; n++) {
		for(i=0; i<((TypeThreadParameter*)data)->logD[n].size; i++)
			((TypeThreadParameter*)data)->logD[n].item[i].val = ((TypeThreadParameter*)data)->min[n] + ((double)i)*((TypeThreadParameter*)data)->step;
		fillLogDistribution(&(((TypeThreadParameter*)data)->logD[n]), &(((TypeThreadParameter*)data)->logCond[n]), ((TypeThreadParameter*)data)->todo[n], ((TypeThreadParameter*)data)->tree, ((TypeThreadParameter*)data)->ffe, ((TypeThreadParameter*)data)->param);
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

//#define MAX_SIZE_TMP 50
//#define INC_BUFFER 50
//#define IS_SEP(c) (c == ' ' || c == '\t' || c == '\r' || c == '\n' || c == ';' || c == ',')

//char **readList(FILE *f) {
	//char c, tmp[MAX_SIZE_TMP+1], **list;
	//int size, sizeBuffer;

	//sizeBuffer = INC_BUFFER;
	//list= (char**) malloc(sizeBuffer*sizeof(char*));
	//size = 0;
	//do {
		//c = getc(f);
	//} while(c!=EOF && IS_SEP(c)); 
	//while(c != EOF) {
		//int i;
		//i = 0;
		//while(i<MAX_SIZE_TMP && c !=EOF && !IS_SEP(c)) {
			//tmp[i] = c;
			//c = getc(f);
			//i++;
		//}
		//tmp[i++] = '\0';
		//if(i == MAX_SIZE_TMP) {
			//fprintf(stderr, "Ident too long (%s) ...", tmp);
			//exit(1);
		//}
		//if(i>1) {
			//if(size >= sizeBuffer) {
				//sizeBuffer += INC_BUFFER;
				//list = (char**) realloc((void *) list, sizeBuffer*sizeof(char*));
			//}
			//list[size] = (char*) malloc((strlen(tmp)+1)*sizeof(char));
			//strcpy(list[size], tmp);
			//size++;
		//}
		//while(c!=EOF && IS_SEP(c))
			//c=getc(f);
	//}
	//if(size >= sizeBuffer) {
		//sizeBuffer += INC_BUFFER;
		//list = (char**) realloc((void *) list, sizeBuffer*sizeof(char*));
	//}
	//list[size++] = NULL;
	//return list;
//}
