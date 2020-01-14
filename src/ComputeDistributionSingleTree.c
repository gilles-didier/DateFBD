#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#include <gsl/gsl_rng.h>

#include "Utils.h"
#include "Tree.h"
#include "Fossil.h"
#include "FossilInt.h"
#include "Model.h"
#include "Uncertainty.h"
#include "MinimizeNLOpt.h"
#include "Distribution.h"
#include "DrawTreeCairo.h"
#include "DrawTreePSTricks.h"
#include "DrawTreeTikz.h"
#include "DrawTreeGeneric.h"
#include "DrawFossilInt.h"
#include "DrawDensity.h"
#include "MCMCImportanceSampling.h"

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

#define HELPMESSAGE "\n\nusage: Distribution [options] <input file> [<output file>]\n\nEstimate the diversification rates of the tree contained in the input file.\nThe input file has to be in Newick format with special tags for fossils ages and origin and end of the diversification, \nit returns a text report with the estimates.\n\nOptions are:\n\t-o <options file name>\tload the settings of the optimizer. <options file name> has to be in the format:\n\t\t:SPE [0;1] :EXT [0;1] :FOS [0:1] :TRI 10 :TOL 1.E-7 :ITE 500\n\t-h\tdisplay help\n\n"

//./disp -o -360 -315 -p 0.18 0.17 0.04 -e -s 100 -d 500 Amniota.txt ../data/Cotylosauria100trees.phy ../data/CotylosauriaAges.csv 
//./disx -o -360 -319 -p 1.668294E-01 1.666746E-01 3.939457E-02 -s 2000 -t 8 -d 100 -e Amniota.txt ../data/Cotylosauria100trees.phy ../data/CotylosauriaAges.csv cotySTA3.csv
//./disx -o -360 -319 -p 1.668294E-01 1.666746E-01 3.939457E-02 -s 500 -t 40 -d 0.1 -f 6 -e ../data/Cotylosauria100trees.phy ../data/CotylosauriaAges.csv cotySTA3M10000E

typedef struct DATA_DRAW_FOSSIL_DENSITY {
	TypeDataDrawFossilInt *dataFossil;
	TypeDataDrawDensity *dataDensity;
} TypeDataDrawFossilDensity;


void drawFossilDensity(int n , double x, double y, TypeInfoDrawTreeGeneric *info, void *data) {
	drawFossilInt(n, x, y, info, (void*) ((TypeDataDrawFossilDensity*)data)->dataFossil);
	drawDensity(n, x, y, info, (void*) ((TypeDataDrawFossilDensity*)data)->dataDensity);
}



int main(int argc, char **argv) {	
	char *inputFileNameTree, *inputFileNameFossil = NULL, *outputPrefix = PREFIX, outputDistribution[STRING_SIZE], option[256], format = '1';
	FILE *fi, *fo;
	int i, j, nSamp = 100, nBurn = 1000, nGap = 10, node = NOSUCH, outDens = 1, maxT = 2, *todo, tab_todo[MAX_NODE], nTodo=0;
	double minTimeIntInf = NO_TIME, minTimeIntSup = NO_TIME, maxTimeIntInf = 0., maxTimeIntSup = 0., step = 0.1, treeScaleStep=10., al = 0.75, prop = 0.25, probSpe = 0.33, probExt = 0.33;
	TypeModelParam param = {.birth=1.3, .death = 1., .fossil = 1., .sampling = 1.}, windSize = {.birth=0.1, .death = 0.1, .fossil = 0.02, .sampling = 1.}, init = {.birth=0.5, .death = 0.5, .fossil = 0.1, .sampling = 1.};

	
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
//				bltoabsTime(tree);
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
				//TypeInfoDrawTreeGeneric info;
				//double *timeSave;
				//int n;
				//timeSave = tree->time;
				//tree->time = (double*) malloc(tree->size*sizeof(double));
				//for(n=0; n<tree->size; n++)
					//tree->time[n] = timeSave[n];
				//fillUnknownTimes(tree->minTime, tree->maxTime, tree);
				//setFunctPDF(&(info.funct));
				//drawTreeFileGenericDebug("tree.pdf", tree, &info, NULL);
				//free((void*) tree->time);
				//tree->time = timeSave;
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
		if(option['a']) {
			option['a'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%le", &probSpe) == 1)
				i++;
			else
				exitProg(ErrorArgument, "2 values are expected after -a");
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &probExt) == 1)
				i++;
			else
				exitProg(ErrorArgument, "2 values are expected after -a");
		}
		if(option['f']) {
			option['f'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &prop) == 1)
				i++;
			else
				exitProg(ErrorArgument, "a number is required after option -t");
		}
		if(option['w']) {
			option['w'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%le", &windSize.birth) == 1)
				i++;
			else
				exitProg(ErrorArgument, "3 values are expected after -a");
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &windSize.death) == 1)
				i++;
			else
				exitProg(ErrorArgument, "3 values are expected after -a");
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &windSize.fossil) == 1)
				i++;
			else
				exitProg(ErrorArgument, "3 values are expected after -a");
		}
		if(option['i']) {
			option['i'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%le", &init.birth) == 1)
				i++;
			else
				exitProg(ErrorArgument, "3 values are expected after -a");
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &init.death) == 1)
				i++;
			else
				exitProg(ErrorArgument, "3 values are expected after -a");
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &init.fossil) == 1)
				i++;
			else
				exitProg(ErrorArgument, "3 values are expected after -a");
		}
		if(option['q']) {
			option['q'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &treeScaleStep) == 1)
				i++;
			else
				exitProg(ErrorArgument, "a number is required after option -t");
		}
		if(option['e']) {
			option['e'] = 0;
			outDens = 0;
		}
		if(option['s']) {
			option['s'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%d", &nSamp) == 1)
				i++;
			else
				exitProg(ErrorArgument, "a number is expected after -s");
		}
		if(option['b']) {
			option['b'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%d", &nBurn) == 1)
				i++;
			else
				exitProg(ErrorArgument, "a number is expected after -b");
		}
		if(option['g']) {
			option['g'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%d", &nGap) == 1)
				i++;
			else
				exitProg(ErrorArgument, "a number is expected after -b");
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
		TypeDistribution *d;
		int i, n;
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
printf("Fossil %d %s\n", fos->sizeFossil, inputFileNameFossil);
		if(getMaxFossilIntTime(fos) > 0.)
			negateFossilInt(fos);
		fixStatus(tree, fos);
//fprintTreeFossilInt(stdout, tree, fos);
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
		d = MCMCSamplingMultipleDistSingleTree(todo, nTodo, tree, fos, step, al, nBurn, nGap, nSamp, prop, &windSize, &init, probSpe, probExt);
		for(n=0; n<nTodo; n++) {
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
			info.param.scaleStep = treeScaleStep;
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
